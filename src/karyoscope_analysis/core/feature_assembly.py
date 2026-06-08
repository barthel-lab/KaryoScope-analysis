"""Engine B overlap graph + clustering (OLC layout stage).

See ``docs/audit/rearrangement_detection.md`` (Engine B). Builds on the feature aligner
(:mod:`karyoscope_analysis.core.feature_align`):

1. **Overlap graph** — for each read pair (after a feature-Jaccard prefilter) the
   best-orientation local alignment is kept as an edge only if it is a **proper overlap**
   (dovetail or containment) clearing a minimum overlap length and normalized identity.
   Internal-only matches (usually shared repeats) are rejected — this is the anti-chaining
   safeguard that makes connected-components clustering trustworthy.
2. **Clusters** — connected components of that graph (singletons kept). Each edge carries a
   relative orientation; a **parity union-find** assigns every read an orientation relative
   to its cluster seed (the longest read), flagging rare orientation conflicts.

Consensus (seed-anchored majority over the layout) builds on this next.
"""

from __future__ import annotations

import bisect
import multiprocessing as mp
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from itertools import islice

from karyoscope_analysis.core.feature_align import (
    Segment,
    SubScore,
    align_best_orientation,
    align_local,
    classify_overlap,
    feature_jaccard,
    reverse_segments,
)
from karyoscope_analysis.core.io.bed import Interval


def _structural(feature: str) -> str:
    """The structural layer of a (possibly composite) label: ``chr13:aSat`` -> ``aSat``."""
    return feature.split(":", 1)[1] if ":" in feature else feature


def idf_weights(reads: Mapping[str, Sequence[Segment]], *, floor: float = 0.1) -> dict[str, float]:
    """Per-feature weight ``max(floor, 1 - document_frequency)`` over the read set.

    A feature present in (almost) every read carries little discriminative information and is
    down-weighted toward ``floor``; a feature seen in few reads keeps weight ~1. This is the
    Engine B answer to repeat-driven chaining: ubiquitous interspersed repeats (LINE/SINE/…)
    stop driving overlaps, so edges must rest on distinctive structure. ``floor`` keeps even
    ubiquitous features slightly informative (and avoids fully disconnecting the graph).

    Keyed by the **structural layer** so it is meaningful on ``chromosome:structural`` overlays
    (otherwise every chromosome's copy of a feature is counted as a distinct rare feature).
    """
    n = len(reads)
    if n == 0:
        return {}
    document_frequency: Counter[str] = Counter()
    for segments in reads.values():
        for feature in {_structural(feat) for feat, _ in segments}:
            document_frequency[feature] += 1
    return {feat: max(floor, 1.0 - count / n) for feat, count in document_frequency.items()}


def _weight(weight: Mapping[str, float] | None, feature: str) -> float:
    """Weight for a feature, looked up by its **structural layer** (so weighting works on
    ``chromosome:structural`` composite labels). Unknown features default to ``1.0``."""
    return weight.get(_structural(feature), 1.0) if weight is not None else 1.0


def _weighted_sub_score(sub_score: SubScore, weight: Mapping[str, float] | None) -> SubScore:
    """Scale a substitution scorer by the (less distinctive of the) two features' weights."""
    if weight is None:
        return sub_score

    def scored(fa: str, fb: str) -> float:
        return sub_score(fa, fb) * min(_weight(weight, fa), _weight(weight, fb))

    return scored


def _weighted_overlap(
    columns: Sequence[tuple[int, int]],
    a: Sequence[Segment],
    b: Sequence[Segment],
    weight: Mapping[str, float] | None,
) -> float:
    """Distinctive shared content: ``sum over aligned columns of min(weight) * min(len)``."""
    total = 0.0
    for i, j in columns:
        fa, la = a[i]
        fb, lb = b[j]
        total += min(_weight(weight, fa), _weight(weight, fb)) * min(la, lb)
    return total


def _distinctive_overlap(
    columns: Sequence[tuple[int, int]],
    a: Sequence[Segment],
    b: Sequence[Segment],
    weight: Mapping[str, float] | None,
    distinctive_weight: float,
) -> float:
    """Raw bp of matched columns whose feature is *distinctive* (weight ≥ ``distinctive_weight``).

    The anti-arm-chaining criterion: an overlap built only of filler (low-weight ``arm``/``ct``)
    contributes 0 here, so it fails ``min_distinctive_bp`` even though its *weighted* overlap is
    large (a huge arm block times a small weight). Distinctiveness uses the structural-layer weight
    (via :func:`_weight`); with no weighting every feature counts (so this reduces to raw bp).
    """
    total = 0.0
    for i, j in columns:
        fa, la = a[i]
        fb, lb = b[j]
        if min(_weight(weight, fa), _weight(weight, fb)) >= distinctive_weight:
            total += min(la, lb)
    return total


@dataclass(frozen=True)
class OverlapEdge:
    """An accepted proper overlap between two reads."""

    a: str
    b: str
    score: float
    identity: float  # score / (match_score * weighted_overlap), in (.., 1]
    overlap_bp: float  # weighted overlap (distinctive shared bp; raw bp when unweighted)
    kind: str  # "dovetail" | "containment"
    flipped: bool  # B is in the opposite orientation to A


def _all_pairs(n: int):
    """Yield every ``(ai, bi)`` index pair with ``ai < bi`` (the all-vs-all fallback)."""
    for ai in range(n):
        for bi in range(ai + 1, n):
            yield ai, bi


def _block_key(feature: str) -> str | None:
    """Blocking key for a feature: the *specific chromosome* of a composite label, else the label.

    ``chr5:p_arm`` -> ``chr5`` (so all chr5 reads bucket together regardless of which arm —
    robust to structural-label differences); a bare label -> itself; an *ambiguous* chromosome
    layer (``autosome``/``categorized``/…) -> ``None`` (not a key — an edge needs the same
    *specific* chromosome, so ambiguous content can't seed a candidate pair).
    """
    chrom, sep, _struct = feature.partition(":")
    if not sep:
        return feature
    return chrom if chrom.startswith("chr") else None


def _candidate_pairs(
    ids: Sequence[str], reads: Mapping[str, Sequence[Segment]], block_min_bp: float
) -> list[tuple[int, int]]:
    """Read-index pairs sharing a blocking key with ≥ ``block_min_bp`` bp in both.

    Buckets reads by :func:`_block_key` (specific chromosome for composite labels) and returns
    the union of within-bucket pairs, so only plausibly-overlapping reads are aligned — avoiding
    the O(N²) all-vs-all scan and the cross-chromosome blow-up from ambiguous composites.
    """
    index: dict[str, list[int]] = {}
    for idx, read_id in enumerate(ids):
        bp: dict[str, float] = {}
        for feature, length in reads[read_id]:
            key = _block_key(feature)
            if key is not None:
                bp[key] = bp.get(key, 0.0) + length
        for key, total in bp.items():
            if total >= block_min_bp:
                index.setdefault(key, []).append(idx)  # idx ascending -> pairs stay (lo, hi)
    pairs: set[tuple[int, int]] = set()
    for members in index.values():
        for x in range(len(members)):
            for y in range(x + 1, len(members)):
                pairs.add((members[x], members[y]))
    return sorted(pairs)


def _memoized(scorer: SubScore) -> SubScore:
    """Cache the (pure) substitution scorer — it is called once per DP cell but has only a few
    hundred distinct ``(feature, feature)`` argument pairs across an entire run."""
    cache: dict[tuple[str, str], float] = {}

    def scored(fa: str, fb: str) -> float:
        key = (fa, fb)
        value = cache.get(key)
        if value is None:
            value = scorer(fa, fb)
            cache[key] = value
        return value

    return scored


@dataclass(frozen=True)
class _EdgeParams:
    """The numeric thresholds for accepting an overlap edge (bundled to pass to workers)."""

    gap_factor: float
    match_score: float
    min_overlap_bp: float
    min_identity: float
    min_jaccard: float
    min_distinctive_bp: float
    distinctive_weight: float


def _edge_for_pair(
    ai: int,
    bi: int,
    ids: Sequence[str],
    reads: Mapping[str, Sequence[Segment]],
    scorer: SubScore,
    weight: Mapping[str, float] | None,
    params: _EdgeParams,
) -> OverlapEdge | None:
    """Align read ``ai`` vs ``bi`` and return the accepted :class:`OverlapEdge`, or ``None``."""
    a = reads[ids[ai]]
    b = reads[ids[bi]]
    if params.min_jaccard > 0.0 and feature_jaccard(a, b) < params.min_jaccard:
        return None
    aln = align_best_orientation(a, b, sub_score=scorer, gap_factor=params.gap_factor)
    kind = classify_overlap(aln, len(a), len(b))
    if kind not in ("dovetail", "containment"):
        return None
    b_used = reverse_segments(b) if aln.reversed_b else b
    overlap_bp = _weighted_overlap(aln.columns, a, b_used, weight)
    if overlap_bp < params.min_overlap_bp:
        return None
    identity = aln.score / (params.match_score * overlap_bp) if overlap_bp else 0.0
    if identity < params.min_identity:
        return None
    if params.min_distinctive_bp > 0.0 and (
        _distinctive_overlap(aln.columns, a, b_used, weight, params.distinctive_weight)
        < params.min_distinctive_bp
    ):
        return None  # overlap rests only on filler (e.g. shared arm) -> not an edge
    return OverlapEdge(ids[ai], ids[bi], aln.score, identity, overlap_bp, kind, aln.reversed_b)


#: Per-worker state, populated by fork inheritance (copy-on-write — never pickled).
_WORKER: dict = {}


def _edges_for_chunk(chunk: Sequence[tuple[int, int]]) -> list[OverlapEdge]:
    """Worker entry point: edges for a chunk of read-index pairs, using forked-in state."""
    g = _WORKER
    edges = []
    for ai, bi in chunk:
        edge = _edge_for_pair(ai, bi, g["ids"], g["reads"], g["scorer"], g["weight"], g["params"])
        if edge is not None:
            edges.append(edge)
    return edges


def _iter_chunks(items, size: int):
    """Yield ``items`` in lists of at most ``size`` (one parallel task each)."""
    it = iter(items)
    while batch := list(islice(it, size)):
        yield batch


def build_overlap_graph(
    reads: Mapping[str, Sequence[Segment]],
    *,
    sub_score: SubScore,
    gap_factor: float,
    match_score: float = 1.0,
    min_overlap_bp: float = 1,
    min_identity: float = 0.8,
    min_jaccard: float = 0.0,
    min_distinctive_bp: float = 0.0,
    distinctive_weight: float = 0.15,
    block_min_bp: float = 0.0,
    workers: int = 1,
    weight: Mapping[str, float] | None = None,
) -> list[OverlapEdge]:
    """Compute the proper-overlap edges among ``reads`` (``{read_id: segments}``).

    A pair becomes an edge iff its best-orientation alignment is a dovetail or containment,
    has at least ``min_overlap_bp`` of (weighted) overlap, normalized identity ≥
    ``min_identity``, and at least ``min_distinctive_bp`` bp of matched *distinctive* features
    (weight ≥ ``distinctive_weight``). With ``weight`` (e.g. genome-frequency), match rewards
    and the overlap length are scaled per feature, so the overlap rests on *distinctive* shared
    content; ``min_distinctive_bp`` additionally rejects overlaps explained only by filler
    (e.g. a shared chromosome arm) — the fix for arm-chaining. The ``min_jaccard`` feature-set
    prefilter prunes obviously-disjoint pairs cheaply (0 = off).

    ``block_min_bp`` (0 = off) enables a **blocking index** so the alignment is not run on all
    O(N²) pairs: reads are bucketed by the features they carry at least ``block_min_bp`` of, and
    only reads sharing such a "major" feature are compared. Since a ``min_overlap_bp`` edge
    requires substantial shared content (with composite labels, of the same chromosome), reads
    sharing no major feature can't form one — so this scales to whole samples (a heuristic seed,
    like minimizer seeding: an edge resting solely on many sub-``block_min_bp`` features is missed).
    """
    scorer = _memoized(_weighted_sub_score(sub_score, weight))
    ids = sorted(reads)
    params = _EdgeParams(
        gap_factor, match_score, min_overlap_bp, min_identity, min_jaccard,
        min_distinctive_bp, distinctive_weight,
    )
    pair_iter = (
        _candidate_pairs(ids, reads, block_min_bp)
        if block_min_bp > 0.0
        else _all_pairs(len(ids))
    )

    if workers and workers > 1:  # align candidate pairs across processes (fork: shared state)
        pairs = list(pair_iter)
        _WORKER.update(ids=ids, reads=reads, scorer=scorer, weight=weight, params=params)
        try:
            chunksize = max(1, len(pairs) // (workers * 8))
            ctx = mp.get_context("fork")
            edges = []
            with ctx.Pool(workers) as pool:
                for chunk_edges in pool.imap_unordered(
                    _edges_for_chunk, _iter_chunks(pairs, chunksize)
                ):
                    edges.extend(chunk_edges)
        finally:
            _WORKER.clear()
        edges.sort(key=lambda e: (e.a, e.b))  # match the serial path's order (determinism)
        return edges

    edges = []
    for ai, bi in pair_iter:
        edge = _edge_for_pair(ai, bi, ids, reads, scorer, weight, params)
        if edge is not None:
            edges.append(edge)
    return edges


class _ParityDSU:
    """Union-find tracking each element's binary parity (orientation) to its root.

    Union by rank, no path compression (so ``parity`` stays simple); depth is ``O(log n)``.
    """

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}
        self.parity: dict[str, int] = {}  # parity relative to parent
        self.conflicts: set[str] = set()  # roots where an inconsistent edge was seen

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            self.parity[x] = 0

    def find(self, x: str) -> tuple[str, int]:
        parity = 0
        while self.parent[x] != x:
            parity ^= self.parity[x]
            x = self.parent[x]
        return x, parity

    def union(self, x: str, y: str, rel: int) -> bool:
        """Require parity ``rel`` (0 same / 1 flipped) between ``x`` and ``y``.

        Returns ``False`` (and records a conflict) if that contradicts the current state.
        """
        rx, px = self.find(x)
        ry, py = self.find(y)
        if rx == ry:
            if (px ^ py) != rel:
                self.conflicts.add(rx)
                return False
            return True
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
            px, py = py, px
        self.parent[ry] = rx
        self.parity[ry] = rel ^ px ^ py
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        return True


@dataclass(frozen=True)
class Cluster:
    """A connected-component cluster of reads, oriented to a seed (longest read)."""

    members: tuple[str, ...]  # sorted by descending length, then id
    seed: str
    reversed_relative_to_seed: dict[str, bool]
    size: int
    orientation_conflict: bool


def _label_propagation(
    node_ids: Sequence[str], edges: Sequence[OverlapEdge], *, max_iter: int = 100
) -> dict[str, str]:
    """Weighted, deterministic label propagation → ``{read_id: community_label}``.

    Each read starts in its own community and repeatedly adopts the community with the greatest
    total edge weight (overlap bp) among its neighbours, tie-broken by the lexicographically
    smallest label (so it is deterministic and order-independent). Communities never cross
    connected components; an isolated read keeps its own label. This subdivides the dense
    per-chromosome / per-haplotype groups while a sparse "bridge" read (a noisy multi-chromosome
    hub, or a single translocation read linking two groups) joins only its strongest neighbour
    instead of merging everything into one component.
    """
    adj: dict[str, dict[str, float]] = {n: {} for n in node_ids}
    for edge in edges:
        if edge.a == edge.b:
            continue
        adj[edge.a][edge.b] = adj[edge.a].get(edge.b, 0.0) + edge.overlap_bp
        adj[edge.b][edge.a] = adj[edge.b].get(edge.a, 0.0) + edge.overlap_bp

    label = {n: n for n in node_ids}
    order = sorted(node_ids)
    for _ in range(max_iter):
        changed = False
        for node in order:
            neighbours = adj[node]
            if not neighbours:
                continue
            score: dict[str, float] = defaultdict(float)
            for nbr, weight in neighbours.items():
                score[label[nbr]] += weight
            top = max(score.values())
            best = min(lbl for lbl, sc in score.items() if sc == top)
            if best != label[node]:
                label[node] = best
                changed = True
        if not changed:
            break
    return label


def cluster_reads(
    reads: Mapping[str, Sequence[Segment]],
    edges: Sequence[OverlapEdge],
    *,
    communities: bool = False,
) -> list[Cluster]:
    """Group reads into clusters, oriented to each cluster's longest read.

    By default a cluster is a **connected component** of the overlap graph. With
    ``communities=True``, each component is subdivided by :func:`_label_propagation`, so a sparse
    bridge (a noisy multi-chromosome hub, or a lone translocation read) no longer transitively
    merges otherwise-distinct groups into one mega-cluster. Orientation parities come from the
    component-wide union-find either way, so they stay consistent within each sub-community.
    """
    dsu = _ParityDSU()
    for read_id in reads:
        dsu.add(read_id)
    for edge in edges:
        dsu.union(edge.a, edge.b, 1 if edge.flipped else 0)

    groups: dict[str, list[str]] = {}
    if communities:
        label = _label_propagation(list(reads), edges)
        for read_id in reads:
            groups.setdefault(label[read_id], []).append(read_id)
    else:
        for read_id in reads:
            root, _ = dsu.find(read_id)
            groups.setdefault(root, []).append(read_id)

    def total_bp(read_id: str) -> int:
        return sum(length for _, length in reads[read_id])

    clusters: list[Cluster] = []
    for members in groups.values():
        seed = max(members, key=lambda r: (total_bp(r), r))
        _, seed_parity = dsu.find(seed)
        seed_root, _ = dsu.find(seed)
        oriented = {}
        for member in members:
            _, parity = dsu.find(member)
            oriented[member] = bool(parity ^ seed_parity)
        ordered = tuple(sorted(members, key=lambda r: (-total_bp(r), r)))
        clusters.append(
            Cluster(
                members=ordered,
                seed=seed,
                reversed_relative_to_seed=oriented,
                size=len(members),
                orientation_conflict=seed_root in dsu.conflicts,
            )
        )
    clusters.sort(key=lambda c: (-c.size, c.seed))
    return clusters


@dataclass(frozen=True)
class LaidOutRead:
    """A cluster member placed in the consensus coordinate frame (for plotting)."""

    read_id: str
    is_seed: bool
    reversed: bool  # oriented relative to the seed
    segments: tuple[Interval, ...]  # (start, end, feature) in consensus coordinates


@dataclass(frozen=True)
class ConsensusPosition:
    """One interval of a cluster's union-spanning consensus."""

    start: int
    end: int
    feature: str  # majority feature over the reads covering this interval
    support: int  # reads voting the majority feature
    coverage: int  # total reads covering this interval


@dataclass(frozen=True)
class ClusterLayout:
    """A cluster laid out in consensus coordinates: placed reads + the union consensus."""

    seed: str
    width: int  # span of the consensus frame (leftmost edge → rightmost edge), shifted to 0
    placed: tuple[LaidOutRead, ...]
    consensus: tuple[ConsensusPosition, ...]


def _cumulative_bp(segments: Sequence[Segment]) -> list[int]:
    """Prefix sums: ``out[k]`` = bp before segment ``k`` (``out[0] == 0``)."""
    out = [0]
    for _, length in segments:
        out.append(out[-1] + length)
    return out


def _anchor_map(anchors: list[tuple[int, int]]) -> Callable[[int], float]:
    """A piecewise-linear member-bp → consensus-bp map through ``(member, consensus)`` anchors.

    Anchors (the read→seed aligned-column starts) are monotonic in both coordinates. Between
    anchors the map interpolates linearly; outside, it extrapolates with slope 1 (bp-preserving),
    so a read's overhangs extend the consensus frame beyond the seed.
    """
    ms = [m for m, _ in anchors]
    cs = [c for _, c in anchors]

    def mapped(p: int) -> float:
        if p <= ms[0]:
            return cs[0] - (ms[0] - p)
        if p >= ms[-1]:
            return cs[-1] + (p - ms[-1])
        k = bisect.bisect_right(ms, p) - 1
        m0, m1, c0, c1 = ms[k], ms[k + 1], cs[k], cs[k + 1]
        return c0 if m1 == m0 else c0 + (p - m0) * (c1 - c0) / (m1 - m0)

    return mapped


def _union_consensus(placed: Sequence[LaidOutRead], width: int) -> tuple[ConsensusPosition, ...]:
    """Majority-vote consensus over the union grid (all placed-segment boundaries)."""
    breaks = {0, width}
    for read in placed:
        for s, e, _f in read.segments:
            breaks |= {s, e}
    bps = sorted(b for b in breaks if 0 <= b <= width)
    if len(bps) < 2:
        return ()
    votes: list[Counter[str]] = [Counter() for _ in range(len(bps) - 1)]
    for read in placed:
        for s, e, f in read.segments:
            for k in range(bisect.bisect_left(bps, s), bisect.bisect_left(bps, e)):
                votes[k][f] += 1
    out: list[ConsensusPosition] = []
    for k, counter in enumerate(votes):
        if not counter:
            continue
        top = max(counter.values())
        feature = min(f for f, n in counter.items() if n == top)
        nxt = ConsensusPosition(bps[k], bps[k + 1], feature, top, sum(counter.values()))
        if out and out[-1].feature == feature and out[-1].end == nxt.start:  # coalesce
            prev = out[-1]
            out[-1] = ConsensusPosition(
                prev.start, nxt.end, feature,
                max(prev.support, nxt.support), max(prev.coverage, nxt.coverage),
            )
        else:
            out.append(nxt)
    return tuple(out)


def consensus_layout(
    reads: Mapping[str, Sequence[Segment]],
    cluster: Cluster,
    *,
    sub_score: SubScore,
    gap_factor: float,
    weight: Mapping[str, float] | None = None,
) -> ClusterLayout:
    """Lay a cluster out in consensus coordinates: place every read so matched features stack.

    Each member is oriented to the seed and locally aligned to it; its aligned columns become
    anchors for a piecewise-linear coordinate map (:func:`_anchor_map`), so a member segment lands
    on the consensus coordinate of the seed segment it matches (aligned features stack vertically;
    overhangs extend the frame). The whole layout is shifted so the leftmost edge is 0, the
    consensus spans the **union** of all reads, and :func:`_union_consensus` majority-votes the
    feature at each grid interval. A member that doesn't align to the seed is placed at its own
    coordinates (offset 0) — a fallback for transitive links.
    """
    scorer = _memoized(_weighted_sub_score(sub_score, weight))
    seed_segments = reads[cluster.seed]
    seed_cum = _cumulative_bp(seed_segments)

    raw: list[tuple[str, bool, bool, list[Interval]]] = []
    for member in cluster.members:
        reversed_ = cluster.reversed_relative_to_seed[member]
        segments = reverse_segments(reads[member]) if reversed_ else list(reads[member])
        cum = _cumulative_bp(segments)
        if member == cluster.seed:
            placed_segs = [(cum[i], cum[i + 1], segments[i][0]) for i in range(len(segments))]
        else:
            aln = align_local(segments, seed_segments, sub_score=scorer, gap_factor=gap_factor)
            if aln.columns:
                mp = _anchor_map([(cum[mi], seed_cum[sj]) for mi, sj in aln.columns])
                placed_segs = [
                    (round(mp(cum[i])), round(mp(cum[i + 1])), segments[i][0])
                    for i in range(len(segments))
                ]
            else:
                placed_segs = [(cum[i], cum[i + 1], segments[i][0]) for i in range(len(segments))]
        raw.append((member, member == cluster.seed, reversed_, placed_segs))

    starts = [s for *_, segs in raw for s, _e, _f in segs]
    ends = [e for *_, segs in raw for _s, e, _f in segs]
    lo = min(starts) if starts else 0
    width = (max(ends) if ends else 0) - lo
    placed = tuple(
        LaidOutRead(rid, is_seed, rev, tuple((s - lo, e - lo, f) for s, e, f in segs))
        for rid, is_seed, rev, segs in raw
    )
    return ClusterLayout(
        seed=cluster.seed, width=width, placed=placed, consensus=_union_consensus(placed, width)
    )


def assemble(
    reads: Mapping[str, Sequence[Segment]],
    *,
    sub_score: SubScore,
    gap_factor: float,
    match_score: float = 1.0,
    min_overlap_bp: float = 1,
    min_identity: float = 0.8,
    min_jaccard: float = 0.0,
    min_distinctive_bp: float = 0.0,
    distinctive_weight: float = 0.15,
    block_min_bp: float = 0.0,
    workers: int = 1,
    communities: bool = False,
    weight: Mapping[str, float] | None = None,
) -> tuple[list[Cluster], list[OverlapEdge]]:
    """Build the overlap graph and cluster — returns ``(clusters, edges)``."""
    edges = build_overlap_graph(
        reads,
        sub_score=sub_score,
        gap_factor=gap_factor,
        match_score=match_score,
        min_overlap_bp=min_overlap_bp,
        min_identity=min_identity,
        min_jaccard=min_jaccard,
        min_distinctive_bp=min_distinctive_bp,
        distinctive_weight=distinctive_weight,
        block_min_bp=block_min_bp,
        workers=workers,
        weight=weight,
    )
    return cluster_reads(reads, edges, communities=communities), edges
