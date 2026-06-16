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
import heapq
import multiprocessing as mp
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import count, islice, pairwise

from karyoscope_analysis.core.feature_align import (
    Alignment,
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
    filler: frozenset[str] | None,
) -> float:
    """Raw bp of matched columns whose feature is *distinctive* (not filler).

    The anti-chaining criterion: an overlap built only of filler contributes 0 here, so it fails
    the ``min_overlap_bp`` size gate even if its weighted overlap is large. When ``filler`` is given, a
    feature is distinctive iff its structural layer is **not** in that set (used to exclude the
    read-set-ubiquitous telomere + arm/ct, which a genome-frequency *weight* can't catch because
    telomere is genome-rare). Otherwise distinctiveness falls back to ``weight ≥ distinctive_weight``
    (and with no weighting every feature counts, i.e. raw bp).
    """
    total = 0.0
    for i, j in columns:
        fa, la = a[i]
        fb, lb = b[j]
        if filler is not None:
            distinctive = _structural(fa) not in filler and _structural(fb) not in filler
        else:
            distinctive = min(_weight(weight, fa), _weight(weight, fb)) >= distinctive_weight
        if distinctive:
            total += min(la, lb)
    return total


#: A matched stretch must total at least this many bp to anchor a transition (a junction), so a
#: sliver can't fabricate the contact between two stretch types that a confident overlap requires.
MIN_STRETCH_BP = 200


def _overlap_transitions(
    columns: Sequence[tuple[int, int]],
    a: Sequence[Segment],
    b: Sequence[Segment],
    min_stretch_bp: float,
) -> int:
    """Number of *transitions* between stretch types within the matched overlap.

    A "stretch type" is a read's composite ``chromosome:feature`` label; a transition is a change
    from one to another along the matched columns (a shared junction — a structural boundary or a
    chromosome breakpoint). Consecutive matched columns of one label are pooled into a stretch and
    counted only if they total ≥ ``min_stretch_bp`` (so a sliver doesn't fabricate a junction). A
    confident overlap has ≥1 transition: it crosses a real junction, not a single uniform shared
    stretch (a coincidental shared repeat).
    """
    stretches: list[str] = []  # confident matched stretches, in alignment order
    cur_label: str | None = None
    cur_bp = 0.0
    for i, j in columns:
        label = a[i][0]
        if label == cur_label:
            cur_bp += min(a[i][1], b[j][1])
            continue
        if cur_label is not None and cur_bp >= min_stretch_bp:
            stretches.append(cur_label)
        cur_label, cur_bp = label, min(a[i][1], b[j][1])
    if cur_label is not None and cur_bp >= min_stretch_bp:
        stretches.append(cur_label)
    return sum(1 for x, y in pairwise(stretches) if x != y)


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
    distinctive_weight: float
    filler: frozenset[str] | None
    require_transition: bool


def _flanks_are_filler(
    aln: Alignment, a: Sequence[Segment], b: Sequence[Segment], filler: frozenset[str] | None
) -> bool:
    """Whether every segment *outside* the aligned span (both reads) is filler.

    An overlap that aligns all the distinctive content but leaves only filler (e.g. a chromosome
    arm) unaligned at the ends is biologically a proper overlap — it just gets classed ``internal``
    because a down-weighted arm failed to align. Distinctive unaligned flanks (a different satellite
    on each read) are *not* filler, so genuine internal-repeat matches stay rejected.
    """
    if filler is None or aln.is_empty:
        return False
    flanks = (
        [a[k][0] for k in range(aln.a_start)]
        + [a[k][0] for k in range(aln.a_end + 1, len(a))]
        + [b[k][0] for k in range(aln.b_start)]
        + [b[k][0] for k in range(aln.b_end + 1, len(b))]
    )
    return all(_structural(f) in filler for f in flanks)


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
    if kind == "none":
        return None
    b_used = reverse_segments(b) if aln.reversed_b else b
    # Proper overlap = dovetail/containment, OR an "internal" overlap whose only *unaligned* flanks
    # are filler (a down-weighted arm that didn't align): biologically a containment, just with the
    # arm clipped. This keeps the anti-chaining gate (distinctive-flanked internal matches — shared
    # repeats between otherwise-unrelated reads — are still rejected).
    if kind not in ("dovetail", "containment") and not _flanks_are_filler(
        aln, a, b_used, params.filler
    ):
        return None
    weighted_bp = _weighted_overlap(aln.columns, a, b_used, weight)
    distinctive_bp = _distinctive_overlap(
        aln.columns, a, b_used, weight, params.distinctive_weight, params.filler
    )
    # The size gate is on *distinctive* (non-filler) matched bp, not the weighted total: a small but
    # real shared structure (e.g. chr20 ITS+TAR1+gSat) shouldn't be rejected just because genome-
    # frequency weighting shrinks its weighted size, and — the anti-chaining half — a filler-only
    # overlap (a shared telomere/arm) scores 0 distinctive bp here, so it never clears the gate.
    if distinctive_bp < params.min_overlap_bp:
        return None
    identity = aln.score / (params.match_score * weighted_bp) if weighted_bp else 0.0
    if identity < params.min_identity:
        return None
    if params.require_transition and _overlap_transitions(aln.columns, a, b_used, MIN_STRETCH_BP) < 1:
        return None  # overlap is a single uniform stretch (no shared junction) -> not confident
    return OverlapEdge(ids[ai], ids[bi], aln.score, identity, weighted_bp, kind, aln.reversed_b)


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
    distinctive_weight: float = 0.15,
    filler_features: frozenset[str] | None = None,
    require_transition: bool = False,
    block_min_bp: float = 0.0,
    workers: int = 1,
    weight: Mapping[str, float] | None = None,
) -> list[OverlapEdge]:
    """Compute the proper-overlap edges among ``reads`` (``{read_id: segments}``).

    A pair becomes an edge iff its best-orientation alignment is a dovetail or containment,
    normalized identity ≥ ``min_identity``, and at least ``min_overlap_bp`` bp of matched
    **distinctive** content (not filler — distinctiveness is "not in ``filler_features``", else
    weight ≥ ``distinctive_weight``). Measuring the size gate on distinctive bp is the anti-chaining
    rule: an overlap explained only by a shared chromosome arm/telomere scores 0 here and is
    rejected. With ``weight`` (e.g. genome-frequency) the match rewards are scaled per feature so the
    overlap *identity* rests on distinctive content. The ``min_jaccard`` feature-set prefilter prunes
    obviously-disjoint pairs cheaply (0 = off).

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
        distinctive_weight, filler_features, require_transition,
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
    """A cluster of reads with a seed (the longest read). Per-read orientation is decided later by
    the layout (:func:`consensus_layout`, off the backbone); here we only flag whether the overlap
    graph's relative orientations were *inconsistent* within the cluster (``orientation_conflict``)."""

    members: tuple[str, ...]  # sorted by descending length, then id
    seed: str
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
        seed_root, _ = dsu.find(seed)  # the parity DSU still flags inconsistent-orientation components
        ordered = tuple(sorted(members, key=lambda r: (-total_bp(r), r)))
        clusters.append(
            Cluster(
                members=ordered,
                seed=seed,
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


def _structural(feature: str) -> str:
    """The structural layer of a ``chromosome:feature`` label (or the whole label if non-composite)."""
    return feature.split(":", 1)[1] if ":" in feature else feature


def _refine_by_concordance(
    placed: dict[str, list[float]],
    oriented_segs: Mapping[str, Sequence[Segment]],
    members: Sequence[str],
    landmark_of: LandmarkOf,
    weight: Mapping[str, float],
    iterations: int = 2,
) -> dict[str, list[float]]:
    """Slide each read to the offset that maximizes feature **concordance** with the cluster
    consensus — the direct "do the features line up?" criterion.

    The backbone anchor places most reads right but pins each on a single junction, so a read whose
    junction it misread sits a feature off. Here every read is re-offset to where its features
    actually overlap matching consensus features. Features are compared by their **backbone landmark**
    (``landmark_of``): a satellite array split across chromosome labels matches as one structural
    token, and a translocation keeps its chromosomes distinct. Each matched bp is scored by ``weight``
    (genome-frequency) so lining up a distinctive ITS/bSat counts far more than a ubiquitous telomere.
    The consensus a read is scored against is built **leaving that read out**, so its own lone block
    can't anchor it where no other read agrees (otherwise it would match its own vote and never move),
    while a block another read also carries still counts. A read only re-offsets when that strictly
    beats staying put. Two passes: the consensus is rebuilt from the refined placement.
    """

    def cell_token(ctr: Counter, exclude: str | None = None) -> str | None:
        best, best_n = None, 0
        for tok, n in ctr.items():
            eff = n - (1 if tok == exclude else 0)
            if eff > best_n or (eff == best_n and best is not None and tok < best):
                best, best_n = tok, eff
        return best if best_n > 0 else None

    placed = {r: list(b) for r, b in placed.items()}
    for _ in range(iterations):
        per_read: dict[str, list[tuple[float, float, str, float]]] = {}
        for r in members:
            bounds, oriented = placed[r], oriented_segs[r]
            segs = []
            for k, (feature, _length) in enumerate(oriented):
                token = landmark_of(feature)
                if token is not None:
                    segs.append((bounds[k], bounds[k + 1], token, weight.get(_structural(feature), 1.0)))
            per_read[r] = segs
        cuts = sorted({round(c) for r in members for s, e, *_ in per_read[r] for c in (s, e)})
        if len(cuts) < 2:
            break
        votes: list[Counter[str]] = [Counter() for _ in range(len(cuts) - 1)]
        cells_of: dict[str, list[tuple[int, str]]] = {r: [] for r in members}
        for r in members:
            for s, e, t, _w in per_read[r]:
                for k in range(bisect.bisect_left(cuts, round(s)), bisect.bisect_left(cuts, round(e))):
                    votes[k][t] += 1
                    cells_of[r].append((k, t))
        full_token = [cell_token(ctr) for ctr in votes]

        changed = False
        for r in members:
            own = {k: t for k, t in cells_of[r]}  # leave-one-out: r's own vote per cell it covers
            # consensus blocks carry their support (how many reads vote that token there), so a read
            # aligns to the *best-supported* copy of a feature — a read mis-anchored ~12 kb off in a
            # second register snaps to the majority's register rather than overlapping a few like-
            # misplaced reads equally (the overlap *length* is the same in either register).
            cons: dict[str, list[tuple[float, float, int]]] = defaultdict(list)
            for k in range(len(votes)):
                tok = cell_token(votes[k], own[k]) if k in own else full_token[k]
                if tok is not None:
                    support = votes[k][tok] - (1 if own.get(k) == tok else 0)
                    block = cons[tok]
                    # Coalesce a whole contiguous same-token block and take its peak support, so the
                    # support gradient *within* one block (a big uniform hub's center vs edges) can't
                    # pull a read to slide inside it — only switching to a separate, better-supported
                    # copy of the feature changes the score.
                    if block and block[-1][1] == cuts[k]:
                        block[-1] = (block[-1][0], cuts[k + 1], max(block[-1][2], support))
                    else:
                        block.append((cuts[k], cuts[k + 1], support))

            def concordance(off: float, cons=cons, segs=per_read[r]) -> float:
                score = 0.0
                for s, e, t, w in segs:
                    a0, a1 = s + off, e + off
                    for cs, ce, support in cons.get(t, ()):
                        overlap = min(a1, ce) - max(a0, cs)
                        if overlap > 0:
                            score += overlap * w * support
                return score

            offsets = {0.0}
            for s, _e, t, _w in per_read[r]:
                offsets.update(cs - s for cs, _ce, _sup in cons.get(t, ()))
            best_off = max(offsets, key=concordance)
            if best_off != 0.0 and concordance(best_off) > concordance(0.0):
                placed[r] = [b + best_off for b in placed[r]]
                changed = True
        if not changed:
            break
    return placed


#: A landmark run must total at least this many bp to count as a real block (not a noise sliver of
#: ambiguous annotation) when **anchoring** a read on the backbone — kept high so a small feature
#: can't hijack the anchor coordinate.
MIN_FEATURE_BLOCK_BP = 500
#: Orientation (which way a read is flipped) only needs the *order* of features, not a stable anchor,
#: so it uses a lower threshold — enough to see a real but small distinctive feature (an ITS ~300 bp,
#: a gSat ~230 bp) and orient the read by it, instead of falling back to a coin-flip best-fit.
ORIENT_MIN_BLOCK_BP = 150
#: Two backbone landmarks within this many bp are a *contact* (a breakpoint or a feature junction
#: to anchor on); farther apart they are arm-separated and not a junction.
ADJACENT_GAP_BP = 2000
#: One side of a chromosome junction must be a contiguous block at least this large — a clear
#: chromosome arm/centromere — for the junction to mark a translocation (vs a shared subtelomere,
#: where both sides are small satellite blocks).
MAJOR_CHROM_BP = 5000
#: The *other* side of a chromosome junction must reach at least this — enough that a real (if short)
#: translocation partner counts, but a satellite sliver mis-assigned across chromosomes does not.
JUNCTION_PARTNER_BP = 1000
#: A clean translocation may carry at most this many chromosomes *beyond* the junction pair as
#: substantial recurring blocks (one extra = a real 3-way fusion). More than this means the reads
#: thread many chromosomes (a noisy centromeric-satellite "hub" cluster) — not a simple
#: translocation, so it lays out on the structural backbone instead.
MAX_OFF_TARGET_CHROMS = 1

#: A landmark function maps a composite ``chromosome:feature`` label to the backbone token it
#: contributes, or ``None`` to skip it.
LandmarkOf = Callable[[str], "str | None"]


def _chromosome_landmark(acrocentrics: frozenset[str]) -> LandmarkOf:
    """Backbone token = the chromosome (``chrN``), acrocentrics collapsed to one ``acrocentric``.

    Used only for a **clean translocation** — a cluster whose reads carry a large, obvious block from
    each of several chromosomes (see :func:`_translocation_chromosomes`). There the chromosomes *are* the
    structure and their breakpoints are the thing to line up.
    """

    def landmark(feature: str) -> str | None:
        chrom = feature.split(":", 1)[0] if ":" in feature else ""
        if not chrom.startswith("chr"):
            return None
        return "acrocentric" if chrom in acrocentrics else chrom

    return landmark


def _structural_landmark(
    structureless: frozenset[str], telomeres: frozenset[str] = frozenset()
) -> LandmarkOf:
    """Backbone token = the structural feature alone, chromosome **stripped** — the default backbone.

    Stripping the chromosome means a satellite array annotated to several chromosome labels (a bSat
    split across ``chr18``/``chr22``/``submetacentric``) is one ``bSat`` landmark, not a scramble, so
    chimeric clusters — lots of small signal from many chromosomes — lay out by their feature
    skeleton (``telomere → TAR1``, ``bSat → ITS``). ``structureless`` arms/ct contribute nothing;
    telomeres and satellites are kept (they are filler for *clustering* but real layout landmarks).

    Telomere subtypes collapse to one ``telomere`` token (like the acrocentric and split-satellite
    collapses): canonical and noncanonical telomere are one contiguous chromosome-terminal region, so
    keeping them distinct splits one ``TAR1 → telomere`` junction into competing ``TAR1 → canonical``
    / ``TAR1 → noncanonical`` junctions that anchor reads to different slots.
    """

    def landmark(feature: str) -> str | None:
        structural = feature.split(":", 1)[1] if ":" in feature else feature
        if structural in structureless:
            return None
        return "telomere" if structural in telomeres else structural

    return landmark


def _translocation_chromosomes(
    reads: Mapping[str, Sequence[Segment]],
    members: Sequence[str],
    acrocentrics: frozenset[str],
    min_block_bp: int,
) -> set[str]:
    """The chromosomes of a **clean translocation**, or an empty set if the cluster isn't one.

    A clean translocation has a **consistent chromosome junction**: two **directly consecutive**
    chromosome blocks on different chromosomes — one large (≥ ``min_block_bp``, a clear arm or
    centromere) and the other ≥ ``JUNCTION_PARTNER_BP`` — and the *same chromosome pair* recurs in at
    least half the reads. "Consecutive" means nothing real lies between them: a ≥ MIN_FEATURE_BLOCK_BP
    run of a third chromosome (e.g. a subtelomere's ``canonical_telomere``, annotated to chr2/chr20,
    sitting between a chr4 arm and an acrocentric centromere) **breaks** the adjacency — that is a
    chr4-end-then-acrocentric structure, not a chr4-acrocentric fusion. That recurring, sized,
    consecutive junction tells a real chr11-chr13 fusion from:

    * a **chimera** whose reads cross many chromosomes but never agree on a pair (no pair reaches half);
    * a **shared subtelomere** (a TAR1 array annotated chr4 then acrocentric): both sides are small
      satellite blocks, so neither reaches ``min_block_bp``;
    * a **mis-assigned satellite sliver** (a bSat split across labels): the sliver < ``JUNCTION_PARTNER_BP``;
    * a **fuzzy junction**: a chromosome's *run* is large only by accumulating small satellite blocks
      (a ``chr4`` that is just a ~900 bp ``TAR1`` in a ``telomere → TAR1`` subtelomere). The size test
      is on each chromosome's **largest single feature block** — a real arm/centromere chunk — not the
      run total, so a satellite-only partner doesn't qualify as a translocation arm and the cluster
      lays out on its structural features (``TAR1``/``cenSat``), whose landmarks are sharp where the
      chromosome boundary is smeared.

    A pair can recur and still be a chimera: cluster reads also threading many *other* substantial
    chromosomes (a centromeric-satellite read whose chromosome layer is scattered slivers, the noisy
    "hub" reads) is not an *obvious, simple* translocation. So if more than
    ``MAX_OFF_TARGET_CHROMS`` chromosomes **beyond** the junction pair have a substantial block in ≥2
    reads, this isn't a clean translocation either — it lays out far better on the structural
    backbone (its distinctive features, e.g. ``TAR1 → ITS``). One extra is allowed: a real 3-way
    fusion (chr11-chr13-chr19) whose third junction falls just under the half-of-reads bar.

    All of those fall back to the structural backbone. Acrocentrics collapse to ``acrocentric`` first.
    """
    pair_reads: Counter = Counter()
    chrom_reads: Counter = Counter()  # reads in which each chromosome has a substantial block
    for r in members:
        # Coalesce into chromosome runs, tracking each run's largest single feature block, so the
        # junction test can require a real arm/centromere chunk (not an accumulated satellite run).
        runs: list[list] = []  # [chromosome, start, end, largest_feature_block_bp]
        cur_c: str | None = None
        cur_f: str | None = None
        cstart = bstart = pos = 0.0
        cur_max = 0.0  # largest feature block seen so far in the current chromosome run
        for feature, length in reads[r]:
            chrom: str | None = feature.split(":", 1)[0] if ":" in feature else ""
            if chrom in acrocentrics:
                chrom = "acrocentric"
            elif not chrom.startswith("chr"):
                chrom = None  # an ambiguous group label / non-composite feature is not a chromosome
            feat = feature.split(":", 1)[1] if ":" in feature else feature
            if chrom != cur_c:
                if cur_c is not None:
                    runs.append([cur_c, cstart, pos, max(cur_max, pos - bstart)])
                cur_c, cstart, bstart, cur_max = chrom, pos, pos, 0.0
            elif feat != cur_f:  # same chromosome, a new feature block begins
                cur_max = max(cur_max, pos - bstart)
                bstart = pos
            cur_f = feat
            pos += length
        if cur_c is not None:
            runs.append([cur_c, cstart, pos, max(cur_max, pos - bstart)])
        # a real (≥ MIN_FEATURE_BLOCK_BP) chromosome run counts; a sliver is absorbed into the gap so
        # it can't break adjacency, but it also can't bridge two chromosomes into a false junction.
        real = [run for run in runs if run[2] - run[1] >= MIN_FEATURE_BLOCK_BP]
        chrom_reads.update({run[0] for run in real})
        seen: set[frozenset] = set()
        for (ca, _sa, ea, big_a), (cb, sb, _eb, big_b) in pairwise(real):
            if (
                ca != cb
                and min(big_a, big_b) >= JUNCTION_PARTNER_BP  # both have a real block (not a sliver)
                and max(big_a, big_b) >= min_block_bp  # one is a clear arm/centromere
                and sb - ea < ADJACENT_GAP_BP
            ):
                seen.add(frozenset({ca, cb}))
        pair_reads.update(seen)
    threshold = max(2, len(members) / 2)
    pairs = [p for p, n in pair_reads.items() if n >= threshold]
    trans = set().union(*pairs) if pairs else set()
    # Reject chimeras: too many off-target chromosomes recurring substantially => not a simple
    # translocation, lay out on the structural backbone instead.
    recurring = {c for c, n in chrom_reads.items() if n >= 2}
    if len(recurring - trans) > MAX_OFF_TARGET_CHROMS:
        return set()
    return trans


def _landmark_sequence(
    segments: Sequence[Segment], landmark_of: LandmarkOf, min_bp: int
) -> list[str]:
    """The *substantial* backbone tokens a read spans, in read order (consecutive dups collapsed).

    A contiguous run of one landmark is kept only if it totals ≥ ``min_bp``, so a sliver of
    mis-assigned chromosome (or a speck of a feature) doesn't read as structure — skipping it
    re-joins the landmarks on either side.
    """
    runs: list[list] = []  # [token, bp] contiguous same-token runs
    for feature, length in segments:
        token = landmark_of(feature)
        if token is None:
            continue
        if runs and runs[-1][0] == token:
            runs[-1][1] += length
        else:
            runs.append([token, length])
    seq: list[str] = []
    for token, bp in runs:
        if bp >= min_bp and (not seq or seq[-1] != token):
            seq.append(token)
    return seq


def _landmark_blocks_at(
    segments: Sequence[Segment], landmark_of: LandmarkOf, min_bp: int
) -> list[tuple[str, float, float]]:
    """Landmark blocks ``(token, start, end)`` in **read coordinates**, ≥ ``min_bp`` each.

    Unlike :func:`_landmark_sequence` this keeps positions, so the gap between consecutive blocks
    (filler between them) is recoverable — two same-token blocks separated by a long arm stay
    separate (an isolated distal copy of a satellite is not merged with the proximal one)."""
    blocks: list[tuple[str, float, float]] = []
    cur: str | None = None
    cstart = pos = 0.0
    for feature, length in segments:
        token = landmark_of(feature)
        if token != cur:
            if cur is not None:
                blocks.append((cur, cstart, pos))
            cur, cstart = token, pos
        pos += length
    if cur is not None:
        blocks.append((cur, cstart, pos))
    return [(t, a, b) for t, a, b in blocks if b - a >= min_bp]


def _landmark_order(sequences: Iterable[Sequence[str]]) -> dict[str, int]:
    """Rank backbone tokens along the cluster's structure from the per-read landmark sequences.

    Reads give *undirected* adjacencies (weighted by count) between consecutive tokens
    (orientation-blind); for a rearrangement (or a structural backbone) these form a path (e.g.
    ``chr11 - chr13 - chr19``, or ``bSat - ITS - TAR1``). Walk it from a loose endpoint, following
    the strongest adjacency, to assign each token a position. The endpoint is the node with the
    least *total* adjacency weight — a middle token (two strong neighbours) has more, so this keeps
    it in the middle even when a few reads skip it and create a spurious shortcut edge (which would
    fool a fewest-neighbours start). Returns ``{token: rank}`` (``{}`` if none); disconnected extras
    are appended deterministically.
    """
    adjacency: dict[str, Counter] = defaultdict(Counter)
    nodes: set[str] = set()
    for seq in sequences:
        nodes.update(seq)
        for a, b in pairwise(seq):
            if a != b:
                adjacency[a][b] += 1
                adjacency[b][a] += 1
    if not nodes:
        return {}
    start = min(nodes, key=lambda n: (sum(adjacency[n].values()), n))  # the loosest endpoint
    order: list[str] = []
    visited: set[str] = set()
    current: str | None = start
    while current is not None and current not in visited:
        order.append(current)
        visited.add(current)
        nxt = [(cnt, n) for n, cnt in adjacency[current].items() if n not in visited]
        current = max(nxt)[1] if nxt else None
    order.extend(sorted(n for n in nodes if n not in visited))  # disconnected tokens
    return {token: i for i, token in enumerate(order)}


def _orient_reads(
    reads: Mapping[str, Sequence[Segment]],
    members: Sequence[str],
    seed: str,
    rank: Mapping[str, int],
    landmark_of: LandmarkOf,
    min_bp: int,
    scorer: SubScore,
    gap_factor: float,
) -> dict[str, bool]:
    """Orient a cluster's reads relative to the seed (``{read: reversed_vs_seed}``).

    The backbone is the orientation signal: given the landmark ``rank`` along the cluster's
    structure (:func:`_landmark_order` — chromosomes for a translocation, distinctive features for a
    single-chromosome cluster), **flip any read whose landmarks run in the descending direction** —
    so every read reads the backbone in one consistent order. A read is oriented by the first vs last
    of its **core** landmarks (those in more than half the reads): a rare distal tail (a
    ``… → telomere``/``mon`` a few reads reach) is not core, so it isn't the first/last landmark
    compared. When that is indecisive because a *core* landmark sits at both ends — a read carrying
    the same satellite twice, e.g. the ``TAR1 → ITS`` the cluster shares plus a second isolated
    ``TAR1`` 169 kb away — the read is oriented by the **run holding its most conserved junction**
    (the adjacent, gap < ``ADJACENT_GAP_BP``, landmark contact most reads share): the gap puts the
    isolated copy in its own run, so the conserved ``TAR1 → ITS`` decides. A read with neither falls
    back to all of its landmarks (so a real but rare partner — a third chromosome in a few reads —
    still orients it), then to whichever orientation aligns better to the seed.
    """
    blocks = {r: _landmark_blocks_at(reads[r], landmark_of, min_bp) for r in members}
    present: Counter = Counter()
    for r in members:
        present.update({t for t, _s, _e in blocks[r]} & set(rank))
    core = {t for t, c in present.items() if c * 2 > len(members)}  # in more than half the reads
    # Conserved junction: how many reads carry each adjacent (contact) landmark pair.
    junction_reads: Counter = Counter()
    for r in members:
        seen = {
            frozenset((ta, tb))
            for (ta, _a, ea), (tb, sb, _b) in pairwise(blocks[r])
            if ta != tb and sb - ea < ADJACENT_GAP_BP and ta in rank and tb in rank
        }
        junction_reads.update(seen)

    def conserved_run_core(rblocks: list[tuple[str, float, float]]) -> list[int]:
        """Core-landmark ranks of the contact-run carrying the read's most-conserved junction."""
        runs: list[list[tuple[str, float, float]]] = []
        for blk in rblocks:
            if runs and blk[1] - runs[-1][-1][2] < ADJACENT_GAP_BP:
                runs[-1].append(blk)
            else:
                runs.append([blk])

        def run_score(run: list[tuple[str, float, float]]) -> int:
            return max(
                (junction_reads[frozenset((a[0], b[0]))]
                 for a, b in pairwise(run) if a[0] != b[0] and a[0] in core and b[0] in core),
                default=0,
            )

        best = max(runs, key=run_score, default=[])
        return [rank[t] for t, _s, _e in best if t in core]

    orient: dict[str, bool] = {}
    undetermined: list[str] = []
    for r in members:
        ranks = [rank[t] for t, _s, _e in blocks[r] if t in core]  # core landmarks across the read
        if len({*ranks}) < 2 or ranks[0] == ranks[-1]:  # indecisive (a core landmark at both ends)
            ranks = conserved_run_core(blocks[r])  # -> the run holding the most-conserved junction
        if len({*ranks}) < 2 or ranks[0] == ranks[-1]:  # still indecisive -> all of the read's landmarks
            ranks = [rank[t] for t, _s, _e in blocks[r] if t in rank]
        if len({*ranks}) >= 2 and ranks[0] != ranks[-1]:  # decisive end-to-end order -> orient by it
            orient[r] = ranks[-1] < ranks[0]
        else:
            undetermined.append(r)

    # Reads with no order signal of their own: orient by best feature alignment to a *determined*
    # read (one the backbone oriented), preferring the one with the most landmarks. Aligning to a
    # well-oriented reference — not a weak seed that is itself a lone symmetric TAR1 — lets the
    # flanking arm break the tie (which end of the TAR1 the arm/telomere sits on). Falls back to the
    # seed when no read was determined.
    ref = max((r for r in members if r in orient), key=lambda r: len(blocks[r]), default=seed)
    ref_oriented = reverse_segments(reads[ref]) if orient.get(ref, False) else list(reads[ref])
    for r in undetermined:
        if r == ref:
            orient[r] = False
            continue
        fwd = align_local(reads[r], ref_oriented, sub_score=scorer, gap_factor=gap_factor)
        rev = align_local(
            reverse_segments(reads[r]), ref_oriented, sub_score=scorer, gap_factor=gap_factor
        )
        orient[r] = rev.score > fwd.score
    if orient.get(seed):  # frame the cluster so the seed is unflipped (a global, relative-preserving flip)
        orient = {r: not flipped for r, flipped in orient.items()}
    return orient


def _place_fixed(
    child_oriented: Sequence[Segment],
    parent_oriented: Sequence[Segment],
    parent_bounds: Sequence[float],
    scorer: SubScore,
    gap_factor: float,
) -> tuple[list[float], float]:
    """Map an already-oriented read's segment boundaries onto consensus coords via an oriented,
    already-placed parent; returns ``(cons_bounds, alignment_score)``. Orientation is fixed by
    :func:`_orient_reads`, so this only aligns forward to find the coordinate offset."""
    aln = align_local(child_oriented, parent_oriented, sub_score=scorer, gap_factor=gap_factor)
    cum = _cumulative_bp(child_oriented)
    if aln.columns:
        mp = _anchor_map([(cum[ui], parent_bounds[pj]) for ui, pj in aln.columns])
        return [mp(x) for x in cum], aln.score
    return [float(x) for x in cum], aln.score


def _token_chromosome(token: str) -> str | None:
    """The chromosome part of a combined ``chromosome:feature`` backbone token (``None`` if none)."""
    return token.split(":", 1)[0] if ":" in token else None


def _substantial_blocks(
    oriented: Sequence[Segment],
    cum: Sequence[float],
    rank: Mapping[str, int],
    landmark_of: LandmarkOf,
    min_bp: int,
) -> list[tuple[str, float, float]]:
    """Contiguous landmark blocks ``(token, start_cum, end_cum)`` totalling ≥ ``min_bp``.

    Same filter as :func:`_landmark_sequence`, so a speck of a landmark (a 36 bp ``dhor`` that is a
    real block in *another* read, hence in ``rank``) can't act as an anchor here.
    """
    blocks: list[tuple[str, float, float]] = []
    i, n = 0, len(oriented)
    while i < n:
        token = landmark_of(oriented[i][0])
        if token not in rank:
            i += 1
            continue
        bp, j = 0.0, i
        while j < n and landmark_of(oriented[j][0]) == token:
            bp += oriented[j][1]
            j += 1
        if bp >= min_bp:
            blocks.append((token, cum[i], cum[j]))
        i = j
    return blocks


def _junction_placement(
    oriented_segs: Mapping[str, Sequence[Segment]],
    members: Sequence[str],
    rank: Mapping[str, int],
    landmark_of: LandmarkOf,
) -> dict[str, list[float]]:
    """Place reads by anchoring on a **conserved backbone boundary** (chromosomes, or features).

    Landmarks are laid out left to right in ``rank`` order, each given a slot as wide as its longest
    observed block. A read is pinned on the **most cluster-conserved boundary** it carries:

    * a *junction* — a contact between two consecutive landmark blocks (gap < ADJACENT_GAP_BP) —
      preferred, ties broken toward a chromosome change; the right landmark's start pins to its slot;
    * else a *breakpoint* — a landmark block's edge that abuts filler (≥ ADJACENT_GAP_BP of arm/ct
      beyond it), e.g. a ``bSat → q_arm`` boundary; the block's filler-facing edge pins to its slot
      edge, so reads line up on the breakpoint and a variable-length feature extends away from it
      rather than drifting on its far edge;
    * else the read's lowest-rank landmark start.

    Anchoring on the *shared* boundary (rather than how much of a landmark a read captured, or a
    read-specific internal boundary) is what lines every read up.
    """
    max_len: dict[str, float] = defaultdict(float)
    for r in members:
        per: dict[str, float] = defaultdict(float)
        for feature, length in oriented_segs[r]:
            token = landmark_of(feature)
            if token in rank:
                per[token] += length
        for token, total in per.items():
            max_len[token] = max(max_len[token], total)
    slot_start: dict[str, float] = {}
    cursor = 0.0
    for token in sorted(rank, key=lambda c: rank[c]):
        slot_start[token] = cursor
        cursor += max_len[token]
    slot_end = {t: slot_start[t] + max_len[t] for t in slot_start}

    read_data: dict[str, tuple[list[float], list[tuple[str, float, float]], float]] = {}
    for r in members:
        cum = _cumulative_bp(oriented_segs[r])
        blocks = _substantial_blocks(oriented_segs[r], cum, rank, landmark_of, MIN_FEATURE_BLOCK_BP)
        read_data[r] = (cum, blocks, cum[-1] if cum else 0.0)
    junction_freq: Counter = Counter()  # reads carrying each landmark-pair junction
    breakpoint_freq: Counter = Counter()  # reads carrying each (token, side) landmark->filler edge
    for _cum, blocks, total in read_data.values():
        for (ta, _sa, ea), (tb, sb, _eb) in pairwise(blocks):
            if ta != tb and sb - ea < ADJACENT_GAP_BP:
                junction_freq[frozenset({ta, tb})] += 1
        for k, (token, start, end) in enumerate(blocks):
            left_end = blocks[k - 1][2] if k > 0 else 0.0
            right_start = blocks[k + 1][1] if k < len(blocks) - 1 else total
            if start - left_end >= ADJACENT_GAP_BP:
                breakpoint_freq[(token, "L")] += 1
            if right_start - end >= ADJACENT_GAP_BP:
                breakpoint_freq[(token, "R")] += 1

    placed: dict[str, list[float]] = {}
    for r in members:
        cum, blocks, total = read_data[r]
        junction = None  # (freq, is_chromosome_change, block_start, slot)
        breakpoint_ = None  # (freq, read_coord, slot_coord)
        lowest: tuple[int, float, str] | None = None
        for k, (token, start, end) in enumerate(blocks):
            if lowest is None or rank[token] < lowest[0]:
                lowest = (rank[token], start, token)
            if k > 0 and blocks[k - 1][0] != token and start - blocks[k - 1][2] < ADJACENT_GAP_BP:
                prev = blocks[k - 1][0]
                ca, cb = _token_chromosome(prev), _token_chromosome(token)
                key = (junction_freq[frozenset({prev, token})], int(bool(ca and cb and ca != cb)))
                if junction is None or key > junction[:2]:
                    junction = (*key, start, slot_start[token])
            left_end = blocks[k - 1][2] if k > 0 else 0.0
            right_start = blocks[k + 1][1] if k < len(blocks) - 1 else total
            if start - left_end >= ADJACENT_GAP_BP:
                f = breakpoint_freq[(token, "L")]
                if breakpoint_ is None or f > breakpoint_[0]:
                    breakpoint_ = (f, start, slot_start[token])
            if right_start - end >= ADJACENT_GAP_BP:
                f = breakpoint_freq[(token, "R")]
                if breakpoint_ is None or f > breakpoint_[0]:
                    breakpoint_ = (f, end, slot_end[token])
        if junction is not None:
            anchor_read, anchor_cons = junction[2], junction[3]
        elif breakpoint_ is not None:
            anchor_read, anchor_cons = breakpoint_[1], breakpoint_[2]
        elif lowest is not None:
            anchor_read, anchor_cons = lowest[1], slot_start[lowest[2]]
        else:
            anchor_read = anchor_cons = 0.0
        offset = anchor_cons - anchor_read
        placed[r] = [x + offset for x in cum]
    return placed


def consensus_layout(
    reads: Mapping[str, Sequence[Segment]],
    cluster: Cluster,
    *,
    neighbors: Mapping[str, Sequence[str]] | None = None,
    sub_score: SubScore,
    gap_factor: float,
    filler: frozenset[str] | None = None,
    structureless: frozenset[str] | None = None,
    acrocentric_chromosomes: frozenset[str] | None = None,
    weight: Mapping[str, float] | None = None,
) -> ClusterLayout:
    """Lay a cluster out in consensus coordinates: place every read so matched features stack.

    Orientation is decided first (:func:`_orient_reads`), then reads are placed along the cluster's
    **backbone** — the sequence of :func:`_combined_landmark` tokens, the composite
    ``chromosome:feature`` using *both* the chromosome and the structural feature (acrocentrics
    collapsed; ``structureless`` arms/ct dropped; telomeres + satellites kept). So the layout
    anchors on chromosome breakpoints **and** structural junctions (a ``telomere → TAR1`` subtelomere
    boundary, a ``bSat → ITS`` contact).

    With ≥ 2 backbone landmarks the reads are placed by :func:`_junction_placement` — anchored on
    their backbone **junctions**, so every shared breakpoint/junction lines up and the backbone
    reads in one consistent order (robust where alignment isn't: landmark identity, not a ~uniform
    shared satellite, fixes the coordinate). A cluster with < 2 landmarks falls back to a maximum-
    spanning-tree walk over the overlap graph (``neighbors``): from the seed each read is attached
    via its *strongest* overlap to an already-placed read (Prim's), in its fixed orientation, mapped
    via a piecewise-linear :func:`_anchor_map`. The layout is then shifted so the leftmost edge is 0,
    the consensus spans the **union** of all reads, and :func:`_union_consensus` majority-votes the
    feature at each grid interval.
    """
    scorer = _memoized(_weighted_sub_score(sub_score, weight))
    members = list(cluster.members)
    member_set = set(members)
    seed = cluster.seed
    filler = filler if filler is not None else frozenset()
    # The layout backbone excludes only structureless features (arms/ct/...); unlike the clustering
    # filler it *keeps* telomeres + satellites as landmarks. Defaults to filler when not given.
    structureless = structureless if structureless is not None else filler
    # Telomeres are the filler features kept as landmarks (filler minus the structureless arms/ct);
    # collapse their subtypes (canonical/noncanonical) to one ``telomere`` backbone token.
    telomeres = filler - structureless
    if neighbors is None:  # star: every member placed directly against the seed
        neighbors = {seed: [m for m in members if m != seed]}

    # Pick the backbone. A *clean translocation* (≥2 chromosomes each with a large, obvious block)
    # lays out on its chromosomes — their breakpoints are the structure. Every other cluster (one
    # chromosome, or a chimera with small signal scattered over many chromosomes) lays out on its
    # structural features, chromosome stripped — so a satellite array split across chromosome labels
    # is one landmark, not a scramble. This is the principled version of the chromosome-vs-structural
    # split: the chromosome backbone is an opt-in for the easy case, structural is the default.
    acro = acrocentric_chromosomes or frozenset()
    if len(_translocation_chromosomes(reads, members, acro, MAJOR_CHROM_BP)) >= 2:
        landmark_of = _chromosome_landmark(acro)
        orient_min = MIN_FEATURE_BLOCK_BP  # chromosome blocks are large; ignore small slivers
    else:
        landmark_of = _structural_landmark(structureless, telomeres)
        orient_min = ORIENT_MIN_BLOCK_BP  # a small ITS/gSat is enough to orient a structural read

    sequences = {
        r: _landmark_sequence(reads[r], landmark_of, MIN_FEATURE_BLOCK_BP) for r in members
    }
    rank = _landmark_order(sequences.values())
    # Orient the rank in the seed's reading direction (the seed is the unflipped reference), so the
    # rank increases left-to-right exactly as reads are placed — otherwise "lowest-rank landmark"
    # would point the wrong way for a read missing the proximal landmark.
    seed_ranks = [rank[t] for t in sequences[seed] if t in rank]
    if len(seed_ranks) >= 2 and seed_ranks[-1] < seed_ranks[0]:
        top = max(rank.values())
        rank = {t: top - v for t, v in rank.items()}

    # Orientation reads off a *finer* backbone than placement: a small but real distinctive feature
    # (a ~400 bp ITS, a ~230 bp gSat next to a TAR1) is enough to tell which way a structural read
    # runs, even though it is too small to anchor on. Without it, reads carrying one big landmark
    # flip on a coin-toss best-fit.
    orient_seqs = {r: _landmark_sequence(reads[r], landmark_of, orient_min) for r in members}
    orient_rank = _landmark_order(orient_seqs.values())
    seed_oranks = [orient_rank[t] for t in orient_seqs[seed] if t in orient_rank]
    if len(seed_oranks) >= 2 and seed_oranks[-1] < seed_oranks[0]:
        top = max(orient_rank.values())
        orient_rank = {t: top - v for t, v in orient_rank.items()}

    orient = _orient_reads(
        reads, members, seed, orient_rank, landmark_of, orient_min, scorer, gap_factor
    )
    oriented_segs = {
        r: (reverse_segments(reads[r]) if orient[r] else list(reads[r])) for r in members
    }

    if rank:
        # Anchor on the backbone: a junction between landmarks (≥2), or — for a single-distinctive-
        # feature cluster (e.g. TAR1 only) — that one landmark's block, so the distinctive feature
        # lines up instead of the layout drifting on a big variable telomere/arm (the MST failure).
        placed: dict[str, list[float]] = _junction_placement(
            oriented_segs, members, rank, landmark_of
        )
    else:
        # No backbone landmark at all: Prim's — commit the highest-scoring overlap next.
        placed = {seed: [float(x) for x in _cumulative_bp(oriented_segs[seed])]}
        tie = count()
        heap: list[tuple[float, int, str, list[float]]] = []

        def offer(parent: str) -> None:
            p_oriented, p_bounds = oriented_segs[parent], placed[parent]
            for child in neighbors.get(parent, ()):
                if child not in member_set or child in placed:
                    continue
                bounds, score = _place_fixed(
                    oriented_segs[child], p_oriented, p_bounds, scorer, gap_factor
                )
                heapq.heappush(heap, (-score, next(tie), child, bounds))

        offer(seed)
        while heap:
            _neg, _t, child, bounds = heapq.heappop(heap)
            if child in placed:  # stale entry (already attached via a stronger overlap)
                continue
            placed[child] = bounds
            offer(child)
        for member in members:  # unreachable via edges: place at their own coordinates
            if member not in placed:
                placed[member] = [float(x) for x in _cumulative_bp(oriented_segs[member])]

    # Refine: slide each read to maximize feature concordance with the consensus, fixing reads the
    # backbone anchor pinned on the wrong junction (the alignment, done directly rather than by proxy).
    if len(members) >= 2:
        placed = _refine_by_concordance(
            placed, oriented_segs, members, landmark_of, weight or {}
        )

    raw: list[tuple[str, bool, bool, list[Interval]]] = []
    for rid in cluster.members:
        oriented, bounds = oriented_segs[rid], placed[rid]
        segs = [(round(bounds[k]), round(bounds[k + 1]), oriented[k][0]) for k in range(len(oriented))]
        raw.append((rid, rid == seed, orient[rid], segs))
    # order reads top-to-bottom by where they start in the consensus (leftmost first), id-tiebroken
    raw.sort(key=lambda row: (min((s for s, _e, _f in row[3]), default=0), row[0]))

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
    distinctive_weight: float = 0.15,
    filler_features: frozenset[str] | None = None,
    require_transition: bool = False,
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
        distinctive_weight=distinctive_weight,
        filler_features=filler_features,
        require_transition=require_transition,
        block_min_bp=block_min_bp,
        workers=workers,
        weight=weight,
    )
    return cluster_reads(reads, edges, communities=communities), edges
