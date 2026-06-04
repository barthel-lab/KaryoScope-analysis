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

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from karyoscope_analysis.core.feature_align import (
    Segment,
    SubScore,
    align_best_orientation,
    align_local,
    classify_overlap,
    feature_jaccard,
    reverse_segments,
)


def idf_weights(reads: Mapping[str, Sequence[Segment]], *, floor: float = 0.1) -> dict[str, float]:
    """Per-feature weight ``max(floor, 1 - document_frequency)`` over the read set.

    A feature present in (almost) every read carries little discriminative information and is
    down-weighted toward ``floor``; a feature seen in few reads keeps weight ~1. This is the
    Engine B answer to repeat-driven chaining: ubiquitous interspersed repeats (LINE/SINE/…)
    stop driving overlaps, so edges must rest on distinctive structure. ``floor`` keeps even
    ubiquitous features slightly informative (and avoids fully disconnecting the graph).
    """
    n = len(reads)
    if n == 0:
        return {}
    document_frequency: Counter[str] = Counter()
    for segments in reads.values():
        for feature in {feat for feat, _ in segments}:
            document_frequency[feature] += 1
    return {feat: max(floor, 1.0 - count / n) for feat, count in document_frequency.items()}


def _weight(weight: Mapping[str, float] | None, feature: str) -> float:
    return weight.get(feature, 1.0) if weight is not None else 1.0


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


def build_overlap_graph(
    reads: Mapping[str, Sequence[Segment]],
    *,
    sub_score: SubScore,
    gap_factor: float,
    match_score: float = 1.0,
    min_overlap_bp: float = 1,
    min_identity: float = 0.8,
    min_jaccard: float = 0.0,
    weight: Mapping[str, float] | None = None,
) -> list[OverlapEdge]:
    """Compute the proper-overlap edges among ``reads`` (``{read_id: segments}``).

    A pair becomes an edge iff its best-orientation alignment is a dovetail or containment,
    has at least ``min_overlap_bp`` of (weighted) overlap, and normalized identity ≥
    ``min_identity``. With ``weight`` (e.g. from :func:`idf_weights`), match rewards and the
    overlap length are scaled per feature, so the overlap must rest on *distinctive* shared
    content — the fix for repeat-driven chaining. The ``min_jaccard`` feature-set prefilter
    prunes obviously-disjoint pairs cheaply (0 = off).
    """
    scorer = _weighted_sub_score(sub_score, weight)
    ids = sorted(reads)
    edges: list[OverlapEdge] = []
    for ai in range(len(ids)):
        a_id = ids[ai]
        a = reads[a_id]
        for bi in range(ai + 1, len(ids)):
            b_id = ids[bi]
            b = reads[b_id]
            if min_jaccard > 0.0 and feature_jaccard(a, b) < min_jaccard:
                continue
            aln = align_best_orientation(a, b, sub_score=scorer, gap_factor=gap_factor)
            kind = classify_overlap(aln, len(a), len(b))
            if kind not in ("dovetail", "containment"):
                continue
            b_used = reverse_segments(b) if aln.reversed_b else b
            overlap_bp = _weighted_overlap(aln.columns, a, b_used, weight)
            if overlap_bp < min_overlap_bp:
                continue
            identity = aln.score / (match_score * overlap_bp) if overlap_bp else 0.0
            if identity < min_identity:
                continue
            edges.append(
                OverlapEdge(a_id, b_id, aln.score, identity, overlap_bp, kind, aln.reversed_b)
            )
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


def cluster_reads(
    reads: Mapping[str, Sequence[Segment]], edges: Sequence[OverlapEdge]
) -> list[Cluster]:
    """Group reads into connected components, oriented to each component's longest read."""
    dsu = _ParityDSU()
    for read_id in reads:
        dsu.add(read_id)
    for edge in edges:
        dsu.union(edge.a, edge.b, 1 if edge.flipped else 0)

    components: dict[str, list[str]] = {}
    for read_id in reads:
        root, _ = dsu.find(read_id)
        components.setdefault(root, []).append(read_id)

    def total_bp(read_id: str) -> int:
        return sum(length for _, length in reads[read_id])

    clusters: list[Cluster] = []
    for root, members in components.items():
        seed = max(members, key=lambda r: (total_bp(r), r))
        _, seed_parity = dsu.find(seed)
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
                orientation_conflict=root in dsu.conflicts,
            )
        )
    clusters.sort(key=lambda c: (-c.size, c.seed))
    return clusters


@dataclass(frozen=True)
class ConsensusPosition:
    """One position of a cluster's seed-anchored consensus."""

    feature: str  # majority feature at this seed position
    length: int  # the seed segment's length (the backbone coordinate)
    support: int  # votes for the majority feature (incl. the seed)
    coverage: int  # total votes at this position (incl. the seed)


@dataclass(frozen=True)
class ClusterConsensus:
    """A cluster's seed-anchored consensus structural sequence."""

    seed: str
    positions: tuple[ConsensusPosition, ...]

    def segments(self) -> list[Segment]:
        """The consensus as a ``(feature, length)`` sequence (the backbone structure)."""
        return [(p.feature, p.length) for p in self.positions]


def _majority(counter: Counter[str], prefer: str) -> str:
    """Most-voted feature; ties broken toward ``prefer`` (the seed), then lexicographically."""
    best = max(counter.values())
    tied = sorted(feat for feat, count in counter.items() if count == best)
    return prefer if prefer in tied else tied[0]


def cluster_consensus(
    reads: Mapping[str, Sequence[Segment]],
    cluster: Cluster,
    *,
    sub_score: SubScore,
    gap_factor: float,
    weight: Mapping[str, float] | None = None,
) -> ClusterConsensus:
    """Seed-anchored consensus: per seed position, the majority feature over members.

    Each member is oriented to the seed and locally aligned to it; every aligned column
    votes its feature at the corresponding seed position (the seed votes for itself). The
    consensus uses the seed's segment lengths as the backbone coordinate (v1). Members linked
    only transitively (not overlapping the seed) contribute no votes — a v1 limitation.
    """
    scorer = _weighted_sub_score(sub_score, weight)
    seed_segments = reads[cluster.seed]
    votes: list[Counter[str]] = [Counter() for _ in seed_segments]
    for position, (feature, _length) in enumerate(seed_segments):
        votes[position][feature] += 1  # the seed votes for itself

    for member in cluster.members:
        if member == cluster.seed:
            continue
        member_segments = (
            reverse_segments(reads[member])
            if cluster.reversed_relative_to_seed[member]
            else list(reads[member])
        )
        aln = align_local(member_segments, seed_segments, sub_score=scorer, gap_factor=gap_factor)
        for member_idx, seed_idx in aln.columns:
            votes[seed_idx][member_segments[member_idx][0]] += 1

    positions = tuple(
        ConsensusPosition(
            feature=_majority(votes[position], prefer=feature),
            length=length,
            support=votes[position][_majority(votes[position], prefer=feature)],
            coverage=sum(votes[position].values()),
        )
        for position, (feature, length) in enumerate(seed_segments)
    )
    return ClusterConsensus(seed=cluster.seed, positions=positions)


def assemble(
    reads: Mapping[str, Sequence[Segment]],
    *,
    sub_score: SubScore,
    gap_factor: float,
    match_score: float = 1.0,
    min_overlap_bp: float = 1,
    min_identity: float = 0.8,
    min_jaccard: float = 0.0,
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
        weight=weight,
    )
    return cluster_reads(reads, edges), edges
