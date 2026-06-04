"""Feature-sequence local alignment for Engine B (OLC clustering).

See ``docs/audit/rearrangement_detection.md`` (Engine B). A read is represented as an
ordered ``(feature, length)`` segment sequence (the coalesced overlay output). This module
aligns two such sequences with a Smith-Waterman **local** alignment:

* substitution scores from a ``sub_score(a, b)`` callable -- tiered exact > sibling >
  unrelated, with ``novel`` neutral (built from the hierarchy by
  :func:`hierarchy_substitution`, but the aligner itself is hierarchy-agnostic);
* match reward weighted by ``min(len_a, len_b)`` -- bp-weighted, so big shared blocks
  dominate and specks contribute ~nothing, and length-change-lenient (the length difference
  is neither rewarded nor penalized);
* a linear per-bp gap -- skipping a segment costs ``gap_factor * its length``, so skipping a
  large structural block is penalized in proportion;
* the better of the two orientations of B (a read and its reverse get the same score).

The result classifies the overlap as dovetail / containment / internal -- the input to the
proper-overlap edge rule of the overlap graph (built in a later module). Pure; no I/O.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, replace

from karyoscope_analysis.core.feature_vocab import NOVEL, FeatureHierarchy
from karyoscope_analysis.core.io.bed import Interval

#: A read segment: ``(feature, length_bp)``.
Segment = tuple[str, int]
#: Substitution scorer: ``(feature_a, feature_b) -> base score`` (per matched bp).
SubScore = Callable[[str, str], float]


def to_segments(intervals: Sequence[Interval]) -> list[Segment]:
    """Convert overlay intervals ``[(start, end, feature)]`` to ``[(feature, length)]``."""
    return [(feature, end - start) for start, end, feature in intervals]


def reverse_segments(segments: Sequence[Segment]) -> list[Segment]:
    """The read in its reverse orientation (segment order reversed; labels unchanged)."""
    return list(reversed(segments))


@dataclass(frozen=True)
class Alignment:
    """A local alignment of read A against read B (B possibly reversed)."""

    score: float
    columns: tuple[tuple[int, int], ...]  # aligned (i, j) segment-index pairs, A-order
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    reversed_b: bool  # whether B was reversed; spans/columns are in the as-aligned frame

    @property
    def is_empty(self) -> bool:
        return not self.columns


_EMPTY = Alignment(0.0, (), -1, -1, -1, -1, False)


def align_local(
    a: Sequence[Segment], b: Sequence[Segment], *, sub_score: SubScore, gap_factor: float
) -> Alignment:
    """Smith-Waterman local alignment of two feature-segment sequences (B as given).

    ``sub_score(fa, fb)`` is the per-bp base score; a matched column scores
    ``sub_score * min(len_a, len_b)``. A gap (skipping a segment of length L) costs
    ``gap_factor * L``. Returns the highest-scoring local alignment (empty if none is
    positive).
    """
    n, m = len(a), len(b)
    # h[i][j] = best local score ending at A[i-1], B[j-1]; ptr: 0 stop, 1 diag, 2 up, 3 left.
    h = [[0.0] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    best = 0.0
    best_i = best_j = 0
    for i in range(1, n + 1):
        fa, la = a[i - 1]
        for j in range(1, m + 1):
            fb, lb = b[j - 1]
            diag = h[i - 1][j - 1] + sub_score(fa, fb) * min(la, lb)
            up = h[i - 1][j] - gap_factor * la
            left = h[i][j - 1] - gap_factor * lb
            cell = 0.0
            move = 0
            if diag > cell:
                cell, move = diag, 1
            if up > cell:
                cell, move = up, 2
            if left > cell:
                cell, move = left, 3
            h[i][j] = cell
            ptr[i][j] = move
            if cell > best:
                best, best_i, best_j = cell, i, j

    columns: list[tuple[int, int]] = []
    i, j = best_i, best_j
    while i > 0 and j > 0 and ptr[i][j] != 0:
        move = ptr[i][j]
        if move == 1:
            columns.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif move == 2:
            i -= 1
        else:
            j -= 1
    if not columns:
        return _EMPTY
    columns.reverse()
    return Alignment(
        score=best,
        columns=tuple(columns),
        a_start=columns[0][0],
        a_end=columns[-1][0],
        b_start=columns[0][1],
        b_end=columns[-1][1],
        reversed_b=False,
    )


def align_best_orientation(
    a: Sequence[Segment], b: Sequence[Segment], *, sub_score: SubScore, gap_factor: float
) -> Alignment:
    """Align A against B and against reversed B; return the higher-scoring alignment."""
    forward = align_local(a, b, sub_score=sub_score, gap_factor=gap_factor)
    reverse = align_local(a, reverse_segments(b), sub_score=sub_score, gap_factor=gap_factor)
    if reverse.score > forward.score:
        return replace(reverse, reversed_b=True)
    return forward


def classify_overlap(aln: Alignment, n_a: int, n_b: int) -> str:
    """Classify an alignment as ``"containment"``, ``"dovetail"``, ``"internal"``, or ``"none"``.

    ``n_b`` is the length of B in the frame the alignment was computed (reversal is
    orientation-symmetric for this classification). A *proper* overlap is a dovetail or a
    containment; an internal-only match (typically a shared repeat) is rejected by the graph.
    """
    if aln.is_empty:
        return "none"
    a_full = aln.a_start == 0 and aln.a_end == n_a - 1
    b_full = aln.b_start == 0 and aln.b_end == n_b - 1
    if a_full or b_full:
        return "containment"
    a_suffix = aln.a_end == n_a - 1
    a_prefix = aln.a_start == 0
    b_suffix = aln.b_end == n_b - 1
    b_prefix = aln.b_start == 0
    if (a_suffix and b_prefix) or (a_prefix and b_suffix):
        return "dovetail"
    return "internal"


def is_proper_overlap(aln: Alignment, n_a: int, n_b: int) -> bool:
    """Whether the alignment is a proper (dovetail or containment) overlap."""
    return classify_overlap(aln, n_a, n_b) in ("containment", "dovetail")


def hierarchy_substitution(
    hierarchy: FeatureHierarchy,
    *,
    match: float = 1.0,
    partial: float = 0.5,
    mismatch: float = -1.0,
) -> SubScore:
    """Build a tiered per-bp substitution scorer from the feature hierarchy.

    Exact match -> ``match``; ``novel`` vs anything -> ``0`` (neutral); two distinct features
    in the same coarse group (satellite / arm / ct / telomere-type) -> ``partial``; otherwise
    -> ``mismatch``.
    """
    groups: list[frozenset[str]] = [
        hierarchy.satellite_features,
        hierarchy.arm_features,
        hierarchy.ct_features,
        hierarchy.canonical_telomere | hierarchy.noncanonical_telomere | hierarchy.its_tar1,
    ]
    groups = [g for g in groups if g]

    def score(fa: str, fb: str) -> float:
        if fa == fb:
            return match
        if fa == NOVEL or fb == NOVEL:
            return 0.0
        for group in groups:
            if fa in group and fb in group:
                return partial
        return mismatch

    return score


def feature_jaccard(a: Iterable[Segment], b: Iterable[Segment]) -> float:
    """Jaccard overlap of the two reads' feature *sets* — a cheap O(N^2) prefilter."""
    fa = {feat for feat, _ in a}
    fb = {feat for feat, _ in b}
    if not fa or not fb:
        return 0.0
    return len(fa & fb) / len(fa | fb)
