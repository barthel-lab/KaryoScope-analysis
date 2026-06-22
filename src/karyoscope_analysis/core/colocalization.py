"""Alignment-free colocalization measurement (rearrangement detection, Engine A).

See ``docs/audit/rearrangement_detection.md`` for the model and the full Engine A spec.
This module is the **measurement layer only**: for one read's coalesced feature-segment
partition (an ``overlay-annotations`` output), the *minimum bp gap* between every
co-present pair of distinct feature types — the building block for the per-pair
colocalization rate and the experiment-vs-control differential test that live on top of it.

The measurement is:

* in **bp** (physical distance; invariant to annotation granularity), not segment count;
* **orientation-invariant** — pairs are unordered, and a read and its reverse give the
  same gaps (proximity is symmetric); and
* **denoise-aware** — a feature is treated as absent unless it has an interval at least
  ``min_occurrence_bp`` long, so annotation specks don't anchor phantom proximities (the
  gap itself is still measured in true bp, including any filtered specks lying between).

It is pure and non-statistical; reference calibration and the differential test build on it.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from karyoscope_analysis.core.io.bed import Interval, iter_annotation_rows

#: An unordered feature pair, stored canonically (sorted) so ``(A, B) == (B, A)``.
FeaturePair = tuple[str, str]


def feature_pair(a: str, b: str) -> FeaturePair:
    """Canonical (sorted) unordered pair of feature names."""
    return (a, b) if a <= b else (b, a)


def read_span(intervals: Sequence[Interval]) -> int:
    """Total bp spanned by a read's intervals (its length)."""
    if not intervals:
        return 0
    return max(end for _, end, _ in intervals) - min(start for start, _, _ in intervals)


def min_gaps(
    intervals: Sequence[Interval], *, min_occurrence_bp: int = 0
) -> dict[FeaturePair, int]:
    """Minimum bp gap between every co-present pair of distinct feature types on one read.

    Args:
        intervals: a read's feature-segment partition ``[(start, end, feature), ...]``,
            sorted by ``start`` (the C4 invariant; e.g. an ``overlay-annotations`` output).
        min_occurrence_bp: a feature counts as present only via intervals at least this
            long; shorter intervals are ignored when locating it (but their bp still count
            toward distances between other features).

    Returns:
        ``{(A, B): gap}`` over unordered distinct pairs that both occur, where ``gap`` is
        the smallest bp distance between an occurrence of one and the nearest occurrence of
        the other (``0`` when they are adjacent). Adjacency is the distance-0 special case.

    A single left-to-right sweep keeps the most-recent qualifying end per feature; the
    nearest preceding occurrence of every other feature gives that feature's smallest gap to
    the current interval, and the min over all intervals is the global closest approach.
    ``O(segments * distinct-features)``.
    """
    last_end: dict[str, int] = {}
    gaps: dict[FeaturePair, int] = {}
    for start, end, feature in intervals:
        if end - start < min_occurrence_bp:
            continue
        for other, other_end in last_end.items():
            if other == feature:
                continue
            gap = start - other_end
            if gap < 0:
                gap = 0  # defensive; a C4 partition never overlaps
            pair = (feature, other) if feature <= other else (other, feature)
            prev = gaps.get(pair)
            if prev is None or gap < prev:
                gaps[pair] = gap
        last_end[feature] = end  # most recent qualifying occurrence
    return gaps


def iter_read_gaps(
    path: str, *, min_occurrence_bp: int = 0
) -> Iterator[tuple[str, int, dict[FeaturePair, int]]]:
    """Stream ``(seq_id, span, {pair: min_gap})`` per read from an overlay-annotation BED.

    Reads the file once, holding one read's intervals at a time (the same low-memory,
    single-pass approach as the rest of the pipeline). ``.gz`` handled transparently.
    """
    current_seq: str | None = None
    buf: list[Interval] = []
    for seq_id, start, end, feature in iter_annotation_rows(path):
        if seq_id != current_seq:
            if current_seq is not None:
                yield (
                    current_seq,
                    read_span(buf),
                    min_gaps(buf, min_occurrence_bp=min_occurrence_bp),
                )
            current_seq = seq_id
            buf = []
        buf.append((start, end, feature))
    if current_seq is not None:
        yield current_seq, read_span(buf), min_gaps(buf, min_occurrence_bp=min_occurrence_bp)
