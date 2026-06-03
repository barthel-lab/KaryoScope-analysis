"""Pure per-sequence interval algebra for annotation overlay/merge.

All functions operate on half-open ``[start, end)`` integer intervals on a single
sequence (the per-``seq_id`` coordinate frame). No I/O and no third-party
dependencies — these are the workhorses behind ``overlay-annotations`` (decisions
M2/M5: a single sweep-line implementation, no pyranges).

A valued interval is a ``(start, end, value)`` tuple; an unvalued one is
``(start, end)``. ``value`` is typically a feature name (``str``).
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import pairwise
from typing import TypeVar

T = TypeVar("T")


def coalesce(intervals: Sequence[tuple[int, int, T]]) -> list[tuple[int, int, T]]:
    """Merge adjacent/overlapping intervals that share the same value.

    Input need not be sorted. Touching (``prev_end == start``) or overlapping
    intervals with an equal value are fused; differing values are left split.
    Returns a new list sorted by start.
    """
    items = sorted(intervals, key=lambda iv: (iv[0], iv[1]))
    out: list[tuple[int, int, T]] = []
    for start, end, value in items:
        if out and out[-1][2] == value and start <= out[-1][1]:
            ps, pe, pv = out[-1]
            out[-1] = (ps, max(pe, end), pv)
        else:
            out.append((start, end, value))
    return out


def merge_overlapping(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/touching intervals (ignoring value) into covered ranges.

    Useful for coverage/length computations where overlaps must be counted once
    (e.g. translocation per-chromosome coverage).
    """
    out: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def total_covered(intervals: Sequence[tuple[int, int]]) -> int:
    """Total bp covered by ``intervals``, counting overlapping regions once."""
    return sum(end - start for start, end in merge_overlapping(intervals))


def refine(
    tracks: Sequence[Sequence[tuple[int, int, T]]],
) -> list[tuple[int, int, tuple[T, ...]]]:
    """Common refinement of several interval tracks sharing one coordinate span.

    Each track must be a gapless, non-overlapping tiling of the **same** ``[lo, hi)``
    span (the C4 invariant, validated upstream by the BED reader). The output is
    split at the union of all track boundaries; each segment carries a tuple of
    values, one per track, in the input track order.

    Raises:
        ValueError: if a track is empty, has a gap/overlap, contains a non-positive
            interval, or spans a different range than the others.
    """
    if not tracks:
        return []

    sorted_tracks: list[list[tuple[int, int, T]]] = []
    spans: list[tuple[int, int]] = []
    for i, track in enumerate(tracks):
        items = sorted(track, key=lambda iv: iv[0])
        if not items:
            raise ValueError(f"track {i} is empty")
        for start, end, _ in items:
            if end <= start:
                raise ValueError(f"track {i} has a non-positive interval [{start}, {end})")
        for (_, end, _), (next_start, _, _) in pairwise(items):
            if end != next_start:
                raise ValueError(
                    f"track {i} is not a gapless tiling: gap/overlap at {end} -> {next_start}"
                )
        sorted_tracks.append(items)
        spans.append((items[0][0], items[-1][1]))

    lo, hi = spans[0]
    if any(span != (lo, hi) for span in spans):
        raise ValueError(f"tracks span different ranges: {spans}")

    boundaries = sorted({b for items in sorted_tracks for iv in items for b in (iv[0], iv[1])})
    pointers = [0] * len(sorted_tracks)
    out: list[tuple[int, int, tuple[T, ...]]] = []
    for a, b in pairwise(boundaries):
        values: list[T] = []
        for ti, items in enumerate(sorted_tracks):
            while items[pointers[ti]][1] <= a:
                pointers[ti] += 1
            values.append(items[pointers[ti]][2])
        out.append((a, b, tuple(values)))
    return out
