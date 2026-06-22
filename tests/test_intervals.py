"""Tests for the pure per-sequence interval algebra."""

from __future__ import annotations

import pytest

from karyoscope_analysis.core import intervals as iv


def test_coalesce_merges_same_value_touching():
    assert iv.coalesce([(0, 5, "a"), (5, 10, "a")]) == [(0, 10, "a")]


def test_coalesce_merges_same_value_overlapping():
    assert iv.coalesce([(0, 6, "a"), (4, 10, "a")]) == [(0, 10, "a")]


def test_coalesce_keeps_different_values():
    assert iv.coalesce([(0, 5, "a"), (5, 10, "b")]) == [(0, 5, "a"), (5, 10, "b")]


def test_coalesce_sorts_unsorted_input():
    assert iv.coalesce([(5, 10, "a"), (0, 5, "a")]) == [(0, 10, "a")]


def test_coalesce_empty():
    assert iv.coalesce([]) == []


def test_merge_overlapping():
    assert iv.merge_overlapping([(0, 5), (5, 10), (20, 25)]) == [(0, 10), (20, 25)]
    assert iv.merge_overlapping([(0, 5), (3, 10)]) == [(0, 10)]


def test_total_covered_counts_overlap_once():
    assert iv.total_covered([(0, 5), (3, 10)]) == 10
    assert iv.total_covered([(0, 5), (10, 15)]) == 10


def test_refine_two_tracks():
    region = [(0, 5, "arm"), (5, 10, "centromeric")]
    repeat = [(0, 3, "LINE"), (3, 8, "SINE"), (8, 10, "LTR")]
    assert iv.refine([region, repeat]) == [
        (0, 3, ("arm", "LINE")),
        (3, 5, ("arm", "SINE")),
        (5, 8, ("centromeric", "SINE")),
        (8, 10, ("centromeric", "LTR")),
    ]


def test_refine_single_track_is_identity_of_values():
    track = [(0, 4, "x"), (4, 9, "y")]
    assert iv.refine([track]) == [(0, 4, ("x",)), (4, 9, ("y",))]


def test_refine_empty_tracks():
    assert iv.refine([]) == []


def test_refine_rejects_gap():
    with pytest.raises(ValueError, match="gapless"):
        iv.refine([[(0, 5, "a"), (6, 10, "b")]])


def test_refine_rejects_overlap():
    with pytest.raises(ValueError, match="gapless"):
        iv.refine([[(0, 6, "a"), (5, 10, "b")]])


def test_refine_rejects_different_spans():
    with pytest.raises(ValueError, match="span different ranges"):
        iv.refine([[(0, 10, "a")], [(0, 8, "b")]])


def test_refine_rejects_non_positive_interval():
    with pytest.raises(ValueError, match="non-positive"):
        iv.refine([[(5, 5, "a")]])
