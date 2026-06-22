"""Tests for the streaming internals of build-feature-matrix.

The CLI builds the matrix by walking all featureset BEDs in lockstep
(:func:`iter_aligned_groups` -> :func:`build_feature_matrix_streaming`), one sequence at
a time. These cover that the streaming path equals the in-memory path, the aligned-group
reader yields whole sequences, and order/coverage disagreements are errors.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.core import build_feature_matrix as core
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import (
    iter_aligned_groups,
    iter_annotation_rows,
    read_annotation_bed,
)

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def test_streaming_matches_in_memory(h, v2_subset_beds):
    """build_feature_matrix_streaming == build_feature_matrix on the same fixtures."""
    fs6 = list(v2_subset_beds)
    streamed = core.build_feature_matrix_streaming(
        {fs: iter_annotation_rows(v2_subset_beds[fs]) for fs in fs6}, h
    )
    in_memory = core.build_feature_matrix(
        {fs: read_annotation_bed(v2_subset_beds[fs]) for fs in fs6}, h
    )
    assert streamed.columns == in_memory.columns
    assert streamed.rows == in_memory.rows
    assert streamed.thresholds == in_memory.thresholds


def test_streaming_matches_in_memory_with_interspersion(h, v2_subset_beds):
    """Equivalence also holds when an interspersion featureset is requested."""
    fs6 = list(v2_subset_beds)
    streamed = core.build_feature_matrix_streaming(
        {fs: iter_annotation_rows(v2_subset_beds[fs]) for fs in fs6},
        h,
        interspersion_featureset="region",
    )
    in_memory = core.build_feature_matrix(
        {fs: read_annotation_bed(v2_subset_beds[fs]) for fs in fs6},
        h,
        interspersion_featureset="region",
    )
    assert streamed.rows == in_memory.rows


def test_iter_aligned_groups_yields_whole_sequences():
    """Each yielded group holds one sequence's intervals from every stream."""
    streams = {
        "region": iter([("r1", 0, 5, "arm"), ("r1", 5, 10, "bSat"), ("r2", 0, 10, "arm")]),
        "repeat": iter([("r1", 0, 10, "LINE"), ("r2", 0, 4, "LINE"), ("r2", 4, 10, "nonrepeat")]),
    }
    groups = list(iter_aligned_groups(streams))
    assert [seq_id for seq_id, _ in groups] == ["r1", "r2"]
    assert groups[0][1] == {
        "region": [(0, 5, "arm"), (5, 10, "bSat")],
        "repeat": [(0, 10, "LINE")],
    }
    assert groups[1][1]["repeat"] == [(0, 4, "LINE"), (4, 10, "nonrepeat")]


def test_iter_aligned_groups_order_mismatch():
    streams = {
        "region": iter([("A", 0, 10, "arm"), ("B", 0, 10, "arm")]),
        "repeat": iter([("B", 0, 10, "LINE"), ("A", 0, 10, "LINE")]),
    }
    with pytest.raises(ValueError, match="same order"):
        list(iter_aligned_groups(streams))


def test_iter_aligned_groups_early_exhaustion():
    streams = {
        "region": iter([("r1", 0, 10, "arm"), ("r2", 0, 10, "arm")]),
        "repeat": iter([("r1", 0, 10, "LINE")]),  # missing r2
    }
    with pytest.raises(ValueError, match="ran out of sequences"):
        list(iter_aligned_groups(streams))


def test_streaming_unknown_interspersion_featureset(h, v2_subset_beds):
    with pytest.raises(ValueError, match="interspersion featureset"):
        core.build_feature_matrix_streaming(
            {fs: iter_annotation_rows(v2_subset_beds[fs]) for fs in v2_subset_beds},
            h,
            interspersion_featureset="nope",
        )
