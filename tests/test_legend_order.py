"""Tests for the featureset-aware legend sort key (built on KaryoScope's legend_sort_key)."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.core.legend_order import feature_sort_key

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


def _sorted(names):
    key = feature_sort_key(HIERARCHY_TSV)
    return sorted(names, key=key)


def test_chromosomes_natural_within_chromosome_featureset():
    # chr2 before chr10 (natural, not lexical) and chrX/chrY after the numerics.
    assert _sorted(["chr10", "chr2", "chrX", "chr1"]) == ["chr1", "chr2", "chr10", "chrX"]


def test_features_group_by_featureset():
    # Telomere types (subtelomeric) group together; satellites (region) group together;
    # the two groups don't interleave regardless of input order.
    out = _sorted(["aSat", "canonical_telomere", "bSat", "noncanonical_telomere"])
    tel = {"canonical_telomere", "noncanonical_telomere"}
    idx = [i for i, n in enumerate(out) if n in tel]
    assert idx == [0, 1] or idx == [2, 3]  # contiguous block, not interleaved


def test_novel_sorts_last():
    out = _sorted(["novel", "aSat", "canonical_telomere"])
    assert out[-1] == "novel"


def test_composite_label_sorts_by_structural_layer():
    # "chr2:aSat" sorts as aSat (region/centromeric), grouped with bare aSat.
    key = feature_sort_key(HIERARCHY_TSV)
    assert key("chr2:aSat")[0] == key("aSat")[0]  # same featureset rank
