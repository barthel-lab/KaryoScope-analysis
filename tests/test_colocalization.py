"""Tests for the alignment-free colocalization measurement layer (Engine A)."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.core import colocalization as co
from karyoscope_analysis.core.io.bed import Interval

# read1 over [0, 35): arm, telomere, arm, rDNA
READ = [
    (0, 10, "arm"),
    (10, 12, "telomere"),
    (12, 30, "arm"),
    (30, 35, "rDNA"),
]


def _reverse(intervals: list[Interval]) -> list[Interval]:
    """Reflect a read's coordinates (its reverse orientation), keeping feature labels."""
    span = co.read_span(intervals)
    return [(span - end, span - start, feat) for start, end, feat in reversed(intervals)]


def test_feature_pair_is_canonical():
    assert co.feature_pair("b", "a") == ("a", "b")
    assert co.feature_pair("a", "b") == ("a", "b")


def test_read_span():
    assert co.read_span(READ) == 35
    assert co.read_span([(100, 110, "x"), (110, 130, "y")]) == 30
    assert co.read_span([]) == 0


def test_min_gaps_basic():
    gaps = co.min_gaps(READ)
    # adjacent pairs -> 0; rDNA and telomere are 18 bp apart (12->30)
    assert gaps == {
        ("arm", "telomere"): 0,
        ("arm", "rDNA"): 0,
        ("rDNA", "telomere"): 18,
    }


def test_min_gaps_is_intervening_distance():
    # A and B separated by a C block -> gap is the bp between them.
    ivals = [(0, 100, "A"), (100, 160, "C"), (160, 200, "B")]
    assert co.min_gaps(ivals)[("A", "B")] == 60  # 160 - 100
    assert co.min_gaps(ivals)[("A", "C")] == 0
    assert co.min_gaps(ivals)[("B", "C")] == 0


def test_min_gaps_takes_closest_approach():
    # A occurs twice; the gap to B is the nearest approach (right side here).
    ivals = [(0, 5, "A"), (5, 60, "filler"), (60, 65, "B"), (65, 100, "A")]
    # left A[0,5]->B[60,65] = 55 ; right A[65,100] is adjacent to B -> 0
    assert co.min_gaps(ivals)[("A", "B")] == 0


def test_min_occurrence_bp_filters_short_features():
    # telomere is only 2 bp -> dropped at threshold 3, so its pairs vanish.
    gaps = co.min_gaps(READ, min_occurrence_bp=3)
    assert gaps == {("arm", "rDNA"): 0}


def test_orientation_invariance():
    assert co.min_gaps(READ) == co.min_gaps(_reverse(READ))
    ivals = [(0, 100, "A"), (100, 160, "C"), (160, 200, "B")]
    assert co.min_gaps(ivals) == co.min_gaps(_reverse(ivals))


def test_iter_read_gaps_streams_a_bed(tmp_path: Path):
    bed = tmp_path / "overlay.bed"
    bed.write_text(
        "read1\t0\t10\tarm\n"
        "read1\t10\t12\ttelomere\n"
        "read1\t12\t30\tarm\n"
        "read2\t0\t5\trDNA\n"
        "read2\t5\t40\tarm\n"
        "read2\t40\t45\taSat\n"
    )
    out = list(co.iter_read_gaps(str(bed)))
    assert [seq_id for seq_id, _, _ in out] == ["read1", "read2"]

    (_, span1, gaps1), (_, span2, gaps2) = out
    assert span1 == 30
    assert gaps1 == {("arm", "telomere"): 0}
    assert span2 == 45
    assert gaps2 == {
        ("arm", "rDNA"): 0,
        ("aSat", "arm"): 0,
        ("aSat", "rDNA"): 35,  # rDNA[0,5] -> aSat[40,45]
    }
