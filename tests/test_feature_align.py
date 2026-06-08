"""Tests for the Engine B feature-sequence aligner."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.core import feature_align as fa
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

# A simple exact-or-mismatch scorer for the alignment-mechanics tests.
EXACT = lambda x, y: 1.0 if x == y else -1.0  # noqa: E731


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def test_to_and_reverse_segments():
    intervals = [(0, 10, "arm"), (10, 30, "bSat")]
    assert fa.to_segments(intervals) == [("arm", 10), ("bSat", 20)]
    assert fa.reverse_segments([("arm", 10), ("bSat", 20)]) == [("bSat", 20), ("arm", 10)]


def test_identical_sequences_are_contained_and_full_score():
    seq = [("X", 100), ("Y", 50)]
    aln = fa.align_local(seq, seq, sub_score=EXACT, gap_factor=0.01)
    assert aln.score == pytest.approx(150.0)  # 100 + 50 matched bp
    assert aln.columns == ((0, 0), (1, 1))
    assert fa.classify_overlap(aln, 2, 2) == "containment"
    assert fa.is_proper_overlap(aln, 2, 2)


def test_containment_of_a_block():
    a = [("X", 100)]
    b = [("W", 50), ("X", 100), ("Z", 50)]
    aln = fa.align_local(a, b, sub_score=EXACT, gap_factor=0.01)
    assert aln.score == pytest.approx(100.0)
    assert aln.columns == ((0, 1),)  # only the X block matches
    assert fa.classify_overlap(aln, len(a), len(b)) == "containment"


def test_dovetail_suffix_prefix():
    a = [("P", 30), ("Q", 100), ("R", 80)]
    b = [("Q", 100), ("R", 80), ("S", 40)]
    aln = fa.align_local(a, b, sub_score=EXACT, gap_factor=0.01)
    assert aln.score == pytest.approx(180.0)  # Q + R
    assert aln.a_end == 2 and aln.b_start == 0  # suffix of A meets prefix of B
    assert fa.classify_overlap(aln, len(a), len(b)) == "dovetail"


def test_unrelated_reads_do_not_align():
    aln = fa.align_local([("X", 100)], [("Y", 100)], sub_score=EXACT, gap_factor=0.01)
    assert aln.is_empty
    assert aln.score == 0.0
    assert fa.classify_overlap(aln, 1, 1) == "none"
    assert not fa.is_proper_overlap(aln, 1, 1)


def test_length_weighting_favours_large_blocks():
    big = fa.align_local([("X", 100)], [("X", 100)], sub_score=EXACT, gap_factor=0.01)
    small = fa.align_local([("X", 2)], [("X", 2)], sub_score=EXACT, gap_factor=0.01)
    assert big.score > small.score  # 100 vs 2


def test_length_mismatch_is_charged_as_a_gap():
    # matching X(12k) to X(15k): reward min(12k,15k) minus gap_factor*|12k-15k| (the length
    # difference is a gap), so gap_factor controls how strictly lengths must agree.
    strict = fa.align_local([("X", 12000)], [("X", 15000)], sub_score=EXACT, gap_factor=1.0)
    assert strict.score == pytest.approx(12000.0 - 3000.0)  # 12000 - 1.0*3000 = 9000
    lenient = fa.align_local([("X", 12000)], [("X", 15000)], sub_score=EXACT, gap_factor=0.01)
    assert lenient.score == pytest.approx(12000.0 - 30.0)  # small penalty at low gap_factor


def test_gap_bridges_an_inserted_block_when_cheap():
    a = [("X", 100), ("G", 50), ("Y", 100)]
    b = [("X", 100), ("Y", 100)]  # B lacks the G block
    aln = fa.align_local(a, b, sub_score=EXACT, gap_factor=0.01)
    # X (100) - skip G (0.01*50) + Y (100)
    assert aln.score == pytest.approx(199.5)
    assert aln.columns == ((0, 0), (2, 1))


def test_large_gap_factor_blocks_bridging():
    a = [("X", 100), ("G", 50), ("Y", 100)]
    b = [("X", 100), ("Y", 100)]
    aln = fa.align_local(a, b, sub_score=EXACT, gap_factor=10.0)
    # skipping G would cost 500 >> 100, so the best local match is a single block.
    assert aln.score == pytest.approx(100.0)
    assert len(aln.columns) == 1


def test_best_orientation_picks_reverse():
    a = [("X", 100), ("Y", 50)]
    b = [("Y", 50), ("X", 100)]  # reverse of a
    aln = fa.align_best_orientation(a, b, sub_score=EXACT, gap_factor=0.01)
    assert aln.reversed_b
    assert aln.score == pytest.approx(150.0)


# ----------------------------------------------------------------- hierarchy substitution
def test_hierarchy_substitution_tiers(h):
    score = fa.hierarchy_substitution(h, match=1.0, partial=0.5, mismatch=-1.0)
    assert score("bSat", "bSat") == 1.0  # exact
    assert score("bSat", "gSat") == 0.5  # both satellites -> partial
    assert score("canonical_telomere", "noncanonical_telomere") == 0.5  # telomere group
    assert score("arm", "canonical_telomere") == -1.0  # different groups -> mismatch
    assert score("novel", "arm") == 0.0  # novel is neutral
    assert score("arm", "novel") == 0.0


def test_chromosome_aware_substitution(h):
    struct = fa.hierarchy_substitution(h, match=1.0, partial=0.5, mismatch=-1.0)
    score = fa.chromosome_aware_substitution(struct, cross_chromosome_penalty=-2.0)
    # same chromosome: structural layer decides
    assert score("chr6:aSat", "chr6:aSat") == 1.0  # exact
    assert score("chr6:aSat", "chr6:bSat") == 0.5  # sibling satellites -> partial
    # different specific chromosomes: structural match + cross-chromosome penalty
    assert score("chr6:aSat", "chr21:aSat") == pytest.approx(1.0 - 2.0)
    # ambiguous chromosome label (not chr*): penalized (can't confirm same chromosome)
    assert score("chr6:aSat", "autosome:aSat") == pytest.approx(1.0 - 2.0)
    assert score("autosome:aSat", "autosome:aSat") == pytest.approx(1.0 - 2.0)
    # no chromosome layer -> degrades to the structural scorer (no penalty)
    assert score("aSat", "bSat") == 0.5
    assert score("aSat", "arm") == -1.0


def test_feature_jaccard_prefilter():
    a = [("X", 10), ("Y", 10)]
    b = [("Y", 10), ("Z", 10)]
    assert fa.feature_jaccard(a, b) == pytest.approx(1 / 3)  # {X,Y} vs {Y,Z}
    assert fa.feature_jaccard(a, []) == 0.0
