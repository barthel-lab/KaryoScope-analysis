"""Tests for the Engine A differential rearrangement test."""

from __future__ import annotations

import math

import pytest

from karyoscope_analysis.core import rearrangement as rr
from karyoscope_analysis.core.colocalization import FeaturePair

AB: FeaturePair = ("A", "B")  # canonical sorted pairs
CD: FeaturePair = ("C", "D")
EF: FeaturePair = ("E", "F")


def _sample(supported: dict[FeaturePair, int], n_reads: int, span: int = 5000, boundaries=()):
    """Build a SampleColocalization: `supported[pair]` reads carry that pair at gap 0."""
    reads: list[dict[FeaturePair, int]] = [{} for _ in range(n_reads)]
    for pair, count in supported.items():
        for i in range(count):
            reads[i][pair] = 0
    records = [(f"r{i}", span, d) for i, d in enumerate(reads)]
    return rr.aggregate(records, boundaries)


# ----------------------------------------------------------------- bucketing / aggregate
def test_length_bucket():
    assert rr.length_bucket(100, ()) == 0
    assert rr.length_bucket(100, (25_000,)) == 0
    assert rr.length_bucket(25_000, (25_000,)) == 1  # boundary goes to the upper bucket
    assert rr.length_bucket(60_000, (25_000, 50_000)) == 2


def test_aggregate_support_and_rate():
    s = _sample({AB: 20, CD: 10}, 100)
    assert s.total(0) == 100
    assert s.support(0, AB, 0) == 20
    assert s.support(0, AB, 1000) == 20  # all supported reads are at gap 0
    assert s.support(0, ("X", "Y"), 0) == 0
    assert s.pooled_rate(AB, 0) == pytest.approx(0.20)
    assert s.pairs() == {AB, CD}


def test_support_counts_only_within_window():
    # one read at gap 0, one at gap 5000
    records = [("r0", 9999, {AB: 0}), ("r1", 9999, {AB: 5000})]
    s = rr.aggregate(records)
    assert s.support(0, AB, 0) == 1
    assert s.support(0, AB, 4999) == 1
    assert s.support(0, AB, 5000) == 2


# ----------------------------------------------------------------- statistics
def test_cmh_single_stratum_odds_ratio():
    # a=20 b=80 c=2 d=98 -> OR = ad/bc
    p, odds = rr._cmh([(20, 80, 2, 98)])
    assert odds == pytest.approx(20 * 98 / (80 * 2))
    assert p < 0.001  # strong enrichment


def test_cmh_null_is_nonsignificant():
    p, odds = rr._cmh([(10, 90, 10, 90), (10, 90, 10, 90)])
    assert odds == pytest.approx(1.0, abs=1e-9)
    assert p > 0.5


def test_cmh_stratified_combines_evidence():
    # modest enrichment in each of two strata -> combined is significant
    p, odds = rr._cmh([(20, 80, 8, 92), (18, 82, 7, 93)])
    assert odds > 1.5
    assert p < 0.01


def test_log2_ratio_handles_zero_control():
    # control support 0 must not blow up (Haldane correction)
    val = rr._log2_ratio(20, 100, 0, 100)
    assert math.isfinite(val)
    assert val > 0


# ----------------------------------------------------------------- end-to-end detection
def test_detect_flags_enriched_novel_pair_and_respects_gates():
    experiment = _sample({AB: 20, CD: 10, EF: 30}, 100)
    control = _sample({AB: 2, CD: 11, EF: 5}, 100)
    reference = _sample({CD: 10, EF: 35}, 100)  # AB absent in reference; EF high (artifact floor)

    calls = rr.detect_rearrangements(
        experiment, control, reference, windows=(0,), min_support=3, min_log2_ratio=1.0
    )
    by_pair = {c.pair: c for c in calls}

    # AB: novel (absent in reference), strongly enriched -> a call.
    ab = by_pair[AB]
    assert ab.passes
    assert ab.direction == "enriched"
    assert ab.reference_abnormal  # ref_rate == 0
    assert ab.exp_rate == pytest.approx(0.20)
    assert ab.ctrl_rate == pytest.approx(0.02)
    assert ab.q_value <= 0.05

    # CD: no real difference -> fails (effect too small / not significant).
    assert not by_pair[CD].passes

    # EF: statistically enriched vs control AND large effect, but does NOT clear the
    # reference (artifact) floor (ref_rate 0.35 > exp_rate 0.30) -> suppressed.
    ef = by_pair[EF]
    assert ef.q_value <= 0.05  # would be "significant"...
    assert abs(ef.log2_ratio) >= 1.0  # ...with a large effect...
    assert ef.ref_rate == pytest.approx(0.35)
    assert not ef.passes  # ...but the floor gate rejects it
    assert not ef.reference_abnormal


def test_detect_without_reference_uses_zero_floor():
    experiment = _sample({AB: 20}, 100)
    control = _sample({AB: 2}, 100)
    calls = rr.detect_rearrangements(experiment, control, windows=(0,))
    ab = next(c for c in calls if c.pair == AB)
    assert ab.ref_rate == 0.0
    assert ab.reference_abnormal
    assert ab.passes


def test_detect_passing_calls_sorted_first():
    experiment = _sample({AB: 20, CD: 10}, 100)
    control = _sample({AB: 2, CD: 11}, 100)
    calls = rr.detect_rearrangements(experiment, control, windows=(0,))
    passing = [c for c in calls if c.passes]
    # passing calls come first
    assert calls[: len(passing)] == passing
