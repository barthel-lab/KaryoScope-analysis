"""Tests for the per-sequence feature metrics (build-feature-matrix core)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from karyoscope_analysis.core import seq_features as sf
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


# A small C4 partition of read1 over [0, 6).
PARTITION = [(0, 3, "bSat"), (3, 4, "ct"), (4, 6, "bSat")]


def test_coverage():
    assert sf.total_bp(PARTITION) == 6
    assert sf.feature_bp(PARTITION) == {"bSat": 5, "ct": 1}
    assert sf.feature_fraction(PARTITION) == {"bSat": 5 / 6, "ct": 1 / 6}


def test_max_block_bp_gap_bridging():
    cov = np.array([1, 1, 1, 0, 1, 1])  # two runs (len 3, len 2) split by a 1-bp gap
    assert sf.max_block_bp(cov, gap_tol=0) == 3  # not merged
    assert sf.max_block_bp(cov, gap_tol=1) == 6  # gap <= tol -> merged into one run
    assert sf.max_block_bp(np.zeros(5, dtype=int)) == 0


def test_window_densities_short_sequence_broadcasts():
    # span 6 < window 1000 -> whole-sequence fraction broadcast; default gap_tol bridges the gap.
    stats = sf.window_densities(PARTITION, window_size=1000)
    bsat = stats["bSat"]
    assert bsat.dmax == bsat.dmin == bsat.dmedian == pytest.approx(5 / 6)
    assert bsat.dterminal == pytest.approx(5 / 6)
    assert bsat.max_block_bp == 6  # default gap_tol=100 bridges the 1-bp ct gap


def test_window_densities_sliding():
    # bSat coverage = [1,1,1,0,1,1]; window=2 -> window sums [2,2,1,1,2] -> /2.
    stats = sf.window_densities(PARTITION, window_size=2, gap_tol=0)
    bsat = stats["bSat"]
    assert bsat.dmax == pytest.approx(1.0)
    assert bsat.dmin == pytest.approx(0.5)
    assert bsat.dmedian == pytest.approx(1.0)  # median([2,2,1,1,2]) = 2 -> /2
    assert bsat.dfirst == pytest.approx(1.0)  # first 2 bp both covered
    assert bsat.dlast == pytest.approx(1.0)
    assert bsat.max_block_bp == 3  # gap_tol=0 -> longest run is [0,3)


def test_classify_feature(h):
    assert sf.classify_feature("bSat", h) == "satellite"
    assert sf.classify_feature("ct", h) == "ct"
    assert sf.classify_feature("arm", h) == "arm"
    assert sf.classify_feature("p_arm", h) == "arm"
    assert sf.classify_feature("canonical_telomere", h) == "canonical"
    assert sf.classify_feature("noncanonical_telomere", h) == "noncanonical"
    assert sf.classify_feature("ITS", h) == "ITS_TAR1"
    assert sf.classify_feature("TAR1", h) == "ITS_TAR1"
    assert sf.classify_feature("nonrepeat", h) == "other"
    # two-layer: layer-2 telomere overrides a layer-1 arm; layer-1 satellite wins outright
    assert sf.classify_feature("arm:canonical_telomere", h) == "canonical"
    assert sf.classify_feature("bSat:LINE", h) == "satellite"


def test_interspersion(h):
    # span exactly 1000 bp -> per-kb values equal raw transition counts
    intervals = [
        (0, 400, "canonical_telomere"),
        (400, 600, "bSat"),
        (600, 1000, "canonical_telomere"),
    ]
    result = sf.interspersion(intervals, h)
    assert result == {"total": 2.0, "can_ncan": 0.0, "tel_sat": 2.0, "arm_tel": 0.0}


def test_interspersion_empty(h):
    assert sf.interspersion([], h) == {
        "total": 0.0,
        "can_ncan": 0.0,
        "tel_sat": 0.0,
        "arm_tel": 0.0,
    }


def test_adaptive_thresholds_clamps():
    rows = [{"bSat": 0.3}, {"bSat": 0.6}]
    # median(0.3, 0.6)=0.45 ; /3 = 0.15 ; clamp to max 0.05
    assert sf.adaptive_thresholds(rows) == {"bSat": pytest.approx(0.05)}


def test_adaptive_thresholds_floor_and_midrange():
    out = sf.adaptive_thresholds([{"x": 0.003}, {"z": 0.012}, {"z": 0.0}])
    assert out["x"] == pytest.approx(0.001)  # 0.003/3 = 0.001 (== floor)
    assert out["z"] == pytest.approx(0.004)  # median([0.012])/3, within bounds


def test_adaptive_thresholds_all_zero_gets_floor():
    assert sf.adaptive_thresholds([{"y": 0.0}, {"y": 0.0}]) == {"y": pytest.approx(0.001)}
