"""Tests for the hierarchy-derived v2 feature vocabulary.

Run against the committed v2 ``hierarchy.tsv`` fixture, asserting the decisions
locked in docs/audit/DECISIONS.md: derive vocab from the hierarchy (D4.2),
satellites = centromeric subtree minus ``ct`` (D4.4), v2-only with v1 names
rejected (D4.1), and ``novel`` is the only out-of-taxonomy feature allowed (C2).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.core.feature_vocab import NOVEL, FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

_CENTROMERIC_SUBTREE = {
    "aSat",
    "alpha_hor",
    "active_hor",
    "dhor",
    "hor",
    "mixedAlpha",
    "mon",
    "bSat",
    "cenSat",
    "ct",
    "gSat",
    "HSat",
    "HSat1",
    "HSat1A",
    "HSat1B",
    "HSat2",
    "HSat3",
}


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def test_feature_sets(h: FeatureHierarchy) -> None:
    assert h.feature_sets() == {
        "acrocentric",
        "chromosome",
        "gene",
        "region",
        "repeat",
        "subtelomeric",
    }


def test_descendants(h: FeatureHierarchy) -> None:
    d = h.descendants("region", "centromeric")
    assert d == _CENTROMERIC_SUBTREE
    assert "centromeric" not in d
    assert h.descendants("region", "centromeric", include_root=True) == d | {"centromeric"}
    assert "arm" not in d  # different subtree


def test_children(h: FeatureHierarchy) -> None:
    assert h.children("region", "HSat") == {"HSat1", "HSat2", "HSat3"}
    assert h.children("region", "HSat1") == {"HSat1A", "HSat1B"}
    assert h.children("region", "HSat1A") == frozenset()  # leaf


def test_satellite_features_exclude_ct(h: FeatureHierarchy) -> None:
    sats = h.satellite_features
    assert "centromeric" in sats
    assert {"aSat", "bSat", "gSat", "cenSat", "HSat3", "mon", "alpha_hor"} <= sats
    assert "ct" not in sats  # decision D4.4
    assert "arm" not in sats and "p_arm" not in sats


def test_arm_and_ct(h: FeatureHierarchy) -> None:
    assert h.arm_features == {"arm", "p_arm", "q_arm"}
    assert h.ct_features == {"ct"}


def test_telomere_groups(h: FeatureHierarchy) -> None:
    assert h.canonical_telomere == {"canonical_telomere"}
    assert h.noncanonical_telomere == {"noncanonical_telomere"}
    assert h.its_tar1 == {"ITS", "TAR1"}


def test_membership_is_v2_only(h: FeatureHierarchy) -> None:
    assert "bSat" in h
    assert "categorized" in h  # background label is a real feature value
    assert "canonical_telomere" in h
    # v1 names are NOT recognized (decision D4.1)
    assert "bsat" not in h
    assert "bsat_specific" not in h
    assert "arm_multigroup1" not in h


def test_validity_and_novel(h: FeatureHierarchy) -> None:
    assert h.is_valid_feature("bSat")
    assert h.is_valid_feature(NOVEL)  # C2: novel is always allowed
    assert not h.is_valid_feature("bsat")  # v1 rejected
    h.require_valid_feature("bSat")  # no raise
    h.require_valid_feature(NOVEL)  # no raise
    with pytest.raises(ValueError, match="unknown feature"):
        h.require_valid_feature("not_a_real_feature")


def test_per_featureset_features(h: FeatureHierarchy) -> None:
    region = h.features("region")
    assert "bSat" in region
    assert "canonical_telomere" not in region  # that's subtelomeric
    assert h.features("no_such_set") == frozenset()


def test_from_tsv_rejects_bad_columns(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    bad.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="expected a TSV"):
        FeatureHierarchy.from_tsv(bad)
