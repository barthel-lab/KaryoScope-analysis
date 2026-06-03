"""Tests for the built-in overlay-annotations presets.

Loading each preset against the v2 hierarchy fixture validates that it is
well-formed AND v2-clean: any leftover v1 name (e.g. ``telomere_like_multigroup1``,
``arm_multigroup1``) would fail hierarchy validation and break these tests.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.core.annotation_resolution import (
    builtin_preset_names,
    load_builtin_preset,
)
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

EXPECTED_PRESETS = {
    "telomere-satellite",
    "priority",
    "chromosome-acrocentric",
    "telomere-acrocentric",
}


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def test_preset_names():
    assert set(builtin_preset_names()) == EXPECTED_PRESETS


@pytest.mark.parametrize("name", sorted(EXPECTED_PRESETS))
def test_preset_loads_and_is_v2_valid(name, h):
    # Must not raise: validates structure + every featureset/feature/@class against the hierarchy.
    spec = load_builtin_preset(name, h)
    assert spec.name == name
    assert spec.precedence  # non-empty


def test_telomere_acrocentric_composite(h):
    spec = load_builtin_preset("telomere-acrocentric", h)
    assert spec.resolve({"acrocentric": "DJ", "subtelomeric": "TAR1"}) == "DJ_TAR1"
    assert (
        spec.resolve({"acrocentric": "DJ", "subtelomeric": "canonical_telomere"})
        == "canonical_telomere"
    )
    assert spec.resolve({"acrocentric": "DJ", "subtelomeric": "nonsubtelomeric"}) == "DJ"


def test_priority_rules(h):
    spec = load_builtin_preset("priority", h)
    # arm is background -> repeat
    assert (
        spec.resolve({"region": "p_arm", "repeat": "LINE", "subtelomeric": "nonsubtelomeric"})
        == "LINE"
    )
    # ct over nonrepeat stays ct
    assert (
        spec.resolve({"region": "ct", "repeat": "nonrepeat", "subtelomeric": "nonsubtelomeric"})
        == "ct"
    )
    # subtel telomere priority overrides region/repeat
    assert (
        spec.resolve({"region": "bSat", "repeat": "LINE", "subtelomeric": "canonical_telomere"})
        == "canonical_telomere"
    )
    # satellite default -> region
    assert (
        spec.resolve({"region": "bSat", "repeat": "LINE", "subtelomeric": "nonsubtelomeric"})
        == "bSat"
    )


def test_chromosome_acrocentric(h):
    spec = load_builtin_preset("chromosome-acrocentric", h)
    assert spec.resolve({"chromosome": "chr13", "acrocentric": "DJ"}) == "DJ"
    assert spec.resolve({"chromosome": "chr13", "acrocentric": "nonacrocentric"}) == "chr13"


def test_telomere_satellite(h):
    spec = load_builtin_preset("telomere-satellite", h)
    assert (
        spec.resolve({"region": "bSat", "subtelomeric": "canonical_telomere"})
        == "canonical_telomere"
    )
    assert spec.resolve({"region": "bSat", "subtelomeric": "nonsubtelomeric"}) == "bSat"
