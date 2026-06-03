"""Tests for the overlay-annotations resolution engine (precedence + ordered rules)."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.core.annotation_resolution import (
    SpecError,
    load_spec,
    load_spec_file,
)
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


# A faithful telomere-acrocentric-style spec (precedence: acrocentric default).
TELO_ACRO = {
    "name": "telomere-acrocentric",
    "precedence": ["acrocentric", "subtelomeric"],
    "rules": [
        {
            "when": {"subtelomeric": ["canonical_telomere", "noncanonical_telomere"]},
            "emit": "subtelomeric",
        },
        {
            "when": {"acrocentric": ["DJ", "PHR", "rDNA"], "subtelomeric": ["ITS", "TAR1"]},
            "emit": "composite",
        },
        {"when": {"subtelomeric": ["ITS", "TAR1"]}, "emit": "subtelomeric"},
    ],
}


def test_precedence_default(h):
    spec = load_spec(TELO_ACRO, h)
    # nonsubtelomeric subtel + informative acro: no rule matches -> precedence default (acrocentric)
    assert spec.resolve({"acrocentric": "DJ", "subtelomeric": "nonsubtelomeric"}) == "DJ"
    # background everywhere -> acrocentric default
    assert (
        spec.resolve({"acrocentric": "nonacrocentric", "subtelomeric": "nonsubtelomeric"})
        == "nonacrocentric"
    )


def test_passthrough_override(h):
    spec = load_spec(TELO_ACRO, h)
    # canonical telomere overrides the acrocentric default
    assert (
        spec.resolve({"acrocentric": "DJ", "subtelomeric": "canonical_telomere"})
        == "canonical_telomere"
    )


def test_composite_uses_precedence_order(h):
    spec = load_spec(TELO_ACRO, h)
    # acrocentric precedes subtelomeric -> "DJ_TAR1" (matches the legacy label)
    assert spec.resolve({"acrocentric": "DJ", "subtelomeric": "TAR1"}) == "DJ_TAR1"
    # TAR1 over a non-composite acro -> keep subtelomeric (rule 3)
    assert spec.resolve({"acrocentric": "SST1", "subtelomeric": "ITS"}) == "ITS"


def test_first_match_wins(h):
    spec = load_spec(TELO_ACRO, h)
    # canonical_telomere (rule 1) beats the composite rule even with a composite-eligible acro
    assert (
        spec.resolve({"acrocentric": "rDNA", "subtelomeric": "canonical_telomere"})
        == "canonical_telomere"
    )


def test_at_class_matching(h):
    spec = load_spec(
        {
            "name": "arm-test",
            "precedence": ["region", "repeat"],
            "rules": [{"when": {"region": "@arm"}, "emit": "repeat"}],
        },
        h,
    )
    assert spec.resolve({"region": "p_arm", "repeat": "LINE"}) == "LINE"  # @arm covers p_arm
    assert spec.resolve({"region": "arm", "repeat": "SINE"}) == "SINE"
    assert (
        spec.resolve({"region": "bSat", "repeat": "LINE"}) == "bSat"
    )  # not an arm -> default region


def test_literal_emit(h):
    spec = load_spec(
        {
            "name": "ct-test",
            "precedence": ["region", "repeat"],
            "rules": [
                {"when": {"region": ["ct"], "repeat": ["nonrepeat"]}, "emit": {"literal": "ct"}}
            ],
        },
        h,
    )
    assert spec.resolve({"region": "ct", "repeat": "nonrepeat"}) == "ct"
    assert (
        spec.resolve({"region": "ct", "repeat": "LINE"}) == "ct"
    )  # default region == "ct" here too


def test_composite_all_is_basic_overlay(h):
    spec = load_spec(
        {
            "name": "overlay",
            "precedence": ["region", "repeat"],
            "rules": [{"emit": {"composite": "all"}}],  # no `when` -> matches everything
        },
        h,
    )
    assert spec.resolve({"region": "bSat", "repeat": "LINE"}) == "bSat:LINE"


# --- validation -----------------------------------------------------------------
def test_rejects_unknown_featureset(h):
    with pytest.raises(SpecError, match="unknown feature set"):
        load_spec({"name": "x", "precedence": ["not_a_set"]}, h)


def test_rejects_unknown_feature_in_when(h):
    with pytest.raises(SpecError, match="not in feature set"):
        load_spec(
            {
                "name": "x",
                "precedence": ["region"],
                "rules": [{"when": {"region": ["nope"]}, "emit": "region"}],
            },
            h,
        )


def test_rejects_unknown_class(h):
    with pytest.raises(SpecError, match="class @nope"):
        load_spec(
            {
                "name": "x",
                "precedence": ["region"],
                "rules": [{"when": {"region": "@nope"}, "emit": "region"}],
            },
            h,
        )


def test_rejects_when_keys_out_of_order(h):
    with pytest.raises(SpecError, match="precedence order"):
        load_spec(
            {
                "name": "x",
                "precedence": ["region", "repeat"],
                "rules": [{"when": {"repeat": ["LINE"], "region": ["bSat"]}, "emit": "region"}],
            },
            h,
        )


def test_rejects_emit_featureset_not_in_precedence(h):
    with pytest.raises(SpecError, match="not in precedence"):
        load_spec(
            {"name": "x", "precedence": ["region"], "rules": [{"emit": "repeat"}]},
            h,
        )


def test_rejects_malformed_schema(h):
    with pytest.raises(SpecError, match="not well-formed"):
        load_spec({"name": "x"}, h)  # missing required 'precedence'


def test_load_spec_file_yaml(h, tmp_path: Path):
    yaml_text = (
        "name: y\n"
        "precedence: [region, repeat]\n"
        "rules:\n"
        "  - when: {region: '@arm'}\n"
        "    emit: repeat\n"
    )
    path = tmp_path / "spec.yaml"
    path.write_text(yaml_text)
    spec = load_spec_file(path, h)
    assert spec.resolve({"region": "p_arm", "repeat": "LINE"}) == "LINE"
