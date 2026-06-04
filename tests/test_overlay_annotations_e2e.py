"""End-to-end tests for ``overlay-annotations`` against real KS_human_CHM13_v2 data.

Fast tests run by default on the tiny committed ``tests/data/v2_subset`` HeLa fixtures
(see that directory's README for provenance). The ``@pytest.mark.integration`` tests
run the same pipeline on the full per-sample BEDs in ``data/raw_bed/`` and are skipped
when that (large, uncommitted) data isn't present; enable with ``pytest -m integration``.
"""

from __future__ import annotations

from itertools import pairwise
from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core.io.bed import read_annotation_bed

# Allowed feature literals an overlay may introduce that aren't raw input features.
# (The priority preset can emit a bare ``ct`` literal; ``ct`` is also a region feature.)
_EXTRA_LITERALS = {"ct"}


def _bed_args(beds: dict[str, Path]) -> list[str]:
    args: list[str] = []
    for fs, path in beds.items():
        args += ["--bed", f"{fs}={path}"]
    return args


def _seq_ids(path: Path) -> set[str]:
    return set(read_annotation_bed(path))


def _input_features(beds: dict[str, Path]) -> set[str]:
    feats: set[str] = set()
    for path in beds.values():
        for ivals in read_annotation_bed(path).values():
            feats.update(feat for _, _, feat in ivals)
    return feats


def _assert_coalesced(groups: dict[str, list[tuple[int, int, str]]]) -> None:
    """No two adjacent intervals within a sequence share the same feature."""
    for seq_id, ivals in groups.items():
        for (_, _, prev), (_, _, cur) in pairwise(ivals):
            assert prev != cur, (
                f"{seq_id}: adjacent intervals share feature {cur!r} (not coalesced)"
            )


# --------------------------------------------------------------- fixture (fast) tests
def test_priority_preset_on_subset(cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path):
    """The `priority` preset (region+repeat+subtelomeric) resolves to a valid BED."""
    inputs = {fs: v2_subset_beds[fs] for fs in ("region", "repeat", "subtelomeric")}
    out = tmp_path / "overlay.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--preset",
            "priority",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output

    # Output re-reads cleanly => it satisfies the C4 partition invariant.
    resolved = read_annotation_bed(out)
    # seq_ids are conserved exactly (overlay neither drops nor invents sequences).
    assert set(resolved) == _seq_ids(inputs["region"])
    # Every resolved feature comes from the inputs (+ the allowed `ct` literal).
    allowed = _input_features(inputs) | _EXTRA_LITERALS
    assert {feat for ivals in resolved.values() for _, _, feat in ivals} <= allowed
    _assert_coalesced(resolved)


def test_telomere_satellite_preset_on_subset(
    cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path
):
    """A second preset (`telomere-satellite`, region+subtelomeric) also resolves cleanly."""
    inputs = {fs: v2_subset_beds[fs] for fs in ("region", "subtelomeric")}
    out = tmp_path / "overlay.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--preset",
            "telomere-satellite",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    resolved = read_annotation_bed(out)
    assert set(resolved) == _seq_ids(inputs["region"])
    _assert_coalesced(resolved)


def test_default_basic_overlay_composites(
    cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path
):
    """With no preset/spec, the default overlay joins featuresets into `a:b` composites."""
    inputs = {fs: v2_subset_beds[fs] for fs in ("region", "subtelomeric")}
    out = tmp_path / "overlay.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--separator",
            ":",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    resolved = read_annotation_bed(out)
    assert set(resolved) == _seq_ids(inputs["region"])
    feats = {feat for ivals in resolved.values() for _, _, feat in ivals}
    # Every composite label is `region_feature:subtelomeric_feature`.
    region_feats = _input_features({"region": inputs["region"]})
    subtel_feats = _input_features({"subtelomeric": inputs["subtelomeric"]})
    for feat in feats:
        left, sep, right = feat.partition(":")
        assert sep == ":", f"expected a composite label, got {feat!r}"
        assert left in region_feats and right in subtel_feats, feat
    _assert_coalesced(resolved)


def test_seq_id_with_slash_roundtrips(cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path):
    """PacBio `.../ccs` read ids (containing `/`) survive the overlay unchanged."""
    inputs = {fs: v2_subset_beds[fs] for fs in ("region", "repeat", "subtelomeric")}
    out = tmp_path / "overlay.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--preset",
            "priority",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    resolved = read_annotation_bed(out)
    assert any("/ccs" in seq_id for seq_id in resolved), "slash-style read id was lost"


# ----------------------------------------------------------------- integration (full)
@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("sample", ["HeLa", "IMR90", "U2OS"])
def test_priority_preset_full_data(cli_runner, hierarchy_tsv, raw_bed_lookup, sample, tmp_path):
    """The `priority` preset runs on each full sample and conserves the read set."""
    beds = raw_bed_lookup(sample)
    if beds is None:
        pytest.skip(f"full v2 raw BEDs for {sample} not present in data/raw_bed/")
    inputs = {fs: beds[fs] for fs in ("region", "repeat", "subtelomeric")}
    out = tmp_path / "overlay.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--preset",
            "priority",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    resolved = read_annotation_bed(out)  # C4-valid by construction of the read
    assert set(resolved) == _seq_ids(inputs["region"])
    _assert_coalesced(resolved)
