"""End-to-end tests for ``build-feature-matrix`` against real KS_human_CHM13_v2 data.

Fast tests run by default on the tiny committed ``tests/data/v2_subset`` HeLa fixtures.
The ``@pytest.mark.integration`` test runs on the full ``data/raw_bed/`` HeLa BEDs and
is skipped when that (large, uncommitted) data isn't present (``pytest -m integration``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import seq_features as sf
from karyoscope_analysis.core.io.bed import read_annotation_bed

_METRICS = {
    "frac",
    "bp",
    "dmax",
    "dmin",
    "dmedian",
    "dfirst",
    "dlast",
    "dterminal",
    "dterminal_min",
    "max_block_bp",
}
_INTERSPERSION_TYPES = {"total", "can_ncan", "tel_sat", "arm_tel"}


def _read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        cells = line.split("\t")
        rows[cells[0]] = dict(zip(header, cells, strict=True))
    return header, rows


def _assert_valid_schema(header: list[str], featuresets: set[str]) -> None:
    """Every column matches the F2 schema (``__`` is the sole delimiter)."""
    assert header[0] == "seq_id"
    for col in header[1:]:
        if col.startswith("interspersion__"):
            assert col.split("__", 1)[1] in _INTERSPERSION_TYPES, col
        elif col.endswith("__total_bp"):
            assert col[: -len("__total_bp")] in featuresets, col
        else:
            fs, metric, _feat = col.split("__", 2)
            assert fs in featuresets, col
            assert metric in _METRICS, col


def _bed_args(beds: dict[str, Path]) -> list[str]:
    args: list[str] = []
    for fs, path in beds.items():
        args += ["--bed", f"{fs}={path}"]
    return args


# --------------------------------------------------------------- fixture (fast) tests
def test_matrix_schema_and_coverage(cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path):
    """Six featuresets -> a per-read matrix whose coverage columns are self-consistent."""
    out = tmp_path / "matrix.tsv"
    result = cli_runner.invoke(
        main,
        [
            "build-feature-matrix",
            *_bed_args(v2_subset_beds),
            "--hierarchy",
            str(hierarchy_tsv),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output

    header, rows = _read_matrix(out)
    featuresets = set(v2_subset_beds)
    _assert_valid_schema(header, featuresets)

    # One row per input sequence, no extras/drops.
    region_groups = read_annotation_bed(v2_subset_beds["region"])
    assert set(rows) == set(region_groups)

    for seq_id, ivals in region_groups.items():
        row = rows[seq_id]
        expected_total = sf.total_bp(ivals)
        # total_bp matches the sequence span (C4 partition).
        assert int(row["region__total_bp"]) == expected_total
        # Per-feature bp sums back to total_bp; fracs sum to 1.
        bp_cols = {k: v for k, v in row.items() if k.startswith("region__bp__")}
        frac_cols = {k: v for k, v in row.items() if k.startswith("region__frac__")}
        assert sum(int(v) for v in bp_cols.values()) == expected_total
        assert sum(float(v) for v in frac_cols.values()) == pytest.approx(1.0)
        # Coverage fractions are in [0, 1].
        assert all(0.0 <= float(v) <= 1.0 for v in frac_cols.values())


def test_matrix_threshold_sidecar(cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path):
    """The adaptive-threshold sidecar is written with a row per featureset/feature."""
    out = tmp_path / "matrix.tsv"
    result = cli_runner.invoke(
        main,
        [
            "build-feature-matrix",
            *_bed_args(v2_subset_beds),
            "--hierarchy",
            str(hierarchy_tsv),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    sidecar = tmp_path / "matrix.adaptive_thresholds.tsv"
    assert sidecar.is_file()
    lines = sidecar.read_text().splitlines()
    assert lines[0] == "featureset\tfeature\tthreshold"
    featuresets_seen = {line.split("\t")[0] for line in lines[1:]}
    # Every input featureset contributes at least one threshold row.
    assert set(v2_subset_beds) <= featuresets_seen
    # Thresholds parse as floats within the configured clamp range.
    for line in lines[1:]:
        _fs, _feat, thr = line.split("\t")
        assert sf.DEFAULT_THRESHOLD_MIN <= float(thr) <= sf.DEFAULT_THRESHOLD_MAX


def test_interspersion_over_overlay(cli_runner, hierarchy_tsv, v2_subset_beds, tmp_path: Path):
    """overlay-annotations -> build-feature-matrix with interspersion over the composite."""
    # Step 1: build a priority overlay (telomere + satellite + background on one track).
    overlay = tmp_path / "overlay.bed.gz"
    ov_inputs = {fs: v2_subset_beds[fs] for fs in ("region", "repeat", "subtelomeric")}
    r1 = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            *_bed_args(ov_inputs),
            "--hierarchy",
            str(hierarchy_tsv),
            "--preset",
            "priority",
            "-o",
            str(overlay),
        ],
    )
    assert r1.exit_code == 0, r1.output

    # Step 2: matrix over the six featuresets + the overlay, interspersion over overlay.
    out = tmp_path / "matrix.tsv"
    beds = {**v2_subset_beds, "overlay": overlay}
    r2 = cli_runner.invoke(
        main,
        [
            "build-feature-matrix",
            *_bed_args(beds),
            "--hierarchy",
            str(hierarchy_tsv),
            "--interspersion-featureset",
            "overlay",
            "-o",
            str(out),
        ],
    )
    assert r2.exit_code == 0, r2.output

    header, rows = _read_matrix(out)
    _assert_valid_schema(header, set(beds))
    # All four interspersion columns are present and non-negative.
    inter_cols = [f"interspersion__{t}" for t in _INTERSPERSION_TYPES]
    assert all(col in header for col in inter_cols)
    for row in rows.values():
        assert all(float(row[col]) >= 0.0 for col in inter_cols)
    # The telomere-rich read shows canonical<->noncanonical telomere transitions.
    tel_read = "ee18d619-9d89-4858-b420-b6fd840107af"
    assert float(rows[tel_read]["interspersion__can_ncan"]) > 0.0


# ----------------------------------------------------------------- integration (full)
@pytest.mark.integration
@pytest.mark.slow
def test_matrix_full_hela(cli_runner, hierarchy_tsv, raw_bed_lookup, tmp_path: Path):
    """The full HeLa matrix builds, with one row per read and a valid schema."""
    beds = raw_bed_lookup("HeLa")
    if beds is None:
        pytest.skip("full v2 raw BEDs for HeLa not present in data/raw_bed/")
    out = tmp_path / "matrix.tsv"
    result = cli_runner.invoke(
        main,
        [
            "build-feature-matrix",
            *_bed_args(beds),
            "--hierarchy",
            str(hierarchy_tsv),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    header, rows = _read_matrix(out)
    _assert_valid_schema(header, set(beds))
    assert set(rows) == set(read_annotation_bed(beds["region"]))
