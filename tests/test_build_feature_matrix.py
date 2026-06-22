"""End-to-end tests for build-feature-matrix (core orchestration + CLI)."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import build_feature_matrix as core
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import read_annotation_bed

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

REGION = "read1\t0\t6\tbSat\nread1\t6\t10\tarm\n"
REPEAT = "read1\t0\t10\tLINE\n"
# A composite overlay-style track (featureset name not in the hierarchy -> C2 skipped).
OVERLAY = (
    "read1\t0\t400\tcanonical_telomere\n"
    "read1\t400\t600\tbSat\n"
    "read1\t600\t1000\tcanonical_telomere\n"
)


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def _write(tmp_path: Path, **tracks: str) -> dict[str, Path]:
    out = {}
    for fs, text in tracks.items():
        p = tmp_path / f"{fs}.bed"
        p.write_text(text)
        out[fs] = p
    return out


def test_core_columns_and_values(h, tmp_path: Path):
    paths = _write(tmp_path, region=REGION, repeat=REPEAT)
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    m = core.build_feature_matrix(beds, h)  # window default 1000 > span -> broadcast

    row = m.rows["read1"]
    assert row["region__total_bp"] == 10
    assert row["region__frac__bSat"] == pytest.approx(0.6)
    assert row["region__bp__bSat"] == 6
    assert row["region__dmax__bSat"] == pytest.approx(0.6)
    assert row["region__max_block_bp__bSat"] == 6
    assert row["repeat__frac__LINE"] == pytest.approx(1.0)
    # thresholds: median([0.6])/3 = 0.2 -> clamp to max 0.05
    assert ("region", "bSat", pytest.approx(0.05)) in m.thresholds


def test_core_interspersion_over_composite(h, tmp_path: Path):
    paths = _write(tmp_path, overlay=OVERLAY)
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    m = core.build_feature_matrix(beds, h, interspersion_featureset="overlay")
    row = m.rows["read1"]
    # span 1000 bp -> per-kb == raw counts; 2 telomere<->satellite transitions
    assert row["interspersion__total"] == 2.0
    assert row["interspersion__tel_sat"] == 2.0


def test_core_unknown_interspersion_featureset(h, tmp_path: Path):
    paths = _write(tmp_path, region=REGION)
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    with pytest.raises(ValueError, match="interspersion featureset"):
        core.build_feature_matrix(beds, h, interspersion_featureset="nope")


def test_core_rejects_v1_feature(h, tmp_path: Path):
    paths = _write(tmp_path, region="read1\t0\t10\tarm_multigroup1\n")
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    with pytest.raises(ValueError, match="unknown feature"):
        core.build_feature_matrix(beds, h)


def test_cli(cli_runner, tmp_path: Path):
    paths = _write(tmp_path, region=REGION, repeat=REPEAT)
    out = tmp_path / "matrix.tsv"
    result = cli_runner.invoke(
        main,
        [
            "build-feature-matrix",
            "--bed",
            f"region={paths['region']}",
            "--bed",
            f"repeat={paths['repeat']}",
            "--hierarchy",
            str(HIERARCHY_TSV),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert out.is_file()
    sidecar = tmp_path / "matrix.adaptive_thresholds.tsv"
    assert sidecar.is_file()

    # Round-trip the matrix header + the single data row.
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "seq_id"
    values = dict(zip(header, lines[1].split("\t"), strict=True))
    assert values["seq_id"] == "read1"
    assert values["region__bp__bSat"] == "6"
    assert float(values["region__frac__bSat"]) == pytest.approx(0.6)

    # Sidecar has the threshold header + region/bSat row.
    sc = sidecar.read_text().splitlines()
    assert sc[0] == "featureset\tfeature\tthreshold"
    assert any(line.startswith("region\tbSat\t") for line in sc[1:])
