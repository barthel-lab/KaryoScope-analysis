"""End-to-end tests for overlay-annotations (core orchestration + CLI)."""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import overlay_annotations as core
from karyoscope_analysis.core.annotation_resolution import load_builtin_preset, load_spec
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import read_annotation_bed

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

# C4-valid tracks for one sequence, all spanning [0, 10).
REGION = "read1\t0\t5\tarm\nread1\t5\t10\tbSat\n"
REPEAT = "read1\t0\t3\tLINE\nread1\t3\t10\tnonrepeat\n"
SUBTEL = "read1\t0\t10\tnonsubtelomeric\n"


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def _write_tracks(tmp_path: Path, **tracks: str) -> dict[str, Path]:
    paths = {}
    for fs, text in tracks.items():
        p = tmp_path / f"{fs}.bed"
        p.write_text(text)
        paths[fs] = p
    return paths


def test_core_priority(h, tmp_path: Path):
    paths = _write_tracks(tmp_path, region=REGION, repeat=REPEAT, subtelomeric=SUBTEL)
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    spec = load_builtin_preset("priority", h)
    rows = list(core.overlay_annotations(beds, spec, h))
    assert rows == [
        ("read1", 0, 3, "LINE"),  # arm -> repeat (LINE)
        ("read1", 3, 5, "nonrepeat"),  # arm -> repeat (nonrepeat)
        ("read1", 5, 10, "bSat"),  # default region
    ]


def test_core_default_overlay_coalesces(h, tmp_path: Path):
    # Two tracks whose composite is constant across the sequence -> coalesced to one row.
    paths = _write_tracks(
        tmp_path,
        region="read1\t0\t10\tbSat\n",
        repeat="read1\t0\t4\tLINE\nread1\t4\t10\tLINE\n",
    )
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    spec = load_spec(
        {
            "name": "overlay",
            "precedence": ["region", "repeat"],
            "rules": [{"emit": {"composite": "all"}}],
        },
        h,
    )
    rows = list(core.overlay_annotations(beds, spec, h))
    assert rows == [("read1", 0, 10, "bSat:LINE")]


def test_core_rejects_mismatched_span(h, tmp_path: Path):
    paths = _write_tracks(
        tmp_path,
        region=REGION,
        repeat=REPEAT,
        subtelomeric="read1\t0\t8\tnonsubtelomeric\n",  # span 0-8 != 0-10
    )
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    spec = load_builtin_preset("priority", h)
    with pytest.raises(ValueError, match="read1"):
        list(core.overlay_annotations(beds, spec, h))


def test_core_rejects_v1_feature(h, tmp_path: Path):
    # arm_multigroup1 is a v1 name -> not in the v2 hierarchy -> C2 error.
    paths = _write_tracks(
        tmp_path,
        region="read1\t0\t10\tarm_multigroup1\n",
        repeat="read1\t0\t10\tLINE\n",
    )
    beds = {fs: read_annotation_bed(p) for fs, p in paths.items()}
    spec = load_spec(
        {
            "name": "o",
            "precedence": ["region", "repeat"],
            "rules": [{"emit": {"composite": "all"}}],
        },
        h,
    )
    with pytest.raises(ValueError, match="unknown feature"):
        list(core.overlay_annotations(beds, spec, h))


def test_cli_overlay(cli_runner, tmp_path: Path):
    paths = _write_tracks(tmp_path, region=REGION, repeat=REPEAT, subtelomeric=SUBTEL)
    out = tmp_path / "out.bed"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            "--bed",
            f"region={paths['region']}",
            "--bed",
            f"repeat={paths['repeat']}",
            "--bed",
            f"subtelomeric={paths['subtelomeric']}",
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--preset",
            "priority",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert read_annotation_bed(out)["read1"] == [
        (0, 3, "LINE"),
        (3, 5, "nonrepeat"),
        (5, 10, "bSat"),
    ]


def test_cli_preset_and_spec_conflict(cli_runner, tmp_path: Path):
    paths = _write_tracks(tmp_path, region=REGION, repeat=REPEAT)
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            "--bed",
            f"region={paths['region']}",
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--preset",
            "priority",
            "--spec",
            str(paths["region"]),
            "-o",
            str(tmp_path / "out.bed"),
        ],
    )
    assert result.exit_code != 0
    assert "mutually exclusive" in result.output


def test_cli_overlay_threads_match_serial(cli_runner, tmp_path: Path):
    # Many sequences with C4 tracks: the parallel overlay must be byte-identical to serial
    # (per-sequence resolution + order-preserving imap), which is what makes --threads safe.
    import random

    rng = random.Random(0)
    region_lines, repeat_lines = [], []
    for r in range(40):
        rid = f"read{r:03d}"
        c1, c2 = rng.randint(2, 8), rng.randint(2, 8)
        region_lines += [f"{rid}\t0\t{c1}\tarm", f"{rid}\t{c1}\t10\tbSat"]
        repeat_lines += [f"{rid}\t0\t{c2}\tLINE", f"{rid}\t{c2}\t10\tnonrepeat"]
    paths = _write_tracks(
        tmp_path,
        region="\n".join(region_lines) + "\n",
        repeat="\n".join(repeat_lines) + "\n",
    )

    def run(out: Path, threads: int) -> str:
        res = cli_runner.invoke(
            main,
            [
                "overlay-annotations",
                "--bed",
                f"region={paths['region']}",
                "--bed",
                f"repeat={paths['repeat']}",
                "--hierarchy",
                str(HIERARCHY_TSV),
                "--threads",
                str(threads),
                "-o",
                str(out),
            ],
        )
        assert res.exit_code == 0, res.output
        return out.read_text()

    serial = run(tmp_path / "serial.bed", 1)
    parallel = run(tmp_path / "parallel.bed", 4)
    assert parallel == serial
    assert serial.strip(), "expected non-empty overlay output"


def test_cli_overlay_rejects_negative_threads(cli_runner, tmp_path: Path):
    paths = _write_tracks(tmp_path, region=REGION, repeat=REPEAT, subtelomeric=SUBTEL)
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            "--bed",
            f"region={paths['region']}",
            "--bed",
            f"repeat={paths['repeat']}",
            "--bed",
            f"subtelomeric={paths['subtelomeric']}",
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--threads",
            "-3",
            "-o",
            str(tmp_path / "out.bed"),
        ],
    )
    assert result.exit_code != 0
    assert "--threads" in result.output
