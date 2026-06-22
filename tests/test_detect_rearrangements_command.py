"""End-to-end test for the detect-rearrangements CLI."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main


def _write_overlay(path: Path, ab_reads: int, cd_reads: int) -> None:
    """A small overlay BED: `ab_reads` reads with A|B adjacent, `cd_reads` with C|D."""
    lines = []
    idx = 0
    for _ in range(ab_reads):
        lines += [f"r{idx}\t0\t100\tA", f"r{idx}\t100\t200\tB"]
        idx += 1
    for _ in range(cd_reads):
        lines += [f"r{idx}\t0\t100\tC", f"r{idx}\t100\t200\tD"]
        idx += 1
    path.write_text("\n".join(lines) + "\n")


def _read_tsv(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def test_detect_rearrangements_cli(cli_runner, tmp_path: Path):
    experiment = tmp_path / "exp.bed"
    control = tmp_path / "ctrl.bed"
    _write_overlay(experiment, ab_reads=20, cd_reads=80)  # (A,B) in 20/100
    _write_overlay(control, ab_reads=2, cd_reads=98)  # (A,B) in 2/100
    out = tmp_path / "calls.tsv"

    result = cli_runner.invoke(
        main,
        [
            "detect-rearrangements",
            "--experiment",
            str(experiment),
            "--control",
            str(control),
            "--window",
            "0",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert out.is_file()

    rows = _read_tsv(out)
    by_pair = {(r["feature_a"], r["feature_b"]): r for r in rows}

    # (A, B) is enriched in experiment -> a passing, enriched call.
    ab = by_pair[("A", "B")]
    assert ab["passes"] == "1"
    assert ab["direction"] == "enriched"
    assert float(ab["exp_rate"]) == 0.2
    assert float(ab["ctrl_rate"]) == 0.02

    # (C, D) is not differentially enriched -> not a call.
    assert by_pair[("C", "D")]["passes"] == "0"
