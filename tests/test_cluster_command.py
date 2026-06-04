"""End-to-end test for the cluster CLI (Engine B)."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


def _read_tsv(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def test_cluster_cli(cli_runner, tmp_path: Path):
    # x and y dovetail (share the B block); z is unrelated -> two clusters.
    bed = tmp_path / "overlay.bed"
    bed.write_text(
        "x\t0\t1000\tarm\nx\t1000\t2000\tbSat\n"
        "y\t0\t1000\tbSat\ny\t1000\t2000\tgSat\n"
        "z\t0\t1000\tchr1\nz\t1000\t2000\tchr2\n"
    )
    out = tmp_path / "clusters.tsv"
    result = cli_runner.invoke(
        main,
        [
            "cluster",
            "--input",
            str(bed),
            "--hierarchy",
            str(HIERARCHY_TSV),
            # default repeat-mask: bSat/arm/gSat are structural (not interspersed repeats),
            # so they keep full weight and the 1000 bp bSat overlap clears the threshold.
            "--min-overlap-bp",
            "500",
            "--min-identity",
            "0.9",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output

    clusters = _read_tsv(out)
    sizes = sorted(int(c["size"]) for c in clusters)
    assert sizes == [1, 2]  # {x,y} and {z}

    consensus = tmp_path / "clusters.consensus.bed"
    layout = tmp_path / "clusters.layout.tsv"
    assert consensus.is_file()
    assert layout.is_file()

    # The 2-read cluster's seed and members appear in the layout.
    layout_rows = _read_tsv(layout)
    big = next(c for c in clusters if int(c["size"]) == 2)
    members = {r["read_id"] for r in layout_rows if r["cluster_id"] == big["cluster_id"]}
    assert members == {"x", "y"}
    seeds = {r["read_id"] for r in layout_rows if r["is_seed"] == "1"}
    assert big["seed"] in seeds


def test_cluster_cli_min_length_filters_reads(cli_runner, tmp_path: Path):
    bed = tmp_path / "overlay.bed"
    bed.write_text("short\t0\t100\tarm\nlong\t0\t5000\tarm\nlong\t5000\t9000\tbSat\n")
    out = tmp_path / "c.tsv"
    result = cli_runner.invoke(
        main,
        [
            "cluster",
            "--input",
            str(bed),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--min-length",
            "1000",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    layout_rows = _read_tsv(tmp_path / "c.layout.tsv")
    assert {r["read_id"] for r in layout_rows} == {"long"}  # short read filtered out


def test_cluster_cli_repeat_mask_breaks_repeat_only_overlap(cli_runner, tmp_path: Path):
    # x and y dovetail ONLY through a shared LINE block (a genome-wide interspersed repeat);
    # the flanks (chr1/chr2) are unrelated, so the LINE block is the only positive overlap.
    bed = tmp_path / "overlay.bed"
    bed.write_text("x\t0\t2000\tchr1\nx\t2000\t4000\tLINE\ny\t0\t2000\tLINE\ny\t2000\t4000\tchr2\n")
    common = [
        "cluster",
        "--input",
        str(bed),
        "--hierarchy",
        str(HIERARCHY_TSV),
        "--min-overlap-bp",
        "1000",
        "--min-identity",
        "0.9",
    ]
    # default repeat-mask: LINE is masked -> the LINE-only overlap carries no weight -> 2 singletons.
    masked = cli_runner.invoke(main, [*common, "-o", str(tmp_path / "masked.tsv")])
    assert masked.exit_code == 0, masked.output
    assert sorted(int(c["size"]) for c in _read_tsv(tmp_path / "masked.tsv")) == [1, 1]
    # uniform: LINE counts -> x,y dovetail through it -> one 2-read cluster.
    uni = cli_runner.invoke(
        main, [*common, "--weight-method", "uniform", "-o", str(tmp_path / "uni.tsv")]
    )
    assert uni.exit_code == 0, uni.output
    assert sorted(int(c["size"]) for c in _read_tsv(tmp_path / "uni.tsv")) == [2]
