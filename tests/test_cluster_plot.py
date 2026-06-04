"""Tests for the colors reader, the SVG cluster renderer, and the cluster-plot CLI."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import cluster_plot as render
from karyoscope_analysis.core.io.colors import load_colors

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"
COLORS_TSV = Path(__file__).resolve().parent / "data" / "colors.tsv"


# ----------------------------------------------------------------- colors reader
def test_load_colors():
    colors = load_colors(COLORS_TSV)
    assert "feature_set" not in colors  # header skipped
    assert colors["aSat"].startswith("#")  # a region satellite has a color
    assert all(c.startswith("#") for c in colors.values())


# ----------------------------------------------------------------- renderer
def test_structural_feature_and_color():
    colors = {"aSat": "#111111"}
    assert render.structural_feature("chr13:aSat") == "aSat"
    assert render.structural_feature("aSat") == "aSat"
    assert render.feature_color("chr13:aSat", colors) == "#111111"  # by structural layer
    # unknown feature -> a stable auto-palette color
    auto = render.feature_color("mystery", {})
    assert auto in render._AUTO_PALETTE
    assert render.feature_color("mystery", {}) == auto  # deterministic


def test_render_cluster_svg():
    placed = [
        render.PlacedRead("seed", True, False, 0, [(0, 100, "aSat"), (100, 200, "bSat")]),
        render.PlacedRead("m", False, False, 100, [(0, 100, "bSat"), (100, 200, "HSat3")]),
    ]
    consensus = [(0, 100, "aSat"), (100, 200, "bSat")]
    colors = {"aSat": "#111111", "bSat": "#222222"}
    svg = render.render_cluster_svg(placed, consensus, colors, title="cluster_0")
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "#111111" in svg and "#222222" in svg  # mapped colors used
    assert "cluster_0" in svg and "(seed)" in svg
    assert svg.count("<rect") >= 6  # consensus(2) + seed(2) + member(2) + legend swatches
    assert "HSat3" in svg  # auto-colored feature appears in the legend


def test_render_clusters_svg_stacks_panels():
    panels = [
        render.ClusterPanel(
            "cluster_0  n=1  chr13",
            [render.PlacedRead("a", True, False, 0, [(0, 100, "chr13:aSat")])],
            [(0, 100, "chr13:aSat")],
        ),
        render.ClusterPanel(
            "cluster_1  n=1  chr21",
            [render.PlacedRead("b", True, False, 0, [(0, 100, "chr21:bSat")])],
            [(0, 100, "chr21:bSat")],
        ),
    ]
    svg = render.render_clusters_svg(panels, {})
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "cluster_0  n=1  chr13" in svg and "cluster_1  n=1  chr21" in svg  # both panels
    assert "aSat" in svg and "bSat" in svg  # one shared legend (structural layer)


def test_render_reversed_member_reflects_coordinates():
    # a reversed 2-segment read: oriented order is [second, first].
    placed = [render.PlacedRead("r", False, True, 0, [(0, 30, "A"), (30, 100, "B")])]
    oriented = render._oriented(placed[0])
    assert [f for _, _, f in oriented] == ["B", "A"]  # order reversed
    assert oriented[0] == (0, 70, "B") and oriented[1] == (70, 100, "A")  # lengths preserved


# ----------------------------------------------------------------- CLI end-to-end
def test_cluster_plot_cli(cli_runner, tmp_path: Path):
    overlay = tmp_path / "overlay.bed"
    overlay.write_text(
        "x\t0\t1000\taSat\nx\t1000\t2000\tbSat\ny\t0\t1000\tbSat\ny\t1000\t2000\tHSat3\n"
    )
    clusters = tmp_path / "clusters.tsv"
    res = cli_runner.invoke(
        main,
        [
            "cluster",
            "--input",
            str(overlay),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--min-overlap-bp",
            "1000",
            "--min-identity",
            "0.9",
            "-o",
            str(clusters),
        ],
    )
    assert res.exit_code == 0, res.output

    svg_out = tmp_path / "plot.svg"
    res2 = cli_runner.invoke(
        main,
        [
            "cluster-plot",
            "--layout",
            str(tmp_path / "clusters.layout.tsv"),
            "--consensus",
            str(tmp_path / "clusters.consensus.bed"),
            "--overlay",
            str(overlay),
            "--colors",
            str(COLORS_TSV),
            "--cluster-id",
            "cluster_0",
            "-o",
            str(svg_out),
        ],
    )
    assert res2.exit_code == 0, res2.output
    svg = svg_out.read_text()
    assert svg.startswith("<svg") and "<rect" in svg
    assert "x" in svg and "y" in svg  # both reads labeled


def test_cluster_plot_cli_all_clusters(cli_runner, tmp_path: Path):
    # Drive the combined mode directly from hand-written layout/consensus (two clusters).
    overlay = tmp_path / "overlay.bed"
    overlay.write_text(
        "a\t0\t1000\taSat\na\t1000\t2000\tbSat\n"
        "b\t0\t1000\tbSat\nb\t1000\t2000\tgSat\n"
        "c\t0\t1000\tct\nc\t1000\t2000\trDNA\n"
        "d\t0\t1000\trDNA\nd\t1000\t2000\tarm\n"
    )
    layout = tmp_path / "c.layout.tsv"
    layout.write_text(
        "cluster_id\tread_id\tis_seed\treversed\toffset\tlength\n"
        "cluster_0\ta\t1\t0\t0\t2000\ncluster_0\tb\t0\t0\t1000\t2000\n"
        "cluster_1\tc\t1\t0\t0\t2000\ncluster_1\td\t0\t0\t1000\t2000\n"
    )
    consensus = tmp_path / "c.consensus.bed"
    consensus.write_text(
        "cluster_id\tstart\tend\tfeature\tsupport\tcoverage\n"
        "cluster_0\t0\t1000\taSat\t1\t1\ncluster_0\t1000\t2000\tbSat\t2\t2\n"
        "cluster_1\t0\t1000\tct\t1\t1\ncluster_1\t1000\t2000\trDNA\t2\t2\n"
    )
    out = tmp_path / "all.svg"
    res = cli_runner.invoke(
        main,
        [
            "cluster-plot",
            "--layout",
            str(layout),
            "--consensus",
            str(consensus),
            "--overlay",
            str(overlay),
            "--colors",
            str(COLORS_TSV),
            "--min-cluster-size",
            "2",
            "-o",
            str(out),
        ],  # no --cluster-id -> all clusters
    )
    assert res.exit_code == 0, res.output
    assert "2 cluster" in res.output  # both clusters rendered
    svg = out.read_text()
    assert svg.startswith("<svg")
    assert "cluster_0" in svg and "cluster_1" in svg  # two panels in one SVG
    assert "aSat" in svg and "rDNA" in svg  # features from both clusters in the shared legend
