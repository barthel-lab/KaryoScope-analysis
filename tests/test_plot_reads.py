"""Tests for the plot-reads SVG core (Phase 3a) and CLI command."""

from __future__ import annotations

import xml.dom.minidom as minidom
from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import plot_reads as pr

COLORS_TSV = Path(__file__).resolve().parent / "data" / "colors.tsv"

_COLORS = {
    "canonical_telomere": "#008000",
    "noncanonical_telomere": "#0000ff",
    "aSat": "#fa0000",
    "bSat": "#0081fa",
}


def _write_bed(path: Path) -> None:
    path.write_text(
        "r1\t0\t5000\tcanonical_telomere\n"
        "r1\t5000\t40000\taSat\n"
        "r2\t0\t3000\tnoncanonical_telomere\n"
        "r2\t3000\t25000\tbSat\n"
    )


# ----------------------------------------------------------------- loading
def test_parse_bed_and_load_specs(tmp_path: Path):
    bed = tmp_path / "HeLa.bed"
    _write_bed(bed)
    reads, order = pr.load_bed_specs([f"HeLa:{bed}"])
    assert order == ["HeLa"]
    assert {r.read_id for r in reads} == {"r1", "r2"}
    r1 = next(r for r in reads if r.read_id == "r1")
    assert r1.length == 40000  # max feature end
    assert r1.sample == "HeLa"


def test_load_specs_missing_file_raises(tmp_path: Path):
    import pytest

    with pytest.raises(FileNotFoundError):
        pr.load_bed_specs([str(tmp_path / "nope.bed")])


# ----------------------------------------------------------------- orientation
def test_orient_telomere_flips_read_with_telomere_at_end():
    # Telomere sits at the end (avg position past midpoint) -> read is flipped.
    read = pr.Read("S", "r", 100, [(0, 50, "aSat"), (90, 100, "canonical_telomere")])
    (oriented,) = pr.orient_reads([read], "telomere")
    first = oriented.features[0]
    assert first[2] == "canonical_telomere"
    assert first[0] == 0  # telomere now at the top


def test_orient_leaves_telomere_top_read_untouched():
    read = pr.Read("S", "r", 100, [(0, 10, "canonical_telomere"), (10, 100, "aSat")])
    (oriented,) = pr.orient_reads([read], "telomere")
    assert list(oriented.features) == list(read.features)


def test_sort_reads_by_sample_then_length():
    reads = [
        pr.Read("B", "b1", 10, []),
        pr.Read("A", "a1", 5, []),
        pr.Read("A", "a2", 50, []),
    ]
    out = pr.sort_reads(reads, ["A", "B"])
    assert [r.read_id for r in out] == ["a2", "a1", "b1"]  # A first, longer first within A


# ----------------------------------------------------------------- legend
def test_filter_legend_features_dedups_specific_and_cleans():
    used = {"chr13_specific", "chr13", "aSat", "noncanonical_telomere"}
    colors = {"chr13": "#111", "aSat": "#222", "noncanonical_telomere": "#333"}
    items = pr.filter_legend_features(used, colors)
    displays = [d for d, _c, _h in items]
    assert "chr13" in displays
    assert displays.count("chr13") == 1  # _specific variant dedup'd
    assert "noncanonical telomere" in displays  # underscores -> spaces


def test_compute_legend_layout_sections_get_distinct_columns():
    items = [
        ("a", "#1", "Telomere"),
        ("b", "#2", "Telomere"),
        ("c", "#3", "Satellite"),
    ]
    layout = pr.compute_legend_layout(items, max_width=1000, font_size=12)
    headers = [it for it in layout.items if it.is_header]
    assert {h.label for h in headers} == {"Telomere", "Satellite"}
    # Each section header starts its own column -> distinct x for the two headers.
    assert len({h.x for h in headers}) == 2


# ----------------------------------------------------------------- rendering
def test_render_vertical_and_horizontal_well_formed():
    reads = [
        pr.Read("HeLa", "r1", 40000, [(0, 5000, "canonical_telomere"), (5000, 40000, "aSat")]),
        pr.Read("HeLa", "r2", 25000, [(0, 3000, "noncanonical_telomere"), (3000, 25000, "bSat")]),
    ]
    cfg = pr.PlotConfig(legend=True, feature_mode="smooth")
    for horizontal in (False, True):
        svg = pr.render(reads, _COLORS, cfg, ["HeLa"], horizontal=horizontal)
        minidom.parseString(svg)  # raises if malformed
        assert svg.startswith("<svg")
        assert "#fa0000" in svg.lower()  # DB feature color used


def test_render_empty_raises():
    import pytest

    with pytest.raises(ValueError):
        pr.render([], _COLORS, pr.PlotConfig(), [])


# ----------------------------------------------------------------- CLI
def test_plot_reads_cli_vertical(cli_runner, tmp_path: Path):
    bed = tmp_path / "HeLa.bed"
    _write_bed(bed)
    out = tmp_path / "out.svg"
    res = cli_runner.invoke(
        main,
        ["plot-reads", "--bed", f"HeLa:{bed}", "--colors", str(COLORS_TSV),
         "--orient", "telomere", "--legend", "-o", str(out)],
    )
    assert res.exit_code == 0, res.output
    assert out.exists()
    minidom.parseString(out.read_text())


def test_plot_reads_cli_length_filter_errors(cli_runner, tmp_path: Path):
    bed = tmp_path / "HeLa.bed"
    _write_bed(bed)
    out = tmp_path / "out.svg"
    res = cli_runner.invoke(
        main,
        ["plot-reads", "--bed", f"HeLa:{bed}", "--colors", str(COLORS_TSV),
         "--min-length", "999999999", "-o", str(out)],
    )
    assert res.exit_code != 0
    assert "no reads to plot" in res.output
    assert not out.exists()
