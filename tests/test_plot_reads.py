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
    (oriented,) = pr.orient_reads([read], {"canonical_telomere"})
    first = oriented.features[0]
    assert first[2] == "canonical_telomere"
    assert first[0] == 0  # telomere now at the top


def test_orient_leaves_telomere_top_read_untouched():
    read = pr.Read("S", "r", 100, [(0, 10, "canonical_telomere"), (10, 100, "aSat")])
    (oriented,) = pr.orient_reads([read], {"canonical_telomere"})
    assert list(oriented.features) == list(read.features)


def test_orient_fallback_used_when_no_primary():
    # No telomere, but a satellite at the end -> flip by the fallback set.
    read = pr.Read("S", "r", 100, [(0, 40, "rDNA"), (80, 100, "bSat")])
    (oriented,) = pr.orient_reads([read], {"canonical_telomere"}, fallback_features={"bSat"})
    assert oriented.features[0][2] == "bSat"
    assert oriented.features[0][0] == 0


def test_sort_reads_by_sample_then_length():
    reads = [
        pr.Read("B", "b1", 10, []),
        pr.Read("A", "a1", 5, []),
        pr.Read("A", "a2", 50, []),
    ]
    out = pr.sort_reads(reads, ["A", "B"])
    assert [r.read_id for r in out] == ["a2", "a1", "b1"]  # A first, longer first within A


# ----------------------------------------------------------------- legend
def test_filter_legend_features_cleans_underscores():
    used = {"chr13", "aSat", "noncanonical_telomere"}
    colors = {"chr13": "#111", "aSat": "#222", "noncanonical_telomere": "#333"}
    items = pr.filter_legend_features(used, colors)
    displays = [d for d, _c, _h in items]
    assert "chr13" in displays
    assert "noncanonical telomere" in displays  # underscores -> spaces
    # A feature with no color is omitted from the legend.
    assert pr.filter_legend_features({"mystery"}, colors) == []


def test_resolve_feature_color_strict():
    import pytest

    colors = {"aSat": "#fa0000"}
    assert pr.resolve_feature_color("aSat", colors) == "#fa0000"
    assert pr.resolve_feature_color("novel", colors) == "#ffffff"  # only novel may be absent
    with pytest.raises(pr.UnknownFeatureError):
        pr.resolve_feature_color("not_in_db", colors)


def test_render_errors_on_unknown_feature():
    import pytest

    reads = [pr.Read("S", "r", 100, [(0, 100, "totally_unknown")])]
    with pytest.raises(pr.UnknownFeatureError, match="totally_unknown"):
        pr.render(reads, {"aSat": "#fa0000"}, pr.PlotConfig(), ["S"])


def test_render_allows_novel_as_white():
    reads = [pr.Read("S", "r", 100, [(0, 50, "novel"), (50, 100, "aSat")])]
    svg = pr.render(reads, {"aSat": "#fa0000"}, pr.PlotConfig(feature_mode="raw"), ["S"])
    minidom.parseString(svg)
    assert "#fa0000" in svg.lower()


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


# ----------------------------------------------------------------- 3b: read-list / grouping
def _write_readlist(path: Path) -> None:
    path.write_text(
        "read_id\tgroup\tsubgroup\tcondition\n"
        "r1\tTumor\tA\thigh\n"
        "r2\tTumor\tB\tlow\n"
        "r3\tNormal\tA\thigh\n"
    )


def test_load_read_list_parses_header_and_rows(tmp_path: Path):
    rl = tmp_path / "list.tsv"
    _write_readlist(rl)
    read_ids, columns, read_data, read_rows = pr.load_read_list(rl)
    assert read_ids == {"r1", "r2", "r3"}
    assert columns == ["group", "subgroup", "condition"]
    assert read_data["r1"] == {"group": "Tumor", "subgroup": "A", "condition": "high"}
    assert len(read_rows) == 3


def test_build_and_apply_grouping():
    rows = [
        ("r1", {"group": "Tumor", "subgroup": "A"}),
        ("r2", {"group": "Tumor", "subgroup": "B"}),
        ("r3", {"group": "Normal", "subgroup": "A"}),
    ]
    read_groups, order = pr.build_grouping(rows, [("group", "Group"), ("subgroup", "Sub")])
    assert read_groups["r1"] == ("Tumor", "A")
    assert order == [("Tumor", "A"), ("Tumor", "B"), ("Normal", "A")]

    reads = [pr.Read("S", rid, 100, []) for rid in ("r1", "r2", "r3")]
    relabeled, sample_order, boundaries = pr.apply_grouping(reads, read_groups, order)
    assert sample_order == ["Tumor — A", "Tumor — B", "Normal — A"]
    assert relabeled[0].sample == "Tumor — A"
    # "Tumor — B" is the last subgroup before the group switches to Normal.
    assert "Tumor — B" in boundaries


def test_assign_heatmap_colors_uses_library_palette():
    from karyoplot.core.colors import TAB20

    read_metadata = {"r1": {"c": "x"}, "r2": {"c": "y"}, "r3": {"c": None}}
    cmap = pr.assign_heatmap_colors(["c"], read_metadata)
    assert cmap["c"]["x"] in TAB20 and cmap["c"]["y"] in TAB20
    assert cmap["c"]["x"] != cmap["c"]["y"]  # distinct values get distinct colors
    assert cmap["c"][None] == pr._HEATMAP_MISSING  # missing -> neutral gray


def test_process_read_list_defaults_and_filter(tmp_path: Path):
    rl = tmp_path / "list.tsv"
    _write_readlist(rl)
    reads = [pr.Read("S", rid, 100, [(0, 100, "aSat")]) for rid in ("r1", "r2", "r3", "r4")]
    data = pr.process_read_list(reads, rl, heatmap=True)
    assert {r.read_id for r in data.reads} == {"r1", "r2", "r3"}  # r4 not in list
    assert [c for c, _ in data.tier_specs] == ["group", "subgroup"]  # default first two cols
    assert data.metadata_columns == ["condition"]  # remaining col -> heatmap

    only_tumor = pr.process_read_list(reads, rl, heatmap=True, filter_groups=["Tumor"])
    assert {r.read_id for r in only_tumor.reads} == {"r1", "r2"}


def test_parse_markers(tmp_path: Path):
    m = tmp_path / "markers.tsv"
    m.write_text("read_id\tstart\tend\nr1\t100\t200\nr1\t300\t400\nr2\t50\t60\n")
    markers = pr.parse_markers(m)
    assert markers["r1"] == [(100, 200), (300, 400)]
    assert markers["r2"] == [(50, 60)]


def test_render_with_heatmap_and_grouping_well_formed():
    reads = [
        pr.Read("Tumor — A", "r1", 40000, [(0, 5000, "canonical_telomere"), (5000, 40000, "aSat")]),
        pr.Read("Normal — A", "r3", 18000, [(0, 2000, "canonical_telomere"), (2000, 18000, "bSat")]),
    ]
    order = ["Tumor — A", "Normal — A"]
    cfg = pr.PlotConfig(
        legend=True,
        label_tiers={"Tumor — A": ("Tumor", "A"), "Normal — A": ("Normal", "A")},
        group_subgroup_order=[("Tumor", "A"), ("Normal", "A")],
        tier_display_names={0: "Group", 1: "Sub"},
        metadata_columns=["condition"],
        read_metadata={"r1": {"condition": "high"}, "r3": {"condition": "high"}},
        heatmap_colors=pr.assign_heatmap_colors(
            ["condition"], {"r1": {"condition": "high"}, "r3": {"condition": "high"}}
        ),
        heatmap_display_names={"condition": "Condition"},
        markers={"r1": [(1000, 1200)]},
    )
    svg = pr.render(reads, _COLORS, cfg, order)
    minidom.parseString(svg)
    assert "Tumor" in svg and "Condition" in svg


def test_plot_reads_cli_heatmap_requires_read_list(cli_runner, tmp_path: Path):
    bed = tmp_path / "S.bed"
    _write_bed(bed)
    res = cli_runner.invoke(
        main,
        ["plot-reads", "--bed", f"S:{bed}", "--colors", str(COLORS_TSV), "--heatmap",
         "-o", str(tmp_path / "out.svg")],
    )
    assert res.exit_code != 0
    assert "require --read-list" in res.output


def test_plot_reads_cli_grouping_heatmap_markers(cli_runner, tmp_path: Path):
    bed = tmp_path / "S.bed"
    bed.write_text(
        "r1\t0\t5000\tcanonical_telomere\nr1\t5000\t40000\taSat\n"
        "r2\t0\t3000\tnoncanonical_telomere\nr2\t3000\t25000\tbSat\n"
        "r3\t0\t2000\tcanonical_telomere\nr3\t2000\t18000\tbSat\n"
    )
    rl = tmp_path / "list.tsv"
    _write_readlist(rl)
    markers = tmp_path / "markers.tsv"
    markers.write_text("r1\t1000\t1200\n")
    out = tmp_path / "out.svg"
    res = cli_runner.invoke(
        main,
        ["plot-reads", "--bed", f"S:{bed}", "--colors", str(COLORS_TSV),
         "--read-list", str(rl), "--heatmap-track", "condition:Condition",
         "--markers", str(markers), "--legend", "-o", str(out)],
    )
    assert res.exit_code == 0, res.output
    svg = out.read_text()
    minidom.parseString(svg)
    assert "Tumor" in svg and "Condition" in svg
