"""Tests for the enrichment heatmap (core selection + the plot-enrichment command)."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import enrichment_plot as ep


def _rows():
    return [
        {
            "cluster_id": "c1",
            "n_total": "10",
            "enriched": "1",
            "log2fc_HeLa": "-inf",
            "log2fc_U2OS": "2.5",
        },
        {
            "cluster_id": "c2",
            "n_total": "5",
            "enriched": "1",
            "log2fc_HeLa": "1.2",
            "log2fc_U2OS": "-inf",
        },
        {
            "cluster_id": "c3",
            "n_total": "3",
            "enriched": "0",  # not enriched
            "log2fc_HeLa": "0.1",
            "log2fc_U2OS": "0.1",
        },
    ]


def test_select_rows_enriched_only_sorted_by_strongest():
    rows = ep.select_rows(_rows(), ["HeLa", "U2OS"], {"c1": "ECTR"})
    assert [r.cluster_id for r in rows] == ["c1", "c2"]  # c3 excluded; c1 (2.5) before c2 (1.2)
    assert rows[0].label == "ECTR"
    assert rows[0].log2fc["HeLa"] == float("-inf")


def test_select_rows_all_clusters_includes_unenriched():
    rows = ep.select_rows(_rows(), ["HeLa", "U2OS"], {}, enriched_only=False)
    assert {r.cluster_id for r in rows} == {"c1", "c2", "c3"}


def test_select_rows_max_clusters():
    rows = ep.select_rows(_rows(), ["HeLa", "U2OS"], {}, max_clusters=1)
    assert [r.cluster_id for r in rows] == ["c1"]


def test_render_heatmap_writes_file(tmp_path: Path):
    rows = ep.select_rows(_rows(), ["HeLa", "U2OS"], {"c1": "ECTR"})
    out = tmp_path / "heat.png"
    ep.render_heatmap(rows, ["HeLa", "U2OS"], str(out))
    assert out.exists() and out.stat().st_size > 0


def test_orient_to_breakpoint_flips_and_finds_breakpoint():
    tel = {"canonical_telomere"}
    # Telomere on the RIGHT -> mirrored to the left; breakpoint = leading telomere block end.
    segs = [(0, 8000, "chr4:q_arm"), (8000, 10000, "chr4:canonical_telomere")]
    nsegs, bp = ep._orient_to_breakpoint(segs, tel)
    assert nsegs[0][2].endswith("canonical_telomere")  # telomere now first (left)
    assert bp == 2000  # 2 kb telomere block at the left


def test_orient_to_breakpoint_no_telomere_breakpoint_zero():
    _nsegs, bp = ep._orient_to_breakpoint([(0, 5000, "chr4:aSat")], {"canonical_telomere"})
    assert bp == 0.0  # no leading telomere -> aligns at start


def test_render_heatmap_with_consensus_panel(tmp_path: Path):
    rows = ep.select_rows(_rows(), ["HeLa", "U2OS"], {"c1": "ECTR", "c2": "subtelomere"})
    consensus = {
        "c1": [(0, 2000, "chr4:canonical_telomere"), (2000, 5000, "chr4:q_arm")],
        "c2": [(0, 3000, "chr7:aSat"), (3000, 4000, "chr7:novel")],  # novel -> white, no crash
    }
    colors = {"canonical_telomere": "#008000", "q_arm": "#808080", "aSat": "#fa0000"}
    out = tmp_path / "heat_consensus.png"
    ep.render_heatmap(rows, ["HeLa", "U2OS"], str(out), consensus=consensus, colors=colors)
    assert out.exists() and out.stat().st_size > 0


def test_plot_enrichment_cli(cli_runner, tmp_path: Path):
    enr = tmp_path / "enrichment.tsv"
    enr.write_text(
        "cluster_id\tn_total\tn_HeLa\tn_U2OS\tlog2fc_HeLa\tlog2fc_U2OS\ttop_group\tprivate\tenriched\n"
        "c1\t10\t0\t10\t-inf\t2.5\tU2OS\t1\t1\n"
        "c2\t5\t4\t1\t1.2\t-0.3\tHeLa\t0\t1\n"
    )
    annot = tmp_path / "annot.tsv"
    annot.write_text(
        "cluster_id\tsize\twidth\tlabel\tchromosomes\tconsensus_signature\n"
        "c1\t10\t1000\tECTR\tchr4+chr18\tcanonical_telomere > q_arm\n"
    )
    out = tmp_path / "heat.png"
    res = cli_runner.invoke(
        main,
        ["plot-enrichment", "--enrichment", str(enr), "--annot", str(annot), "-o", str(out)],
    )
    assert res.exit_code == 0, res.output
    assert out.exists() and out.stat().st_size > 0
    assert "2 clusters x 2 groups" in res.output
