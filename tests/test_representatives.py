"""Tests for the consensus-as-representative catalog and the select-representatives command."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import representatives as rep


def test_consensus_signature_collapses_consecutive_dupes():
    segs = [
        (0, 5, "canonical_telomere"),
        (5, 9, "canonical_telomere"),  # consecutive dupe -> collapsed
        (9, 20, "ITS"),
        (20, 40, "bSat"),
    ]
    assert rep.consensus_signature(segs) == "canonical_telomere > ITS > bSat"


def test_consensus_signature_sorts_by_start():
    segs = [(20, 40, "bSat"), (0, 20, "aSat")]  # given out of order
    assert rep.consensus_signature(segs) == "aSat > bSat"


def test_build_catalog_filters_and_sorts():
    sizes = {"c0": 5, "c1": 1, "c2": 3}
    widths = {"c0": 100, "c1": 50, "c2": 80}
    consensus = {
        "c0": [(0, 50, "aSat"), (50, 100, "bSat")],
        "c1": [(0, 50, "ITS")],
        "c2": [(0, 80, "canonical_telomere")],
    }
    reps = rep.build_catalog(sizes, widths, consensus, min_cluster_size=2)
    # c1 (size 1) dropped; sorted by descending size -> c0 then c2.
    assert [r.cluster_id for r in reps] == ["c0", "c2"]
    assert reps[0].size == 5 and reps[0].n_segments == 2 and reps[0].width == 100
    assert reps[0].signature == "aSat > bSat"


def test_select_representatives_cli(cli_runner, tmp_path: Path):
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text(
        "cluster_id\tsize\tseed\tconsensus_segments\twidth\torientation_conflict\n"
        "cluster_0\t4\ts0\t2\t200\t0\n"
        "cluster_1\t1\ts1\t1\t100\t0\n"  # singleton -> excluded at default min size 2
    )
    consensus = tmp_path / "consensus.bed"
    consensus.write_text(
        "cluster_id\tstart\tend\tfeature\tsupport\tcoverage\n"
        "cluster_0\t0\t100\tcanonical_telomere\t4\t4\n"
        "cluster_0\t100\t200\taSat\t3\t4\n"
        "cluster_1\t0\t100\tbSat\t1\t1\n"
    )
    out = tmp_path / "reps.tsv"
    res = cli_runner.invoke(
        main,
        ["select-representatives", "--clusters", str(clusters), "--consensus", str(consensus),
         "-o", str(out)],
    )
    assert res.exit_code == 0, res.output
    lines = out.read_text().splitlines()
    assert lines[0] == "cluster_id\tsize\tn_segments\twidth\tconsensus_signature"
    assert len(lines) == 2  # header + cluster_0 only (singleton cluster_1 excluded)
    assert lines[1].startswith("cluster_0\t4\t2\t200\t")
    assert "canonical_telomere > aSat" in lines[1]
