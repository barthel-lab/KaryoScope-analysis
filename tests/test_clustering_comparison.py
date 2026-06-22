"""Tests for clustering comparison (ARI/NMI + overlap) and the compare-clusterings command."""

from __future__ import annotations

import math
from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import clustering_comparison as cc


def test_identical_partitions_score_one():
    labels = {"r1": "a", "r2": "a", "r3": "b", "r4": "b"}
    res = cc.compare(labels, dict(labels))
    assert math.isclose(res.ari, 1.0) and math.isclose(res.nmi, 1.0)
    assert res.n_common == 4 and res.n_clusters_1 == 2 and res.n_clusters_2 == 2


def test_relabeled_partitions_still_score_one():
    # Same partition, different cluster names -> ARI/NMI are label-invariant -> still 1.
    a = {"r1": "a", "r2": "a", "r3": "b", "r4": "b"}
    b = {"r1": "X", "r2": "X", "r3": "Y", "r4": "Y"}
    res = cc.compare(a, b)
    assert math.isclose(res.ari, 1.0)


def test_independent_partitions_low_ari():
    a = {"r1": "a", "r2": "a", "r3": "b", "r4": "b"}
    b = {"r1": "x", "r2": "y", "r3": "x", "r4": "y"}  # orthogonal split
    res = cc.compare(a, b)
    assert res.ari < 0.5


def test_compare_no_common_reads_raises():
    import pytest

    with pytest.raises(ValueError, match="share no read"):
        cc.compare({"r1": "a"}, {"r2": "b"})


def test_overlap_pairs_sorted_desc():
    a = {"r1": "a", "r2": "a", "r3": "b"}
    b = {"r1": "x", "r2": "x", "r3": "y"}
    pairs = cc.overlap_pairs(a, b)
    assert pairs[0] == ("a", "x", 2)
    assert ("b", "y", 1) in pairs


def _layout(path: Path, assignments: dict[str, str]) -> None:
    rows = ["cluster_id\tread_id\tis_seed\treversed\tstart\tend\tfeature"]
    rows += [f"{c}\t{r}\t1\t0\t0\t100\taSat" for r, c in assignments.items()]
    path.write_text("\n".join(rows) + "\n")


def test_compare_clusterings_cli(cli_runner, tmp_path: Path):
    a = {"r1": "c0", "r2": "c0", "r3": "c1", "r4": "c1"}
    b = {"r1": "k0", "r2": "k0", "r3": "k1", "r4": "k1"}  # identical partition, renamed
    l1, l2 = tmp_path / "a.layout.tsv", tmp_path / "b.layout.tsv"
    _layout(l1, a)
    _layout(l2, b)
    out = tmp_path / "cmp.txt"
    res = cli_runner.invoke(
        main,
        [
            "compare-clusterings",
            "--layout1",
            str(l1),
            "--layout2",
            str(l2),
            "--label1",
            "run1",
            "--label2",
            "run2",
            "-o",
            str(out),
        ],
    )
    assert res.exit_code == 0, res.output
    assert "ARI): 1.0" in out.read_text()
    overlap = (tmp_path / "cmp.overlap.tsv").read_text()
    assert overlap.splitlines()[0] == "cluster_run1\tcluster_run2\tn_shared"
    assert "c0\tk0\t2" in overlap
