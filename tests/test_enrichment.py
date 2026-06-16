"""Tests for the cluster-enrichment model and the test-enrichment / pool-samples commands."""

from __future__ import annotations

import gzip
import math
from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import enrichment as enrich


def _example():
    """A small pooled clustering: group A has 4 reads, B has 16 (depth asymmetry on purpose)."""
    read_to_group = {f"a{i}": "A" for i in range(1, 5)} | {f"b{i}": "B" for i in range(1, 17)}
    read_to_cluster = {}
    for r in ("a1", "a2", "a4"):
        read_to_cluster[r] = "k_altlike"  # private to A (3 of A's 4 reads)
    for r in ("a3", "b1", "b2", "b3"):
        read_to_cluster[r] = "k_bal"  # mixed
    for i in range(4, 17):
        read_to_cluster[f"b{i}"] = "k_b"  # private to B (13 of B's 16 reads)
    return read_to_cluster, read_to_group


def test_compute_enrichment_private_and_effect():
    read_to_cluster, read_to_group = _example()
    results, groups, totals = enrich.compute_enrichment(read_to_cluster, read_to_group)
    assert groups == ["A", "B"]
    assert totals == {"A": 4, "B": 16}
    by_id = {r.cluster_id: r for r in results}

    alt = by_id["k_altlike"]
    assert alt.private and alt.top_group == "A" and alt.enriched
    # 3 of A's 4 reads (0.75) vs pooled 3/20 (0.15) -> log2(5) ~= 2.32, well over the 2x threshold.
    assert math.isclose(alt.effects["A"], math.log2(0.75 / 0.15), rel_tol=1e-6)
    assert alt.counts == {"A": 3, "B": 0}

    # A B-private cluster is NOT strongly enriched here, because B is the large (deep) group:
    # 13/16 (0.81) vs pooled 13/20 (0.65) -> log2(1.25) ~= 0.32 < 1. Depth asymmetry handled.
    b = by_id["k_b"]
    assert b.private and b.top_group == "B" and not b.enriched

    bal = by_id["k_bal"]
    assert not bal.private and not bal.enriched


def test_singletons_are_not_enriched():
    # A single-read cluster has a large fold-change but no compositional evidence -> not enriched.
    read_to_cluster = {"u1": "solo"}
    read_to_group = {"u1": "U2OS", "h1": "HeLa", "h2": "HeLa"}  # h1/h2 unclustered here
    read_to_cluster |= {"h1": "k", "h2": "k"}
    results, _g, _t = enrich.compute_enrichment(read_to_cluster, read_to_group, min_cluster_size=2)
    solo = next(r for r in results if r.cluster_id == "solo")
    assert solo.private and not solo.enriched  # private, but a singleton -> not called enriched
    assert solo.effects["U2OS"] > 1  # the fold-change is large; the size gate is what stops it
    # min_cluster_size=1 would let it through.
    res1, _g, _t = enrich.compute_enrichment(read_to_cluster, read_to_group, min_cluster_size=1)
    assert next(r for r in res1 if r.cluster_id == "solo").enriched


def test_compute_enrichment_results_sorted_by_top_effect():
    read_to_cluster, read_to_group = _example()
    results, _g, _t = enrich.compute_enrichment(read_to_cluster, read_to_group)
    top_effects = [r.effects[r.top_group] for r in results]
    assert top_effects == sorted(top_effects, reverse=True)  # most-enriched first


def test_enrichment_tsv_shape():
    read_to_cluster, read_to_group = _example()
    results, groups, _t = enrich.compute_enrichment(read_to_cluster, read_to_group)
    tsv = enrich.enrichment_tsv(results, groups)
    header = tsv.splitlines()[0].split("\t")
    assert header == [
        "cluster_id", "n_total", "n_A", "n_B", "frac_A", "frac_B",
        "log2fc_A", "log2fc_B", "top_group", "private", "enriched",
    ]
    assert len(tsv.splitlines()) == 1 + len(results)


# ----------------------------------------------------------------- CLI: pool-samples + test-enrichment
def _write_bed(path: Path, read_prefix: str, feature: str, gz: bool = False) -> None:
    body = f"{read_prefix}_r1\t0\t1000\t{feature}\n{read_prefix}_r2\t0\t1000\t{feature}\n"
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(body)
    else:
        path.write_text(body)


def test_pool_samples_namespaces_and_writes_read_list(cli_runner, tmp_path: Path):
    a, b = tmp_path / "A.bed", tmp_path / "B.bed.gz"
    _write_bed(a, "x", "aSat")
    _write_bed(b, "x", "bSat", gz=True)  # same read names as A -> must be namespaced apart
    pooled = tmp_path / "pooled.bed"
    res = cli_runner.invoke(
        main, ["pool-samples", "--bed", f"A:{a}", "--bed", f"B:{b}", "-o", str(pooled)]
    )
    assert res.exit_code == 0, res.output
    rows = pooled.read_text().splitlines()
    ids = {r.split("\t")[0] for r in rows}
    assert ids == {"A|x_r1", "A|x_r2", "B|x_r1", "B|x_r2"}  # namespaced, no collision
    rl = (tmp_path / "pooled.samples.tsv").read_text().splitlines()
    assert rl[0] == "read_id\tsample"
    assert "A|x_r1\tA" in rl and "B|x_r2\tB" in rl


def test_test_enrichment_cli(cli_runner, tmp_path: Path):
    layout = tmp_path / "c.layout.tsv"
    layout.write_text(
        "cluster_id\tread_id\tis_seed\treversed\tstart\tend\tfeature\n"
        "k0\tU2OS|r1\t1\t0\t0\t100\tITS\n"
        "k0\tU2OS|r2\t0\t0\t0\t100\tITS\n"
        "k1\tHeLa|r3\t1\t0\t0\t100\taSat\n"
    )
    read_list = tmp_path / "samples.tsv"
    read_list.write_text(
        "read_id\tsample\nU2OS|r1\tU2OS\nU2OS|r2\tU2OS\nHeLa|r3\tHeLa\n"
    )
    out = tmp_path / "enrichment.tsv"
    res = cli_runner.invoke(
        main,
        ["test-enrichment", "--layout", str(layout), "--read-list", str(read_list), "-o", str(out)],
    )
    assert res.exit_code == 0, res.output
    txt = out.read_text()
    assert "k0" in txt and "k1" in txt
    # k0 is U2OS-private; the summary should report a U2OS private/enriched cluster.
    assert "U2OS" in res.output
