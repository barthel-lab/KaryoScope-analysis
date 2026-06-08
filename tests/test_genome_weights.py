"""Tests for genome-frequency feature weights (core + genome-weights CLI)."""

from __future__ import annotations

import math
from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import genome_weights as gw

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


def test_tally_feature_bp():
    streams = {"region": iter([("chr1", 0, 100, "arm"), ("chr1", 100, 110, "aSat"),
                               ("chr2", 0, 50, "arm")])}
    assert gw.tally_feature_bp(streams) == {"region": {"arm": 150, "aSat": 10}}


def test_compute_weights_rare_feature_is_one():
    # arm tiles 99% of the partition (uninformative), aSat 1% (distinctive).
    weights = gw.compute_genome_weights({"region": {"arm": 990, "aSat": 10}})
    by_feature = {w.feature: w for w in weights}
    assert by_feature["aSat"].weight == 1.0  # rarest -> max info content -> 1
    assert by_feature["arm"].weight < by_feature["aSat"].weight  # ubiquitous -> down-weighted
    # weight = ic / ic_max = -ln(0.99) / -ln(0.01)
    assert math.isclose(by_feature["arm"].weight, -math.log(0.99) / -math.log(0.01))


def test_compute_weights_common_scale_across_featuresets():
    # the global max info content (rarest feature anywhere) anchors weight 1.
    weights = gw.compute_genome_weights(
        {"region": {"arm": 900, "aSat": 100}, "subtelomeric": {"canonical_telomere": 5, "nonsubtelomeric": 995}}
    )
    top = max(weights, key=lambda w: w.weight)
    assert top.feature == "canonical_telomere" and top.weight == 1.0  # rarest overall
    by = {w.feature: w.weight for w in weights}
    assert by["arm"] < by["aSat"] < by["canonical_telomere"]


def test_structural_weight_map_keeps_largest():
    weights = [
        gw.FeatureWeight("region", "rDNA", 100, 0.1, 2.3, 0.4),
        gw.FeatureWeight("acrocentric", "rDNA", 50, 0.2, 1.6, 0.7),
    ]
    assert gw.structural_weight_map(weights)["rDNA"] == 0.7  # most-distinctive interpretation


def test_genome_weights_cli_and_loader(cli_runner, tmp_path: Path):
    region = tmp_path / "ref.region.bed"
    region.write_text("chr1\t0\t990\tarm\nchr1\t990\t1000\taSat\n")
    subtel = tmp_path / "ref.subtelomeric.bed"
    subtel.write_text("chr1\t0\t5\tcanonical_telomere\nchr1\t5\t1000\tnonsubtelomeric\n")
    out = tmp_path / "weights.tsv"
    res = cli_runner.invoke(
        main,
        [
            "genome-weights",
            "--bed", f"region={region}",
            "--bed", f"subtelomeric={subtel}",
            "--hierarchy", str(HIERARCHY_TSV),
            "-o", str(out),
        ],
    )
    assert res.exit_code == 0, res.output
    loaded = gw.load_structural_weights(out)
    # canonical_telomere is rarest overall -> weight 1; arm is most ubiquitous -> smallest.
    assert loaded["canonical_telomere"] == 1.0
    assert loaded["arm"] < loaded["aSat"] < loaded["canonical_telomere"]


def test_cluster_genome_freq_uses_weights(cli_runner, tmp_path: Path):
    # Two reads share only a heavily down-weighted q_arm block plus distinct telomere ends.
    overlay = tmp_path / "overlay.bed"
    overlay.write_text(
        "x\t0\t1000\tchr1:canonical_telomere\nx\t1000\t6000\tchr1:q_arm\n"
        "y\t0\t5000\tchr2:q_arm\ny\t5000\t6000\tchr2:canonical_telomere\n"
    )
    weights = tmp_path / "w.tsv"
    weights.write_text(
        "feature_set\tfeature\tgenome_bp\tgenome_fraction\tinfo_content\tweight\n"
        "region\tq_arm\t100\t5.000e-01\t0.6931\t0.0300\n"
        "subtelomeric\tcanonical_telomere\t1\t1.000e-03\t6.9078\t0.5000\n"
    )
    out = tmp_path / "clusters.tsv"
    res = cli_runner.invoke(
        main,
        [
            "cluster", "--input", str(overlay), "--hierarchy", str(HIERARCHY_TSV),
            "--weight-method", "genome-freq", "--genome-weights", str(weights),
            "--min-overlap-bp", "1000", "-o", str(out),
        ],
    )
    assert res.exit_code == 0, res.output
    # the shared q_arm (5000 bp x weight 0.03 = 150 weighted bp) is below --min-overlap-bp 1000,
    # so the two reads do NOT merge on arm alone.
    assert "2 clusters" in res.output


def test_cluster_genome_freq_requires_weights_path(cli_runner, tmp_path: Path):
    overlay = tmp_path / "overlay.bed"
    overlay.write_text("x\t0\t1000\tchr1:q_arm\n")
    res = cli_runner.invoke(
        main,
        ["cluster", "--input", str(overlay), "--hierarchy", str(HIERARCHY_TSV),
         "--weight-method", "genome-freq", "-o", str(tmp_path / "c.tsv")],
    )
    assert res.exit_code != 0
    assert "genome-freq requires --genome-weights" in res.output


def test_genome_weights_cli_rejects_unknown_feature(cli_runner, tmp_path: Path):
    region = tmp_path / "ref.region.bed"
    region.write_text("chr1\t0\t100\tarm\nchr1\t100\t200\tnot_a_feature\n")
    out = tmp_path / "weights.tsv"
    res = cli_runner.invoke(
        main,
        ["genome-weights", "--bed", f"region={region}", "--hierarchy", str(HIERARCHY_TSV), "-o", str(out)],
    )
    assert res.exit_code != 0
    assert "not_a_feature" in res.output
