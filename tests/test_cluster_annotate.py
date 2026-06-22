"""Tests for the consensus-signature cluster labeler and the cluster-annotate command."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import cluster_annotate as ca
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


def _h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def test_label_ectr_telomere_both_ends():
    # canonical telomere at both ends of the span -> ECTR.
    segs = [
        (0, 2000, "chr7:canonical_telomere"),
        (2000, 8000, "chr7:q_arm"),
        (8000, 10000, "chr7:canonical_telomere"),
    ]
    assert ca.label_cluster(segs, _h(), ca.LabelConfig()) == "ECTR"


def test_label_subtelomere_one_end_and_type_ii_alt():
    # Telomere at one end only -> subtelomere; a long canonical block -> Type II ALT subtelomere.
    short = [(0, 1000, "chr7:canonical_telomere"), (1000, 10000, "chr7:q_arm")]
    assert ca.label_cluster(short, _h(), ca.LabelConfig()) == "subtelomere"
    longblk = [(0, 7000, "chr7:canonical_telomere"), (7000, 10000, "chr7:q_arm")]
    assert ca.label_cluster(longblk, _h(), ca.LabelConfig()) == "Type II ALT subtelomere"


def test_label_interstitial_telomere():
    segs = [
        (0, 4000, "chr7:q_arm"),
        (4000, 5000, "chr7:canonical_telomere"),
        (5000, 10000, "chr7:q_arm"),
    ]
    assert ca.label_cluster(segs, _h(), ca.LabelConfig()) == "interstitial telomere"


def test_label_interstitial_its_tar1():
    segs = [(0, 4000, "chr7:q_arm"), (4000, 5000, "chr7:ITS"), (5000, 10000, "chr7:q_arm")]
    assert ca.label_cluster(segs, _h(), ca.LabelConfig()) == "interstitial ITS/TAR1"


def test_label_satellite_dominant():
    segs = [(0, 9000, "chr7:aSat"), (9000, 10000, "chr7:q_arm")]  # 90% satellite
    assert ca.label_cluster(segs, _h(), ca.LabelConfig()) == "satellite-dominant"


def test_label_unlabeled_arm_only():
    segs = [(0, 10000, "chr7:q_arm")]
    assert ca.label_cluster(segs, _h(), ca.LabelConfig()) == ""


def test_chromosomes_of_uses_specific_chromosomes_only():
    segs = [
        (0, 100, "chr2:q_arm"),
        (100, 200, "chr13:p_arm"),
        (200, 300, "autosome:arm"),  # grouping node -> not a specific chromosome
    ]
    assert ca.chromosomes_of(segs, _h()) == ["chr2", "chr13"]  # natural order


def test_cluster_annotate_cli(cli_runner, tmp_path: Path):
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text(
        "cluster_id\tsize\tseed\tconsensus_segments\twidth\torientation_conflict\n"
        "cluster_0\t4\ts0\t3\t10000\t0\n"
        "cluster_1\t1\ts1\t1\t100\t0\n"  # singleton -> excluded
    )
    consensus = tmp_path / "consensus.bed"
    consensus.write_text(
        "cluster_id\tstart\tend\tfeature\tsupport\tcoverage\n"
        "cluster_0\t0\t2000\tchr7:canonical_telomere\t4\t4\n"
        "cluster_0\t2000\t8000\tchr7:q_arm\t4\t4\n"
        "cluster_0\t8000\t10000\tchr7:canonical_telomere\t4\t4\n"
        "cluster_1\t0\t100\tchr2:aSat\t1\t1\n"
    )
    out = tmp_path / "annot.tsv"
    res = cli_runner.invoke(
        main,
        [
            "cluster-annotate",
            "--clusters",
            str(clusters),
            "--consensus",
            str(consensus),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "-o",
            str(out),
        ],
    )
    assert res.exit_code == 0, res.output
    lines = out.read_text().splitlines()
    assert lines[0] == "cluster_id\tsize\twidth\tlabel\tchromosomes\tconsensus_signature"
    assert len(lines) == 2  # cluster_0 only (singleton excluded)
    assert "\tECTR\tchr7\t" in lines[1]
