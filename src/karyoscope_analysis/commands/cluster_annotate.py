"""``karyoscope-analysis cluster-annotate`` — label each cluster from its consensus structure.

Assigns a structural biological label (ECTR / subtelomere / Type II ALT subtelomere / interstitial
telomere / interstitial ITS-TAR1 / satellite-dominant) to each cluster, read off the consensus
using **hierarchy-derived feature classes** (no hardcoded feature names). Numeric thresholds
default to a documented human/CHM13 preset and are overridable. Replaces the legacy
`cluster_annotate.py` (per-read density aggregation + SVD feature-importance + sample-specific Type
I ALT relabel are not ported — the consensus subsumes the aggregation, and sample composition is
`test-enrichment`'s job).
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import cluster_annotate as core
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.clusters import read_clusters_table, read_consensus_segments


@click.command(
    name="cluster-annotate",
    help="Label each cluster from its consensus structure (hierarchy-derived).",
)
@click.option(
    "--clusters",
    "clusters_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` clusters.tsv (cluster_id, size, ..., width, ...).",
)
@click.option(
    "--consensus",
    "consensus_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` consensus BED (cluster_id, start, end, feature, ...).",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database hierarchy.tsv — defines the telomere/satellite/ITS-TAR1/chromosome classes.",
)
@click.option(
    "--min-cluster-size",
    default=2,
    show_default=True,
    type=int,
    help="Annotate only clusters with at least this many reads.",
)
@click.option(
    "--end-fraction",
    default=0.15,
    show_default=True,
    type=float,
    help="Fraction of the consensus span counted as an 'end' (for telomere-at-end detection).",
)
@click.option(
    "--satellite-fraction",
    default=0.8,
    show_default=True,
    type=float,
    help="Satellite bp / span needed to call a cluster satellite-dominant.",
)
@click.option(
    "--alt-block-bp",
    default=6000,
    show_default=True,
    type=int,
    help="Contiguous canonical-telomere bp for a subtelomere to be a Type II ALT subtelomere.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output cluster-annotation TSV.",
)
def cmd(
    clusters_path: Path,
    consensus_path: Path,
    hierarchy_path: Path,
    min_cluster_size: int,
    end_fraction: float,
    satellite_fraction: float,
    alt_block_bp: int,
    output: Path,
) -> None:
    """Label each cluster from its consensus and write a cluster-annotation TSV."""
    try:
        cluster_sizes, cluster_widths = read_clusters_table(clusters_path)
        consensus_by_cluster = read_consensus_segments(consensus_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)

    cfg = core.LabelConfig(
        end_fraction=end_fraction,
        satellite_fraction=satellite_fraction,
        alt_block_bp=alt_block_bp,
    )
    rows = core.annotate(
        cluster_sizes,
        cluster_widths,
        consensus_by_cluster,
        hierarchy,
        cfg=cfg,
        min_cluster_size=min_cluster_size,
    )
    output.write_text(core.annotation_tsv(rows))

    counts: dict[str, int] = {}
    for r in rows:
        counts[r.label or "(unlabeled)"] = counts.get(r.label or "(unlabeled)", 0) + 1
    breakdown = ", ".join(
        f"{label}={n}" for label, n in sorted(counts.items(), key=lambda x: -x[1])
    )
    click.echo(f"Annotated {len(rows)} cluster(s) -> {output}")
    if breakdown:
        click.echo(f"  labels: {breakdown}")
