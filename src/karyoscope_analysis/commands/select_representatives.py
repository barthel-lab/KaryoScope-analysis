"""``karyoscope-analysis select-representatives`` — catalog each cluster's representative structure.

The cluster consensus *is* the representative (Engine B computes it). This reads ``cluster``'s
``clusters.tsv`` + ``consensus.bed`` and writes a per-cluster catalog: size, consensus segment
count, width, and a compact consensus *signature* (the ordered feature path) — a readable index of
the distinct structural haplotypes found. Replaces the legacy best-read selection (obsolete now
that the consensus exists and ``cluster-plot`` renders it).
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import representatives as rep
from karyoscope_analysis.core.io.clusters import read_clusters_table, read_consensus_segments


@click.command(
    name="select-representatives",
    help="Catalog each cluster's representative structure (its consensus).",
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
    "--min-cluster-size",
    default=2,
    show_default=True,
    type=int,
    help="Catalog only clusters with at least this many reads.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output representatives catalog TSV.",
)
def cmd(clusters_path: Path, consensus_path: Path, min_cluster_size: int, output: Path) -> None:
    """Write a per-cluster catalog of representative (consensus) structures."""
    try:
        cluster_sizes, cluster_widths = read_clusters_table(clusters_path)
        consensus_by_cluster = read_consensus_segments(consensus_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    reps = rep.build_catalog(
        cluster_sizes, cluster_widths, consensus_by_cluster, min_cluster_size=min_cluster_size
    )
    output.write_text(rep.catalog_tsv(reps))
    click.echo(
        f"Cataloged {len(reps)} representative structure(s) "
        f"(clusters with >= {min_cluster_size} reads) -> {output}"
    )
