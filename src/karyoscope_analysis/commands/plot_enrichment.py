"""``karyoscope-analysis plot-enrichment`` — heatmap of per-cluster enrichment across samples.

Reads ``test-enrichment``'s table (and, optionally, ``cluster-annotate``'s labels) and renders a
clusters x samples heatmap colored by log2 fold-enrichment — the capstone summary showing which
structural haplotypes concentrate in which line.
"""

from __future__ import annotations

import csv
from pathlib import Path

import click

from karyoscope_analysis.core import enrichment_plot as ep


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@click.command(name="plot-enrichment", help="Heatmap of per-cluster enrichment across samples.")
@click.option(
    "--enrichment",
    "enrichment_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`test-enrichment` output TSV.",
)
@click.option(
    "--annot",
    "annot_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="`cluster-annotate` output TSV, to label rows by their structural label.",
)
@click.option(
    "--max-clusters",
    type=int,
    default=None,
    help="Show at most this many clusters (the most-enriched). Default: all enriched.",
)
@click.option(
    "--all-clusters",
    is_flag=True,
    help="Include all (multi-read) clusters, not just the enriched ones.",
)
@click.option("--clamp", type=float, default=4.0, show_default=True, help="Color-scale limit (log2).")
@click.option("--dark", is_flag=True, help="Dark background.")
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output figure (.png/.pdf/.svg by extension).",
)
def cmd(
    enrichment_path: Path,
    annot_path: Path | None,
    max_clusters: int | None,
    all_clusters: bool,
    clamp: float,
    dark: bool,
    output: Path,
) -> None:
    """Render the clusters x samples enrichment heatmap."""
    enrichment = _read_tsv(enrichment_path)
    if not enrichment:
        raise click.ClickException(f"no rows in {enrichment_path}")
    # Groups are the log2fc_<group> columns, in file order.
    groups = [c[len("log2fc_"):] for c in enrichment[0] if c.startswith("log2fc_")]
    if not groups:
        raise click.ClickException(f"{enrichment_path}: no 'log2fc_<group>' columns found")

    labels: dict[str, str] = {}
    if annot_path is not None:
        labels = {r["cluster_id"]: r.get("label", "") for r in _read_tsv(annot_path)}

    rows = ep.select_rows(
        enrichment, groups, labels, enriched_only=not all_clusters, max_clusters=max_clusters
    )
    if not rows:
        raise click.ClickException("no clusters to plot (try --all-clusters)")

    ep.render_heatmap(rows, groups, str(output), clamp=clamp, dark_mode=dark)
    click.echo(f"Rendered enrichment heatmap: {len(rows)} clusters x {len(groups)} groups -> {output}")
