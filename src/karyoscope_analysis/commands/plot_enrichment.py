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
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.clusters import read_consensus_segments
from karyoscope_analysis.core.io.colors import load_colors
from karyoscope_analysis.core.legend_order import feature_sort_key


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
    "--consensus",
    "consensus_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="`cluster` consensus BED. With --colors, draws each row's consensus structure beside the "
    "heatmap (to check the label against the actual structure).",
)
@click.option(
    "--colors",
    "colors_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Database colors.tsv for the consensus-structure panel (needs --consensus).",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Database hierarchy.tsv to order the consensus feature legend (default: beside --colors).",
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
@click.option(
    "--clamp", type=float, default=4.0, show_default=True, help="Color-scale limit (log2)."
)
@click.option(
    "--normalize-consensus",
    is_flag=True,
    help="Normalize each consensus bar to its own width (compare structure only). Default: a "
    "shared absolute bp scale, so cluster lengths are comparable.",
)
@click.option(
    "--align-telomere/--no-align-telomere",
    default=True,
    show_default=True,
    help="Orient each consensus telomere-left and align rows at the telomere->rest breakpoint "
    "(needs the DB hierarchy for the telomere class). Ignored with --normalize-consensus.",
)
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
    consensus_path: Path | None,
    colors_path: Path | None,
    hierarchy_path: Path | None,
    max_clusters: int | None,
    all_clusters: bool,
    clamp: float,
    normalize_consensus: bool,
    align_telomere: bool,
    dark: bool,
    output: Path,
) -> None:
    """Render the clusters x samples enrichment heatmap."""
    enrichment = _read_tsv(enrichment_path)
    if not enrichment:
        raise click.ClickException(f"no rows in {enrichment_path}")
    # Groups are the log2fc_<group> columns, in file order.
    groups = [c[len("log2fc_") :] for c in enrichment[0] if c.startswith("log2fc_")]
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

    # Optional consensus-structure panel.
    consensus = colors = sort_key = telomere = None
    if consensus_path is not None or colors_path is not None:
        if consensus_path is None or colors_path is None:
            raise click.UsageError("--consensus and --colors must be given together")
        consensus = read_consensus_segments(consensus_path)
        colors = load_colors(colors_path)
        hpath = hierarchy_path or colors_path.parent / "hierarchy.tsv"
        if hpath.exists():
            sort_key = feature_sort_key(hpath)
            telomere = set(FeatureHierarchy.from_tsv(hpath).telomere_features)

    ep.render_heatmap(
        rows,
        groups,
        str(output),
        clamp=clamp,
        dark_mode=dark,
        consensus=consensus,
        colors=colors,
        sort_key=sort_key,
        normalize_consensus=normalize_consensus,
        telomere=telomere,
        align_telomere=align_telomere,
    )
    extra = " + consensus" if consensus else ""
    click.echo(
        f"Rendered enrichment heatmap{extra}: {len(rows)} clusters x {len(groups)} groups -> {output}"
    )
