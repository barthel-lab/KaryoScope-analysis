"""``karyoscope-analysis compare-clusterings`` — concordance of two clusterings of the same reads.

Takes two ``cluster`` ``layout.tsv`` files (e.g. different parameters or featuresets), computes the
Adjusted Rand Index and Normalized Mutual Information over their shared reads, and writes a report
plus the cluster-to-cluster read-overlap table.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import clustering_comparison as cc
from karyoscope_analysis.core.io.clusters import read_layout_assignments


@click.command(
    name="compare-clusterings",
    help="Compare two clusterings of the same reads (ARI/NMI + cluster overlap).",
)
@click.option(
    "--layout1",
    "layout1_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="First clustering's layout.tsv (cluster_id, read_id, ...).",
)
@click.option(
    "--layout2",
    "layout2_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Second clustering's layout.tsv.",
)
@click.option("--label1", default="A", show_default=True, help="Display label for clustering 1.")
@click.option("--label2", default="B", show_default=True, help="Display label for clustering 2.")
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Report output path; the cluster-overlap table is written alongside as *.overlap.tsv.",
)
def cmd(
    layout1_path: Path,
    layout2_path: Path,
    label1: str,
    label2: str,
    output: Path,
) -> None:
    """Compare two clusterings and write a concordance report + overlap table."""
    try:
        labels1 = read_layout_assignments(layout1_path)
        labels2 = read_layout_assignments(layout2_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    try:
        result = cc.compare(labels1, labels2)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    report = cc.report(result, label1, label2)
    output.write_text(report)
    overlap_path = output.with_suffix(".overlap.tsv")
    overlap_path.write_text(cc.overlap_tsv(cc.overlap_pairs(labels1, labels2), label1, label2))
    click.echo(report.rstrip("\n"))
    click.echo(f"Wrote {output} and {overlap_path}")
