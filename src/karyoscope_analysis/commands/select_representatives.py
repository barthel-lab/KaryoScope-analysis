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


def _read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise click.ClickException(f"empty file: {path}")
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:] if line]


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
    header, rows = _read_tsv(clusters_path)
    try:
        ci, si, wi = header.index("cluster_id"), header.index("size"), header.index("width")
    except ValueError as e:
        raise click.ClickException(
            f"{clusters_path}: expected 'cluster_id', 'size', 'width' columns, got {header}"
        ) from e
    cluster_sizes = {r[ci]: int(r[si]) for r in rows if len(r) > max(ci, si, wi)}
    cluster_widths = {r[ci]: int(r[wi]) for r in rows if len(r) > max(ci, si, wi)}

    cons_header, cons_rows = _read_tsv(consensus_path)
    try:
        cci, sti, eni, fti = (
            cons_header.index("cluster_id"), cons_header.index("start"),
            cons_header.index("end"), cons_header.index("feature"),
        )
    except ValueError as e:
        raise click.ClickException(
            f"{consensus_path}: expected 'cluster_id', 'start', 'end', 'feature' columns"
        ) from e
    consensus_by_cluster: dict[str, list[rep.Segment]] = {}
    for r in cons_rows:
        if len(r) > max(cci, sti, eni, fti):
            consensus_by_cluster.setdefault(r[cci], []).append((int(r[sti]), int(r[eni]), r[fti]))

    reps = rep.build_catalog(
        cluster_sizes, cluster_widths, consensus_by_cluster, min_cluster_size=min_cluster_size
    )
    output.write_text(rep.catalog_tsv(reps))
    click.echo(
        f"Cataloged {len(reps)} representative structure(s) "
        f"(clusters with >= {min_cluster_size} reads) -> {output}"
    )
