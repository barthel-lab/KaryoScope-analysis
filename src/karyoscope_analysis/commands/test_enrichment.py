"""``karyoscope-analysis test-enrichment`` — per-cluster enrichment across samples/groups.

Consumes ``cluster``'s ``layout.tsv`` (read -> cluster) and a read-list TSV mapping each read to a
sample/group, and writes a per-cluster enrichment table: per-group read counts, depth-normalized
fractions, log2 fold-enrichment vs the pool, and ``private``/``enriched`` flags. Descriptive and
effect-size-first — see :mod:`karyoscope_analysis.core.enrichment` for the model and its
no-replication caveat.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import enrichment as enrich
from karyoscope_analysis.core.plot_reads import load_read_list


def _read_to_cluster(layout_path: Path) -> dict[str, str]:
    """Map ``read_id -> cluster_id`` from a ``layout.tsv`` (all of a read's segments share it)."""
    lines = layout_path.read_text().splitlines()
    if not lines:
        raise click.ClickException(f"empty layout file: {layout_path}")
    header = lines[0].split("\t")
    try:
        ci, ri = header.index("cluster_id"), header.index("read_id")
    except ValueError as e:
        raise click.ClickException(
            f"{layout_path}: expected 'cluster_id' and 'read_id' columns, got {header}"
        ) from e
    out: dict[str, str] = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) > max(ci, ri):
            out[f[ri]] = f[ci]
    return out


@click.command(name="test-enrichment", help="Per-cluster enrichment across samples/groups.")
@click.option(
    "--layout",
    "layout_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` layout.tsv (cluster_id, read_id, ...).",
)
@click.option(
    "--read-list",
    "read_list_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="TSV mapping read_id to sample/group (read_id + a group column; e.g. pool-samples output).",
)
@click.option(
    "--group-col",
    default="sample",
    show_default=True,
    help="Read-list column that partitions reads into groups to test enrichment across.",
)
@click.option(
    "--min-log2-effect",
    default=1.0,
    show_default=True,
    type=float,
    help="Flag a cluster 'enriched' when a group's log2 fold-enrichment reaches this (1.0 = 2x).",
)
@click.option(
    "--min-cluster-size",
    default=2,
    show_default=True,
    type=int,
    help="Minimum reads for a cluster to be callable 'enriched' (singletons carry no "
    "compositional evidence; they're still listed but never flagged enriched).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output enrichment TSV.",
)
def cmd(
    layout_path: Path,
    read_list_path: Path,
    group_col: str,
    min_log2_effect: float,
    min_cluster_size: int,
    output: Path,
) -> None:
    """Write a per-cluster enrichment table across samples/groups."""
    read_to_cluster = _read_to_cluster(layout_path)

    _ids, columns, read_data, _rows = load_read_list(read_list_path)
    if group_col not in columns:
        raise click.ClickException(
            f"--group-col {group_col!r} not in read-list columns {columns}"
        )
    read_to_group = {
        rid: meta[group_col] for rid, meta in read_data.items() if meta.get(group_col)
    }
    if not read_to_group:
        raise click.ClickException(f"no reads have a {group_col!r} value in {read_list_path}")

    matched = sum(1 for r in read_to_cluster if r in read_to_group)
    if matched == 0:
        raise click.ClickException(
            "no read IDs are shared between the layout and the read-list "
            "(check that the read-list IDs match the (pooled) read IDs used for clustering)"
        )

    results, groups, group_totals = enrich.compute_enrichment(
        read_to_cluster, read_to_group,
        min_log2_effect=min_log2_effect, min_cluster_size=min_cluster_size,
    )
    output.write_text(enrich.enrichment_tsv(results, groups))

    n_enriched = sum(1 for r in results if r.enriched)
    n_singleton = sum(1 for r in results if r.n_total < min_cluster_size)
    click.echo(
        f"Tested {len(results)} clusters ({n_singleton} below size {min_cluster_size}, not "
        f"callable) over {matched} reads in groups "
        f"{', '.join(f'{g}={group_totals[g]}' for g in groups)}; {n_enriched} enriched. "
        f"Wrote {output}"
    )
    click.echo(enrich.summarize(results, groups))
