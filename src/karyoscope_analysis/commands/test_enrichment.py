"""``karyoscope-analysis test-enrichment`` — per-cluster enrichment across samples/groups.

Consumes ``cluster``'s ``layout.tsv`` (read -> cluster) and a read-list TSV mapping each read to a
sample/group, and writes a per-cluster enrichment table: per-group read counts, depth-normalized
fractions, log2 fold-enrichment vs the pool, and ``private``/``enriched`` flags. Descriptive and
effect-size-first — see :mod:`karyoscope_analysis.core.enrichment` for the model and its
no-replication caveat.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import click

from karyoscope_analysis.core import enrichment as enrich
from karyoscope_analysis.core.io.clusters import read_layout_assignments


def _read_group_map(path: Path, group_col: str) -> tuple[list[str], dict[str, str]]:
    """Stream a read-list TSV -> ``(columns, {read_id: group_value})`` for one column only.

    test-enrichment needs only each read's group assignment, so this avoids
    ``load_read_list``'s full per-read metadata dict + row list -- bounding memory for a large
    pooled read-list. Header detection matches ``load_read_list``; last value wins per read.
    """
    p = str(path)
    open_func = gzip.open if p.endswith(".gz") else open
    mode = "rt" if p.endswith(".gz") else "r"
    columns: list[str] = []
    gi: int | None = None
    read_to_group: dict[str, str] = {}
    with open_func(p, mode) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            rid = fields[0]
            if rid.lower() in ("sequence", "read", "read_id"):
                if len(fields) > 1:
                    columns = fields[1:]
                    gi = columns.index(group_col) + 1 if group_col in columns else None
                continue
            if gi is not None and gi < len(fields) and fields[gi]:
                read_to_group[rid] = fields[gi]
    return columns, read_to_group


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
    try:
        read_to_cluster = read_layout_assignments(layout_path)
    except ValueError as e:
        raise click.ClickException(str(e)) from e

    columns, read_to_group = _read_group_map(read_list_path, group_col)
    if group_col not in columns:
        raise click.ClickException(f"--group-col {group_col!r} not in read-list columns {columns}")
    if not read_to_group:
        raise click.ClickException(f"no reads have a {group_col!r} value in {read_list_path}")

    matched = sum(1 for r in read_to_cluster if r in read_to_group)
    if matched == 0:
        raise click.ClickException(
            "no read IDs are shared between the layout and the read-list "
            "(check that the read-list IDs match the (pooled) read IDs used for clustering)"
        )

    results, groups, group_totals = enrich.compute_enrichment(
        read_to_cluster,
        read_to_group,
        min_log2_effect=min_log2_effect,
        min_cluster_size=min_cluster_size,
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
