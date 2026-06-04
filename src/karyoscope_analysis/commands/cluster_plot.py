"""``karyoscope-analysis cluster-plot`` — render an Engine B cluster as SVG.

Reads `cluster`'s ``layout.tsv`` + ``consensus.bed`` plus the overlay BED that was clustered,
and draws one cluster as feature-colored read tracks aligned to the seed (see
``core/cluster_plot.py``). Self-contained SVG; the single read-renderer for the package.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import cluster_plot as render
from karyoscope_analysis.core.io.bed import read_annotation_bed
from karyoscope_analysis.core.io.colors import load_colors


def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _dominant_chromosome(consensus: list[tuple[int, int, str]]) -> str:
    """Most-covered chromosome layer of a consensus (``chr:feature`` labels), or ''."""
    by_chrom: dict[str, int] = {}
    for s, e, feature in consensus:
        chrom, sep, _ = feature.partition(":")
        if sep:
            by_chrom[chrom] = by_chrom.get(chrom, 0) + (e - s)
    return max(by_chrom, key=by_chrom.get) if by_chrom else ""  # type: ignore[arg-type]


@click.command(name="cluster-plot", help="Render Engine B cluster(s) as an SVG of read tracks.")
@click.option(
    "--layout",
    "layout_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` layout TSV (cluster_id, read_id, is_seed, reversed, offset, length).",
)
@click.option(
    "--consensus",
    "consensus_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` consensus BED (cluster_id, start, end, feature, ...).",
)
@click.option(
    "--overlay",
    "overlay_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The overlay-annotations BED that was clustered (per-read feature segments).",
)
@click.option(
    "--colors",
    "colors_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database colors.tsv (feature_set, feature, color). Falls back to an auto-palette.",
)
@click.option(
    "--cluster-id",
    default=None,
    help="Render just this cluster (e.g. cluster_0). Omit to render all clusters in one SVG.",
)
@click.option(
    "--min-cluster-size",
    default=2,
    show_default=True,
    type=int,
    help="When rendering all clusters, include only clusters with at least this many reads.",
)
@click.option(
    "--min-segment-bp",
    default=0,
    show_default=True,
    type=int,
    help="Drop feature segments shorter than this when drawing (visual denoise).",
)
@click.option("--width", default=1200, show_default=True, type=int)
@click.option("--row-height", default=12, show_default=True, type=int)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output SVG path.",
)
def cmd(
    layout_path: Path,
    consensus_path: Path,
    overlay_path: Path,
    colors_path: Path | None,
    cluster_id: str | None,
    min_cluster_size: int,
    min_segment_bp: int,
    width: int,
    row_height: int,
    output: Path,
) -> None:
    """Render one cluster, or all clusters stacked in one SVG."""
    overlay = read_annotation_bed(overlay_path)
    colors = load_colors(colors_path) if colors_path else {}

    def keep(segments):
        return [(s, e, f) for s, e, f in segments if e - s >= min_segment_bp]

    # Group layout + consensus rows by cluster, preserving file order (size-descending).
    layout_by_cluster: dict[str, list[dict[str, str]]] = {}
    for r in _read_tsv(layout_path):
        layout_by_cluster.setdefault(r["cluster_id"], []).append(r)
    consensus_by_cluster: dict[str, list[tuple[int, int, str]]] = {}
    for r in _read_tsv(consensus_path):
        consensus_by_cluster.setdefault(r["cluster_id"], []).append(
            (int(r["start"]), int(r["end"]), r["feature"])
        )

    if cluster_id is not None:
        if cluster_id not in layout_by_cluster:
            raise click.ClickException(f"cluster {cluster_id!r} not found in {layout_path}")
        selected = [cluster_id]
    else:
        selected = [c for c, rows in layout_by_cluster.items() if len(rows) >= min_cluster_size]
        if not selected:
            raise click.ClickException(f"no clusters with >= {min_cluster_size} reads")

    panels: list[render.ClusterPanel] = []
    for cid in selected:
        placed: list[render.PlacedRead] = []
        for r in layout_by_cluster[cid]:
            read_id = r["read_id"]
            if read_id not in overlay:
                raise click.ClickException(f"read {read_id!r} ({cid}) absent from {overlay_path}")
            placed.append(
                render.PlacedRead(
                    read_id=read_id,
                    is_seed=r["is_seed"] == "1",
                    reversed=r["reversed"] == "1",
                    offset=int(r["offset"]),
                    segments=keep(overlay[read_id]),
                )
            )
        consensus = keep(consensus_by_cluster.get(cid, []))
        chrom = _dominant_chromosome(consensus_by_cluster.get(cid, []))
        title = f"{cid}  n={len(placed)}" + (f"  {chrom}" if chrom else "")
        panels.append(render.ClusterPanel(title=title, placed=placed, consensus=consensus))

    svg = render.render_clusters_svg(panels, colors, width=width, row_height=row_height)
    output.write_text(svg)
    n_reads = sum(len(p.placed) for p in panels)
    click.echo(f"Rendered {len(panels)} cluster(s), {n_reads} reads, to {output}")
