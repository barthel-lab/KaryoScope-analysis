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


@click.command(name="cluster-plot", help="Render one Engine B cluster as an SVG of read tracks.")
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
@click.option("--cluster-id", required=True, help="Which cluster to render (e.g. cluster_0).")
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
    cluster_id: str,
    min_segment_bp: int,
    width: int,
    row_height: int,
    output: Path,
) -> None:
    """Render one cluster to an SVG."""
    overlay = read_annotation_bed(overlay_path)
    colors = load_colors(colors_path) if colors_path else {}

    def keep(segments):
        return [(s, e, f) for s, e, f in segments if e - s >= min_segment_bp]

    layout_rows = [r for r in _read_tsv(layout_path) if r["cluster_id"] == cluster_id]
    if not layout_rows:
        raise click.ClickException(f"cluster {cluster_id!r} not found in {layout_path}")

    placed: list[render.PlacedRead] = []
    for r in layout_rows:
        read_id = r["read_id"]
        if read_id not in overlay:
            raise click.ClickException(
                f"read {read_id!r} (cluster {cluster_id}) absent from {overlay_path}"
            )
        placed.append(
            render.PlacedRead(
                read_id=read_id,
                is_seed=r["is_seed"] == "1",
                reversed=r["reversed"] == "1",
                offset=int(r["offset"]),
                segments=keep(overlay[read_id]),
            )
        )

    consensus = keep(
        (int(r["start"]), int(r["end"]), r["feature"])
        for r in _read_tsv(consensus_path)
        if r["cluster_id"] == cluster_id
    )

    svg = render.render_cluster_svg(
        placed, consensus, colors, width=width, row_height=row_height, title=cluster_id
    )
    output.write_text(svg)
    click.echo(f"Rendered {cluster_id} ({len(placed)} reads) to {output}")
