"""``karyoscope-analysis plot-reads`` — render per-read feature BEDs as stacked bars (SVG).

Phase 3a: the SVG core. Reads one or more per-read feature BEDs (``read_id  start  end
feature``), optionally reorients them so a feature class sits at the top, and draws each read
as a bar of feature-colored runs — stacked as columns (vertical, default) or rows
(``--horizontal``) — with a scale bar and an optional auto-filtered legend. Rendering goes
through ``karyoplot`` (``rasterize_features`` + drawsvg). Feature colors come from the DB
``colors.tsv``; features absent from it render white (the ``novel`` sentinel).

Heatmap/metadata tracks, read-list grouping and markers (3b) and PNG/animation (3c) land later.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import plot_reads as render
from karyoscope_analysis.core.io.colors import load_colors


@click.command(name="plot-reads", help="Render per-read feature BEDs as stacked colored bars (SVG).")
@click.option(
    "--bed",
    "bed_specs",
    multiple=True,
    required=True,
    help="Per-read feature BED as 'SAMPLE:PATH' (or bare 'PATH'; sample = filename stem). "
    "Repeatable; sample order follows first appearance.",
)
@click.option(
    "--colors",
    "colors_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database colors.tsv (feature_set, feature, color). Features absent from it render white.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output SVG path.",
)
@click.option("--horizontal", is_flag=True, help="Draw reads as horizontal rows (default: vertical columns).")
@click.option(
    "--orient",
    type=click.Choice(["none", "telomere", "chromosome", "satellite"]),
    default="none",
    show_default=True,
    help="Reorient each read so this feature class is at the top (position 0).",
)
@click.option(
    "--feature-mode",
    type=click.Choice(["smooth", "transition", "raw"]),
    default="smooth",
    show_default=True,
    help="Rasterization: smooth (windowed majority vote), transition (min-width), raw (bp rects).",
)
@click.option(
    "--background",
    type=click.Choice(["black", "white"]),
    default="black",
    show_default=True,
    help="Background color (sets text/line color to its complement).",
)
@click.option("--legend", is_flag=True, help="Draw an auto-filtered color legend below the reads.")
@click.option("--no-scale-bar", "scale_bar", flag_value=False, default=True, help="Omit the scale bar.")
@click.option("--no-header", is_flag=True, help="Omit sample labels and separator lines.")
@click.option("--read-border", is_flag=True, help="Outline each read bar.")
@click.option("--bar-width", type=int, default=5, show_default=True, help="Read bar thickness (px).")
@click.option("--read-spacing", type=int, default=5, show_default=True, help="Spacing between reads (px).")
@click.option("--sample-spacing", type=int, default=20, show_default=True, help="Spacing between samples (px).")
@click.option("--ratio", type=float, default=1 / 300, show_default=True, help="bp-to-pixel ratio.")
@click.option("--font-size", type=int, default=11, show_default=True, help="Label/scale-bar text size.")
@click.option("--min-length", type=int, default=None, help="Drop reads shorter than this (bp).")
@click.option("--max-length", type=int, default=None, help="Drop reads longer than this (bp).")
@click.option("--scale-bar-bp", type=int, default=None, help="Scale bar size in bp (default auto: 10kb/5kb).")
def cmd(
    bed_specs: tuple[str, ...],
    colors_path: Path,
    output: Path,
    horizontal: bool,
    orient: str,
    feature_mode: str,
    background: str,
    legend: bool,
    scale_bar: bool,
    no_header: bool,
    read_border: bool,
    bar_width: int,
    read_spacing: int,
    sample_spacing: int,
    ratio: float,
    font_size: int,
    min_length: int | None,
    max_length: int | None,
    scale_bar_bp: int | None,
) -> None:
    """Render per-read feature BEDs as an SVG of stacked colored bars."""
    colors = load_colors(colors_path)
    try:
        reads, sample_order = render.load_bed_specs(list(bed_specs))
    except FileNotFoundError as e:
        raise click.ClickException(str(e)) from e

    if min_length is not None:
        reads = [r for r in reads if r.length >= min_length]
    if max_length is not None:
        reads = [r for r in reads if r.length <= max_length]
    if not reads:
        raise click.ClickException("no reads to plot after length filtering")

    if orient != "none":
        reads = render.orient_reads(reads, orient)
    reads = render.sort_reads(reads, sample_order)

    cfg = render.PlotConfig(
        bar_width=bar_width, read_spacing=read_spacing, sample_spacing=sample_spacing,
        ratio=ratio, background=background, font_size=font_size, feature_mode=feature_mode,
        read_border=read_border, draw_scale_bar=scale_bar, no_header=no_header,
        scale_bar_bp=scale_bar_bp, legend=legend,
    )
    svg = render.render(reads, colors, cfg, sample_order, horizontal=horizontal)
    output.write_text(svg)
    click.echo(f"Rendered {len(reads)} reads from {len(sample_order)} sample(s) to {output}")
