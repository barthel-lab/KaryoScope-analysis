"""``karyoscope-analysis plot-reads`` — render per-read feature BEDs as stacked bars (SVG).

Reads one or more per-read feature BEDs (``read_id  start  end feature``), optionally reorients
them so a feature class sits at the top, and draws each read as a bar of feature-colored runs —
stacked as columns (vertical, default) or rows (``--horizontal``) — with a scale bar and an
optional auto-filtered legend. Rendering goes through ``karyoplot`` (``rasterize_features`` +
drawsvg). Feature colors come from the DB ``colors.tsv``; only ``novel`` may be absent (it renders
white) — any other uncolored feature is an error. ``--orient`` resolves its telomere / chromosome
/ satellite feature class from the DB ``hierarchy.tsv`` (no hardcoded feature lists).

A ``--read-list`` TSV adds two-tier grouping (``--label-tier``), a metadata heatmap above the
reads (``--heatmap``/``--heatmap-track``; categorical colors from karyoplot's shared palette),
and group filtering (``--filter-group``); ``--markers`` draws arrowheads at given bp positions.
``--format png``/``both`` converts the SVG to PNG via ``rsvg-convert``; ``--preset telogator``
reproduces the legacy ``telogator-reads-viz`` defaults (telomere orientation, SVG+PNG). The
panning animation (D7) remains deferred.
"""

from __future__ import annotations

from pathlib import Path

import click
from karyoplot.svg.export import RsvgConvertMissingError, svg_to_png

from karyoscope_analysis.core import plot_reads as render
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.colors import load_colors
from karyoscope_analysis.core.legend_order import feature_sort_key


@click.command(
    name="plot-reads", help="Render per-read feature BEDs as stacked colored bars (SVG)."
)
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
    help="Output path. The extension is set from --format (.svg / .png); 'both' writes each.",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["svg", "png", "both"]),
    default="svg",
    show_default=True,
    help="Output format. PNG is produced by converting the SVG (needs rsvg-convert).",
)
@click.option(
    "--png-scale", type=int, default=4, show_default=True, help="PNG zoom factor (rsvg -z)."
)
@click.option(
    "--preset",
    type=click.Choice(["none", "telogator"]),
    default="none",
    show_default=True,
    help="Convenience preset. 'telogator' defaults --orient telomere and --format both "
    "(reproduces the legacy telogator-reads-viz; explicit flags still win).",
)
@click.option(
    "--horizontal", is_flag=True, help="Draw reads as horizontal rows (default: vertical columns)."
)
@click.option(
    "--orient",
    type=click.Choice(["none", "telomere", "chromosome", "satellite"]),
    default="none",
    show_default=True,
    help="Reorient each read so this feature class is at the top (position 0). Needs the DB "
    "hierarchy (--hierarchy, or hierarchy.tsv beside --colors) to resolve the feature class.",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Database hierarchy.tsv (for --orient). Default: hierarchy.tsv next to --colors.",
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
@click.option(
    "--no-scale-bar", "scale_bar", flag_value=False, default=True, help="Omit the scale bar."
)
@click.option("--no-header", is_flag=True, help="Omit sample labels and separator lines.")
@click.option("--read-border", is_flag=True, help="Outline each read bar.")
@click.option(
    "--bar-width", type=int, default=5, show_default=True, help="Read bar thickness (px)."
)
@click.option(
    "--read-spacing", type=int, default=5, show_default=True, help="Spacing between reads (px)."
)
@click.option(
    "--sample-spacing",
    type=int,
    default=20,
    show_default=True,
    help="Spacing between samples/groups (px).",
)
@click.option(
    "--subgroup-spacing",
    type=int,
    default=None,
    help="Spacing within a group (px; default: sample-spacing).",
)
@click.option("--ratio", type=float, default=1 / 300, show_default=True, help="bp-to-pixel ratio.")
@click.option(
    "--font-size", type=int, default=11, show_default=True, help="Label/scale-bar text size."
)
@click.option("--min-length", type=int, default=None, help="Drop reads shorter than this (bp).")
@click.option("--max-length", type=int, default=None, help="Drop reads longer than this (bp).")
@click.option(
    "--scale-bar-bp", type=int, default=None, help="Scale bar size in bp (default auto: 10kb/5kb)."
)
@click.option(
    "--read-list",
    "read_list_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="TSV (read_id + columns) to filter/group reads and supply heatmap metadata.",
)
@click.option(
    "--label-tier",
    "label_tiers",
    multiple=True,
    help="Read-list column for a label tier as 'COLUMN[:DISPLAY]'; repeatable (group, then subgroup). "
    "Default: the first two columns.",
)
@click.option(
    "--heatmap", is_flag=True, help="Draw a metadata heatmap above the reads (needs --read-list)."
)
@click.option(
    "--heatmap-track",
    "heatmap_tracks",
    multiple=True,
    help="Read-list column for a heatmap row as 'COLUMN[:DISPLAY]'; repeatable (implies --heatmap). "
    "Default: read-list columns after the label tiers.",
)
@click.option(
    "--filter-group",
    "filter_groups",
    multiple=True,
    help="Keep only reads whose group is in these values.",
)
@click.option(
    "--markers",
    "markers_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="TSV (read_id, start, end) drawing left-pointing arrowheads at those bp positions.",
)
@click.option(
    "--marker-scale", type=float, default=1.0, show_default=True, help="Arrowhead size factor."
)
@click.pass_context
def cmd(
    ctx: click.Context,
    bed_specs: tuple[str, ...],
    colors_path: Path,
    output: Path,
    fmt: str,
    png_scale: int,
    preset: str,
    horizontal: bool,
    orient: str,
    hierarchy_path: Path | None,
    feature_mode: str,
    background: str,
    legend: bool,
    scale_bar: bool,
    no_header: bool,
    read_border: bool,
    bar_width: int,
    read_spacing: int,
    sample_spacing: int,
    subgroup_spacing: int | None,
    ratio: float,
    font_size: int,
    min_length: int | None,
    max_length: int | None,
    scale_bar_bp: int | None,
    read_list_path: Path | None,
    label_tiers: tuple[str, ...],
    heatmap: bool,
    heatmap_tracks: tuple[str, ...],
    filter_groups: tuple[str, ...],
    markers_path: Path | None,
    marker_scale: float,
) -> None:
    """Render per-read feature BEDs as an SVG (and/or PNG) of stacked colored bars."""
    if preset == "telogator":  # legacy telogator-reads-viz defaults; explicit flags still win
        default = click.core.ParameterSource.DEFAULT
        if ctx.get_parameter_source("orient") == default:
            orient = "telomere"
        if ctx.get_parameter_source("fmt") == default:
            fmt = "both"

    colors = load_colors(colors_path)
    try:
        reads, sample_order = render.load_bed_specs(list(bed_specs))
    except FileNotFoundError as e:
        raise click.ClickException(str(e)) from e

    if min_length is not None:
        reads = [r for r in reads if r.length >= min_length]
    if max_length is not None:
        reads = [r for r in reads if r.length <= max_length]

    # Read-list: filter + grouping + heatmap metadata.
    rl = None
    if read_list_path is not None:
        rl = render.process_read_list(
            reads,
            read_list_path,
            label_tier_specs=label_tiers,
            heatmap_track_specs=heatmap_tracks,
            heatmap=heatmap or bool(heatmap_tracks),
            filter_groups=filter_groups,
        )
        reads = rl.reads
    elif heatmap or heatmap_tracks or label_tiers or filter_groups:
        raise click.UsageError(
            "--heatmap/--heatmap-track/--label-tier/--filter-group require --read-list"
        )

    if not reads:
        raise click.ClickException("no reads to plot after filtering")

    want_heatmap = bool(rl and rl.metadata_columns and (heatmap or heatmap_tracks))
    if want_heatmap and horizontal:
        click.echo("warning: --heatmap is not supported with --horizontal; disabling heatmap")
        want_heatmap = False

    if orient != "none":
        hpath = hierarchy_path or colors_path.parent / "hierarchy.tsv"
        if not hpath.exists():
            raise click.UsageError(
                f"--orient needs the DB hierarchy; none found at {hpath}. Pass --hierarchy."
            )
        vocab = FeatureHierarchy.from_tsv(hpath)
        primary, fallback = {
            "telomere": (vocab.telomere_features, None),
            "chromosome": (vocab.chromosomes, None),
            "satellite": (vocab.telomere_features, vocab.satellite_features),
        }[orient]
        reads = render.orient_reads(reads, set(primary), set(fallback) if fallback else None)

    # Two-tier grouping relabels reads to composite "group — subgroup" samples.
    label_tier_map: dict[str, tuple[str, str | None]] = {}
    tier_display_names: dict[int, str] = {}
    group_subgroup_order: list[tuple[str, str | None]] = []
    group_boundaries = frozenset()
    if rl and rl.read_groups:
        reads, sample_order, group_boundaries = render.apply_grouping(
            reads, rl.read_groups, rl.group_subgroup_order
        )
        group_subgroup_order = rl.group_subgroup_order
        label_tier_map = dict(zip(sample_order, group_subgroup_order, strict=True))
        tier_display_names = {i: name for i, (_c, name) in enumerate(rl.tier_specs)}

    heatmap_colors = (
        render.assign_heatmap_colors(rl.metadata_columns, rl.read_metadata) if want_heatmap else {}
    )
    markers = render.parse_markers(markers_path) if markers_path is not None else {}
    reads = render.sort_reads(reads, sample_order)

    # KaryoScope-style legend ordering when a hierarchy is available (next to --colors or --hierarchy).
    legend_key = None
    if legend:
        lpath = hierarchy_path or colors_path.parent / "hierarchy.tsv"
        if lpath.exists():
            legend_key = feature_sort_key(lpath)

    cfg = render.PlotConfig(
        bar_width=bar_width,
        read_spacing=read_spacing,
        sample_spacing=sample_spacing,
        subgroup_spacing=subgroup_spacing,
        ratio=ratio,
        background=background,
        font_size=font_size,
        feature_mode=feature_mode,
        read_border=read_border,
        draw_scale_bar=scale_bar,
        no_header=no_header,
        scale_bar_bp=scale_bar_bp,
        legend=legend,
        label_tiers=label_tier_map,
        group_subgroup_order=group_subgroup_order,
        tier_display_names=tier_display_names,
        group_boundaries=group_boundaries,
        metadata_columns=rl.metadata_columns if want_heatmap else (),
        read_metadata=rl.read_metadata if want_heatmap else {},
        heatmap_colors=heatmap_colors,
        heatmap_display_names=dict(rl.heatmap_specs) if want_heatmap else {},
        markers=markers,
        marker_scale=marker_scale,
        legend_sort_key=legend_key,
    )
    try:
        svg = render.render(reads, colors, cfg, sample_order, horizontal=horizontal)
    except render.UnknownFeatureError as e:
        raise click.ClickException(str(e)) from e

    # PNG is produced by converting the SVG; always write the SVG (it's the source of truth).
    svg_path = output.with_suffix(".svg")
    svg_path.write_text(svg)
    written = [svg_path]
    if fmt in ("png", "both"):
        png_path = output.with_suffix(".png")
        try:
            svg_to_png(svg_path, png_path, scale=png_scale, raise_on_error=True)
        except RsvgConvertMissingError as e:
            raise click.ClickException(f"{e} (needed for --format {fmt})") from e
        written.append(png_path)
        if fmt == "png":
            svg_path.unlink()  # png-only: the SVG was just an intermediate
            written = [png_path]

    paths = ", ".join(str(p) for p in written)
    click.echo(f"Rendered {len(reads)} reads from {len(sample_order)} sample(s) to {paths}")
