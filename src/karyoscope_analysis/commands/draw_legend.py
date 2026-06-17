"""``karyoscope-analysis draw-legend`` — standalone SVG legend from the DB palette.

Reads the database ``colors.tsv`` (``feature_set, feature, color``) and renders a standalone
legend grouped by feature set, via :mod:`karyoplot.svg.legend`. Every color comes from the DB —
nothing is hardcoded here. This replaces the legacy ``KaryoScope_draw_legend.py``: feature-set
grouping comes from the ``feature_set`` column rather than a manual ``Header:feat1,feat2`` string.
"""

from __future__ import annotations

from pathlib import Path

import click
from karyoplot.core import theme as themes
from karyoplot.svg.legend import featureset_legend_items, make_legend_drawing, merge_by_color

from karyoscope_analysis.core.io.colors import load_colors, load_colors_by_featureset
from karyoscope_analysis.core.legend_order import feature_sort_key


def _csv_set(value: str | None) -> set[str] | None:
    """Parse a comma-separated option into a set, or ``None`` if empty/unset."""
    if not value:
        return None
    return {s.strip() for s in value.split(",") if s.strip()}


def _group_layout(items: list[tuple[str, str, bool]]) -> tuple[int, int]:
    """Return ``(rows, cols)`` placing each feature-set group in its own column.

    ``cols`` is the number of header rows (one per group); ``rows`` is the size of
    the largest group. With ``rows >= every group`` and ``cols == n_groups``,
    ``make_legend_drawing`` starts each group cleanly at the top of a new column
    and never wraps a group mid-column or drops groups past ``cols``.
    """
    sizes: list[int] = []
    for _label, _color, is_header in items:
        if is_header:
            sizes.append(0)
        elif sizes:
            sizes[-1] += 1
    n_groups = len(sizes)
    return (max(sizes, default=1) or 1, max(n_groups, 1))


@click.command(name="draw-legend", help="Render a standalone SVG legend from the DB color palette.")
@click.option(
    "--colors",
    "colors_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database colors.tsv (feature_set, feature, color).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output SVG path.",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Database hierarchy.tsv for KaryoScope-style feature ordering within each set. "
    "Default: hierarchy.tsv next to --colors if present (else colors.tsv order).",
)
@click.option(
    "--feature-set",
    "feature_sets",
    multiple=True,
    help="Limit to these feature set(s); repeatable. Default: all sets, in colors.tsv order.",
)
@click.option("--include", default=None, help="Comma-separated feature names to include.")
@click.option("--exclude", default=None, help="Comma-separated feature names to exclude.")
@click.option(
    "--merge-same-color",
    is_flag=True,
    help="Collapse features sharing a hex color into one row (drops feature-set grouping).",
)
@click.option(
    "--clean-labels/--raw-labels",
    default=True,
    show_default=True,
    help="Strip suffixes and underscores from labels, or show feature names verbatim.",
)
@click.option(
    "--theme",
    "theme_name",
    type=click.Choice(["dark", "light"]),
    default="dark",
    show_default=True,
    help="Color theme (background/text).",
)
@click.option(
    "--columns",
    type=int,
    default=None,
    help="Number of columns (auto if unset). Used only with --merge-same-color; the grouped "
    "layout always uses one column per feature set so groups can't be truncated.",
)
@click.option(
    "--rows",
    type=int,
    default=None,
    help="Number of rows (auto if unset). Used only with --merge-same-color.",
)
@click.option("--swatch-size", type=int, default=12, show_default=True)
@click.option("--font-size", type=int, default=12, show_default=True)
def cmd(
    colors_path: Path,
    output: Path,
    hierarchy_path: Path | None,
    feature_sets: tuple[str, ...],
    include: str | None,
    exclude: str | None,
    merge_same_color: bool,
    clean_labels: bool,
    theme_name: str,
    columns: int | None,
    rows: int | None,
    swatch_size: int,
    font_size: int,
) -> None:
    """Render a standalone legend SVG from the DB color palette."""
    theme = themes.get(theme_name)

    # KaryoScope-style within-featureset ordering if a hierarchy is available.
    hpath = hierarchy_path or colors_path.parent / "hierarchy.tsv"
    sort_key = feature_sort_key(hpath) if hpath.exists() else None

    if merge_same_color:
        collapsed = load_colors(colors_path)
        merged = merge_by_color(list(collapsed.items()))
        items = [(label, color, False) for label, color in merged]
    else:
        items = featureset_legend_items(
            load_colors_by_featureset(colors_path),
            feature_sets=list(feature_sets) or None,
            include=_csv_set(include),
            exclude=_csv_set(exclude),
            clean_labels=clean_labels,
            sort_key=sort_key,
        )
        # Grouped layout: one column per feature set, with `rows` tall enough for the
        # largest group, so make_legend_drawing never wraps a group mid-column or
        # truncates groups past `cols` (which a small --columns would otherwise do).
        rows, columns = _group_layout(items)

    if not items:
        raise click.ClickException(
            "no legend entries after filtering; check --feature-set/--include/--exclude"
        )

    drawing = make_legend_drawing(
        items,
        theme=theme,
        rows=rows,
        cols=columns,
        swatch_size=swatch_size,
        font_size=font_size,
    )
    drawing.save_svg(str(output))
    click.echo(f"Saved legend: {output} ({drawing.width:.0f} x {drawing.height:.0f} px)")
