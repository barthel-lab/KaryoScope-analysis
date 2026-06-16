"""Render Engine B clusters as SVG read tracks (consensus-coordinate layout).

The single read-renderer for the package. Each read is a horizontal row of feature-colored
rectangles placed in the cluster's **consensus coordinate frame** (computed by
`feature_assembly.consensus_layout` and serialized to `cluster`'s ``layout.tsv``), so matched
features **stack vertically** across reads; the union-spanning consensus track sits on top. A
gap in a row means that read has no feature there (it didn't align / doesn't extend that far) —
length filtering is the caller's choice, off by default, so gaps are unambiguous.

Rendering goes through ``karyoplot`` (the shared plotting library): segment rows are drawn with
:func:`karyoplot.svg.drawing.draw_annotation_track`, the feature/chromosome legend with
:func:`karyoplot.svg.legend.draw_grouped_legend`, onto a ``drawsvg`` Drawing — the same stack the
KaryoScope engine renders with. Engine-B specifics (composite ``chr:feature`` labels, seed
marking, the consensus row) stay here. Colors come from the DB palette: only ``novel`` may be
absent (it renders white); any other uncolored feature raises :class:`UnknownFeatureError`
rather than getting a fallback color. One cluster (:func:`render_cluster_svg`) or many stacked in
one figure (:func:`render_clusters_svg`), sharing a single feature legend.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import drawsvg as draw
from karyoplot.core.fonts import DEFAULT_FONT_FAMILY
from karyoplot.svg.drawing import draw_annotation_track
from karyoplot.svg.legend import draw_grouped_legend
from karyoscope.core.io.features import NOVEL_NAME

from karyoscope_analysis.core.feature_vocab import UnknownFeatureError
from karyoscope_analysis.core.io.bed import Interval

#: Chrome (non-data) text colors. Not feature colors, so not DB-sourced.
_TEXT = "#000000"
_MUTED = "#555555"
_LEGEND_ROW_H = 14  # matches draw_grouped_legend's internal row height


@dataclass(frozen=True)
class PlacedRead:
    """A read placed in the consensus frame: ``segments`` are ``(start, end, feature)`` in
    consensus coordinates (already oriented + positioned)."""

    read_id: str
    is_seed: bool
    reversed: bool
    segments: Sequence[Interval]


@dataclass(frozen=True)
class ClusterPanel:
    """One cluster to draw: a title, its width, placed reads, and its (union) consensus."""

    title: str
    width: int
    placed: Sequence[PlacedRead]
    consensus: Sequence[Interval]


def structural_feature(label: str) -> str:
    """The structural layer of a (possibly composite) label: ``chr13:aSat`` -> ``aSat``."""
    return label.split(":", 1)[1] if ":" in label else label


def feature_color(label: str, colors: Mapping[str, str]) -> str:
    """Color for a label's structural layer, from the DB palette.

    Strict: only ``novel`` may be absent from ``colors`` (it renders white); any other feature
    without a color raises :class:`UnknownFeatureError` (the DB is authoritative).
    """
    feature = structural_feature(label)
    if feature in colors:
        return colors[feature]
    if feature == NOVEL_NAME:
        return "#ffffff"
    raise UnknownFeatureError(feature)


def chromosome_layer(label: str) -> str:
    """The chromosome layer of a composite label: ``chr13:aSat`` -> ``chr13``; ``aSat`` -> ``''``."""
    return label.split(":", 1)[0] if ":" in label else ""


def chromosome_color(label: str, colors: Mapping[str, str]) -> str:
    """Color for a label's chromosome layer, from the DB palette (strict — see :func:`feature_color`)."""
    chrom = chromosome_layer(label)
    if chrom in colors:
        return colors[chrom]
    raise UnknownFeatureError(chrom)


def validate_colors(panels: Sequence[ClusterPanel], colors: Mapping[str, str]) -> None:
    """Raise :class:`UnknownFeatureError` if any segment's feature/chromosome has no DB color."""
    unknown: set[str] = set()
    for panel in panels:
        segs = list(panel.consensus) + [s for r in panel.placed for s in r.segments]
        for _s, _e, label in segs:
            feat = structural_feature(label)
            if feat not in colors and feat != NOVEL_NAME:
                unknown.add(feat)
            chrom = chromosome_layer(label)
            if chrom and chrom not in colors:
                unknown.add(chrom)
    if unknown:
        raise UnknownFeatureError(
            "features have no color in the database palette (only 'novel' may be absent): "
            + ", ".join(sorted(unknown))
        )


def _panel_rows(
    panel: ClusterPanel, *, consensus_track: bool
) -> list[tuple[str, bool, Sequence[Interval]]]:
    rows: list[tuple[str, bool, Sequence[Interval]]] = (
        [("consensus", True, panel.consensus)] if consensus_track else []
    )
    rows += [(r.read_id, r.is_seed, r.segments) for r in panel.placed]
    return rows


def _collect_present(
    panels: Sequence[ClusterPanel],
    colors: Mapping[str, str],
    *,
    chromosome_track: bool,
    consensus_track: bool,
) -> tuple[dict[str, str], dict[str, str]]:
    """Feature -> color and chromosome -> color maps for the legend (across all panels)."""
    present: dict[str, str] = {}
    present_chrom: dict[str, str] = {}
    for panel in panels:
        for _label, _seed, segs in _panel_rows(panel, consensus_track=consensus_track):
            for _s, _e, f in segs:
                present.setdefault(structural_feature(f), feature_color(f, colors))
                if chromosome_track and chromosome_layer(f):
                    present_chrom.setdefault(chromosome_layer(f), chromosome_color(f, colors))
    return present, present_chrom


def _chrom_height(row_height: int, chromosome_track: bool) -> int:
    return max(4, round(row_height * 0.55)) if chromosome_track else 0


def _draw_panel(
    d: draw.Drawing,
    y: float,
    panel: ClusterPanel,
    colors: Mapping[str, str],
    *,
    width: int,
    label_width: int,
    row_height: int,
    chromosome_track: bool,
    consensus_track: bool,
) -> float:
    """Draw one cluster panel (shared consensus x-scale) starting at ``y``; return the next ``y``.

    Each row is the structural-feature track and, directly below it (no gap), a thinner
    chromosome-colored track (when ``chromosome_track``) — so a read's structure and its
    chromosome identity line up. The union consensus is the top row unless ``consensus_track``
    is off (then only the reads are shown).
    """
    plot_width = max(1, width - label_width)
    scale = plot_width / max(1, panel.width)
    chrom_h = _chrom_height(row_height, chromosome_track)

    if panel.title:
        d.append(draw.Text(
            panel.title, 11, 6, y + 11, fill=_TEXT, font_weight="bold",
            font_family=DEFAULT_FONT_FAMILY,
        ))
        y += 17

    for label, is_seed, segs in _panel_rows(panel, consensus_track=consensus_track):
        tag = f"{label}{' (seed)' if is_seed and label != 'consensus' else ''}"
        emphasis = is_seed or label == "consensus"
        d.append(draw.Text(
            tag[:36], 9, 10, y + row_height - 2, fill=_TEXT if emphasis else _MUTED,
            font_family=DEFAULT_FONT_FAMILY,
        ))
        feat_segs = [(s, e, structural_feature(f)) for s, e, f in segs]
        feat_colors = {structural_feature(f): feature_color(f, colors) for _s, _e, f in segs}
        draw_annotation_track(
            d, feat_segs, y, row_height, scale, label_width, "features", "", _TEXT,
            label_width, feat_colors, {}, plot_width, 0, panel.width,
        )
        if chromosome_track:
            chrom_segs = [
                (s, e, chromosome_layer(f)) for s, e, f in segs if chromosome_layer(f)
            ]
            chrom_colors = {
                chromosome_layer(f): chromosome_color(f, colors)
                for _s, _e, f in segs
                if chromosome_layer(f)
            }
            draw_annotation_track(
                d, chrom_segs, y + row_height, chrom_h, scale, label_width, "chromosomes", "",
                _TEXT, label_width, chrom_colors, {}, plot_width, 0, panel.width,
            )
        y += row_height + chrom_h + 2
    return y + 10  # gap after the panel


def _figure_height(
    panels: Sequence[ClusterPanel],
    n_features: int,
    n_chrom: int,
    *,
    row_height: int,
    chromosome_track: bool,
    consensus_track: bool,
) -> float:
    """Total SVG height — computed up front (drawsvg fixes the viewBox at creation time)."""
    chrom_h = _chrom_height(row_height, chromosome_track)
    y: float = 12
    for panel in panels:
        if panel.title:
            y += 17
        n_rows = (1 if consensus_track else 0) + len(panel.placed)
        y += n_rows * (row_height + chrom_h + 2) + 10
    # Legend: a header row + one row per item, per draw_grouped_legend's column layout.
    legend_rows = max(n_features, n_chrom)
    return y + 14 + _LEGEND_ROW_H + legend_rows * _LEGEND_ROW_H + 14


def render_clusters_svg(
    panels: Sequence[ClusterPanel],
    colors: Mapping[str, str],
    *,
    width: int = 1200,
    row_height: int = 11,
    label_width: int = 220,
    chromosome_track: bool = True,
    consensus_track: bool = True,
) -> str:
    """Render one or more cluster panels, stacked, into a single SVG with shared legends."""
    validate_colors(panels, colors)
    present, present_chrom = _collect_present(
        panels, colors, chromosome_track=chromosome_track, consensus_track=consensus_track
    )
    total_h = _figure_height(
        panels, len(present), len(present_chrom),
        row_height=row_height, chromosome_track=chromosome_track, consensus_track=consensus_track,
    )
    d = draw.Drawing(width, total_h, id_prefix="cluster")

    y: float = 12
    for panel in panels:
        y = _draw_panel(
            d, y, panel, colors, width=width, label_width=label_width, row_height=row_height,
            chromosome_track=chromosome_track, consensus_track=consensus_track,
        )

    tracks = ["features"] + (["chromosomes"] if present_chrom else [])
    draw_grouped_legend(
        d, 6, y + 14, _TEXT,
        used_colors={"features": present, "chromosomes": present_chrom},
        track_labels={"features": "Features", "chromosomes": "Chromosomes"},
        tracks=tracks, layout="column", column_width=180,
    )

    svg = d.as_svg()
    return svg[svg.index("<svg"):]  # drop the <?xml?> prolog (keep the string API stable)


def render_cluster_svg(
    placed: Sequence[PlacedRead],
    consensus: Sequence[Interval],
    width: int,
    colors: Mapping[str, str],
    *,
    svg_width: int = 1200,
    row_height: int = 12,
    label_width: int = 220,
    title: str = "",
    chromosome_track: bool = True,
) -> str:
    """Render a single cluster to an SVG (a one-panel :func:`render_clusters_svg`)."""
    return render_clusters_svg(
        [ClusterPanel(title or "cluster", width, placed, consensus)],
        colors,
        width=svg_width,
        row_height=row_height,
        label_width=label_width,
        chromosome_track=chromosome_track,
    )
