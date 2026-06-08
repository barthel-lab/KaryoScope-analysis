"""Render Engine B clusters as SVG read tracks (consensus-coordinate layout).

The single read-renderer for the package. Each read is a horizontal row of feature-colored
rectangles placed in the cluster's **consensus coordinate frame** (computed by
`feature_assembly.consensus_layout` and serialized to `cluster`'s ``layout.tsv``), so matched
features **stack vertically** across reads; the union-spanning consensus track sits on top. A
gap in a row means that read has no feature there (it didn't align / doesn't extend that far) —
length filtering is the caller's choice, off by default, so gaps are unambiguous.

Self-contained (emits raw SVG; no plotting deps). One cluster (:func:`render_cluster_svg`) or
many stacked in one figure (:func:`render_clusters_svg`), sharing a single feature legend.
"""

from __future__ import annotations

import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from karyoscope_analysis.core.io.bed import Interval

#: Fallback palette for features absent from the colors file (deterministic by name).
_AUTO_PALETTE = (
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
    "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
)


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
    """Color for a label: the colors file (by structural layer), else a stable auto-color."""
    feature = structural_feature(label)
    if feature in colors:
        return colors[feature]
    return _AUTO_PALETTE[zlib.crc32(feature.encode()) % len(_AUTO_PALETTE)]


def _esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _draw_panel(
    elements: list[str],
    y: float,
    panel: ClusterPanel,
    colors: Mapping[str, str],
    present: dict[str, str],
    *,
    width: int,
    label_width: int,
    row_height: int,
) -> float:
    """Draw one cluster panel (shared consensus x-scale) starting at ``y``; return the next ``y``."""
    rows: list[tuple[str, bool, Sequence[Interval]]] = [("consensus", True, panel.consensus)]
    rows += [(r.read_id, r.is_seed, r.segments) for r in panel.placed]
    scale = max(1, width - label_width) / max(1, panel.width)

    if panel.title:
        elements.append(
            f'<text x="6" y="{y + 11:.0f}" font-size="11" font-weight="bold">'
            f"{_esc(panel.title)}</text>"
        )
        y += 17

    for label, is_seed, segs in rows:
        tag = f"{label}{' (seed)' if is_seed and label != 'consensus' else ''}"
        emphasis = is_seed or label == "consensus"
        elements.append(
            f'<text x="10" y="{y + row_height - 2:.0f}" font-size="9" '
            f'fill="{"#000" if emphasis else "#555"}">{_esc(tag[:36])}</text>'
        )
        for s, e, f in segs:
            color = feature_color(f, colors)
            present.setdefault(structural_feature(f), color)
            x = label_width + s * scale
            w = max(0.5, (e - s) * scale)
            elements.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{row_height}" fill="{color}" />'
            )
        y += row_height + 2
    return y + 10  # gap after the panel


def _draw_legend(elements: list[str], y: float, present: Mapping[str, str], width: int) -> float:
    elements.append(f'<text x="6" y="{y:.0f}" font-size="10" font-weight="bold">Features</text>')
    lx, ly = 6, y + 10
    for feature, color in sorted(present.items()):
        if lx + 150 > width:
            lx, ly = 6, ly + 16
        elements.append(f'<rect x="{lx}" y="{ly:.0f}" width="11" height="11" fill="{color}" />')
        elements.append(
            f'<text x="{lx + 15}" y="{ly + 10:.0f}" font-size="9">{_esc(feature)}</text>'
        )
        lx += 150
    return ly + 24


def render_clusters_svg(
    panels: Sequence[ClusterPanel],
    colors: Mapping[str, str],
    *,
    width: int = 1200,
    row_height: int = 11,
    label_width: int = 220,
) -> str:
    """Render one or more cluster panels, stacked, into a single SVG with a shared legend."""
    elements: list[str] = []
    present: dict[str, str] = {}
    y: float = 12
    for panel in panels:
        y = _draw_panel(
            elements, y, panel, colors, present,
            width=width, label_width=label_width, row_height=row_height,
        )
    total_h = _draw_legend(elements, y + 4, present, width)
    body = "\n".join(elements)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{int(total_h)}" '
        f'font-family="sans-serif">\n{body}\n</svg>\n'
    )


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
) -> str:
    """Render a single cluster to an SVG (a one-panel :func:`render_clusters_svg`)."""
    return render_clusters_svg(
        [ClusterPanel(title or "cluster", width, placed, consensus)],
        colors,
        width=svg_width,
        row_height=row_height,
        label_width=label_width,
    )
