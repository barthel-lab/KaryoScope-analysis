"""Render an Engine B cluster as an SVG of feature-colored read tracks.

The single read-renderer for the package (collapsing the legacy ``plot-reads`` /
``cluster-plot`` / ``telogator-reads-viz``): each read is a horizontal row of
feature-colored rectangles, oriented and offset into the cluster seed's coordinate frame
(from `cluster`'s ``layout.tsv``), with the consensus track on top and a feature legend.

Self-contained (emits raw SVG; no plotting deps) so it runs anywhere and is easy to test;
the drawing primitives are the natural thing to push down into ``karyoplot.svg`` later. The
legacy dendrogram / cluster-enrichment visuals are intentionally not carried over (they
belonged to the dropped ``cluster_analysis``); their audits remain for future revival.
"""

from __future__ import annotations

import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from karyoscope_analysis.core.io.bed import Interval

#: Fallback palette for features absent from the colors file (deterministic by name).
_AUTO_PALETTE = (
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
    "#86BCB6",
    "#D37295",
)


@dataclass(frozen=True)
class PlacedRead:
    """A read placed in the seed frame: ``offset`` (bp) + raw per-read intervals."""

    read_id: str
    is_seed: bool
    reversed: bool
    offset: int
    segments: Sequence[Interval]  # (start, end, feature) in the read's own frame


def structural_feature(label: str) -> str:
    """The structural layer of a (possibly composite) label: ``chr13:aSat`` -> ``aSat``."""
    return label.split(":", 1)[1] if ":" in label else label


def feature_color(label: str, colors: Mapping[str, str]) -> str:
    """Color for a label: the colors file (by structural layer), else a stable auto-color."""
    feature = structural_feature(label)
    if feature in colors:
        return colors[feature]
    return _AUTO_PALETTE[zlib.crc32(feature.encode()) % len(_AUTO_PALETTE)]


def _oriented(read: PlacedRead) -> list[Interval]:
    """The read's segments in the seed-facing orientation, shifted to its offset."""
    if not read.reversed:
        local = [(s, e, f) for s, e, f in read.segments]
    else:
        local = []
        pos = 0
        for s, e, f in reversed(list(read.segments)):
            local.append((pos, pos + (e - s), f))
            pos += e - s
    return [(read.offset + s, read.offset + e, f) for s, e, f in local]


def _esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_cluster_svg(
    placed: Sequence[PlacedRead],
    consensus: Sequence[Interval],
    colors: Mapping[str, str],
    *,
    width: int = 1200,
    row_height: int = 12,
    label_width: int = 220,
    title: str = "",
) -> str:
    """Render a cluster (placed reads + consensus) to an SVG string."""
    rows: list[tuple[str, bool, list[Interval]]] = [("consensus", True, list(consensus))]
    rows += [(r.read_id, r.is_seed, _oriented(r)) for r in placed]

    coords = [c for _, _, segs in rows for s, e, _f in segs for c in (s, e)]
    x0 = min(coords) if coords else 0
    span = max(1, (max(coords) if coords else 1) - x0)
    scale = max(1, width - label_width) / span

    elements: list[str] = []
    pad_top = 30 if title else 14
    if title:
        elements.append(
            f'<text x="6" y="18" font-size="13" font-weight="bold">{_esc(title)}</text>'
        )

    present: dict[str, str] = {}
    y = pad_top
    for label, is_seed, segs in rows:
        tag = f"{label}{' (seed)' if is_seed and label != 'consensus' else ''}"
        emphasis = is_seed or label == "consensus"
        elements.append(
            f'<text x="6" y="{y + row_height - 2}" font-size="9" '
            f'fill="{"#000" if emphasis else "#555"}">{_esc(tag[:36])}</text>'
        )
        for s, e, f in segs:
            color = feature_color(f, colors)
            present.setdefault(structural_feature(f), color)
            x = label_width + (s - x0) * scale
            w = max(0.5, (e - s) * scale)
            elements.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{row_height}" fill="{color}" />'
            )
        y += row_height + 2

    legend_y = y + 12
    elements.append(f'<text x="6" y="{legend_y}" font-size="10" font-weight="bold">Features</text>')
    lx, ly = 6, legend_y + 10
    for feature, color in sorted(present.items()):
        if lx + 150 > width:
            lx, ly = 6, ly + 16
        elements.append(f'<rect x="{lx}" y="{ly}" width="11" height="11" fill="{color}" />')
        elements.append(f'<text x="{lx + 15}" y="{ly + 10}" font-size="9">{_esc(feature)}</text>')
        lx += 150
    total_h = ly + 24

    body = "\n".join(elements)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{int(total_h)}" '
        f'font-family="sans-serif">\n{body}\n</svg>\n'
    )
