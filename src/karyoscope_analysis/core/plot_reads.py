"""Render per-read feature BEDs as stacked colored bars (SVG).

Migrated from the legacy ``KaryoScope_plot_reads.py`` (SVG core). Each read is a bar of
feature-colored runs; reads stack as columns (vertical, the default) or rows (horizontal).
Feature rasterization (bp -> pixel runs) is delegated to
:mod:`karyoplot.svg.reads` (``rasterize_features``); only the figure assembly — column/row
layout, orientation, scale bar, sample headers, and the auto-filtered legend — lives here.

Feature colors come from the DB palette (``colors.tsv`` collapsed to ``{feature: hex}``);
features absent from the palette render white (``#ffffff``), the KaryoScope ``novel`` sentinel
convention. This module is SVG-only (Phase 3a); heatmap/metadata tracks, read-list grouping,
markers (3b) and PNG/animation (3c) are layered on later.
"""

from __future__ import annotations

import fnmatch
import gzip
import re
from collections import defaultdict, namedtuple
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import drawsvg as draw
from karyoplot.core.fonts import DEFAULT_FONT_FAMILY
from karyoplot.svg.reads import rasterize_features

# Feature-name vocabulary used to decide which end of a read goes to the top (position 0).
# These are feature *names*, not colors — orientation is structural, not a palette.
TELOMERE_FEATURES = {"canonical_telomere", "noncanonical_telomere"}
CHROMOSOME_FEATURES = {f"chr{c}" for c in [*range(1, 23), "X", "Y"]} | {
    f"chr{c}_specific" for c in [*range(1, 23), "X", "Y"]
}
SATELLITE_FEATURES = {
    "active", "inactive", "divergent", "monomeric",
    "hsat1A", "hsat1B", "hsat2", "hsat3", "bsat", "gsat", "censat", "noncentromeric",
}

Feature = tuple[int, int, str]  # (start, end, feature_name)


@dataclass(frozen=True)
class Read:
    """One read: feature intervals in bp, plus the sample it belongs to."""

    sample: str
    read_id: str
    length: int
    features: Sequence[Feature]


@dataclass
class PlotConfig:
    """Layout/style knobs for the SVG renderers (Phase 3a subset)."""

    bar_width: int = 5
    read_spacing: int = 5
    sample_spacing: int = 20
    ratio: float = 1 / 300  # bp -> pixels (height in vertical mode, width in horizontal)
    top_margin: int = 30
    left_margin: int = 30
    bottom_margin: int = 15
    background: str = "black"  # "black" or "white"
    font_size: int = 11
    feature_mode: str = "smooth"  # "smooth" | "transition" | "raw"
    min_feature_width: float = 0.5
    min_width_exclude: Sequence[str] = field(default_factory=lambda: ["novel", "*arm*", "ct*"])
    oversample: int = 1
    read_border: bool = False
    draw_scale_bar: bool = True
    no_header: bool = False
    scale_bar_bp: int | None = None
    legend: bool = False

    @property
    def text_color(self) -> str:
        return "#ffffff" if self.background == "black" else "#000000"


# --------------------------------------------------------------------------- loading

def parse_bed_file(bed_path: str | Path, sample: str) -> list[Read]:
    """Parse a per-read feature BED (``read_id  start  end  feature``) into reads."""
    bed_path = str(bed_path)
    read_features: dict[str, list[Feature]] = defaultdict(list)
    open_func = gzip.open if bed_path.endswith(".gz") else open
    mode = "rt" if bed_path.endswith(".gz") else "r"
    with open_func(bed_path, mode) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                read_features[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
    return [
        Read(sample, rid, max(e for _s, e, _f in feats), feats)
        for rid, feats in read_features.items()
    ]


def load_bed_specs(bed_specs: Sequence[str]) -> tuple[list[Read], list[str]]:
    """Load BEDs from ``SAMPLE:PATH`` (or bare ``PATH``) specs; return (reads, sample_order)."""
    all_reads: list[Read] = []
    sample_order: list[str] = []
    for spec in bed_specs:
        if ":" in spec and not spec.startswith("/"):
            sample, bed_path = spec.split(":", 1)
        else:
            bed_path = spec
            sample = Path(bed_path).name.split(".")[0]
        if sample not in sample_order:
            sample_order.append(sample)
        if not Path(bed_path).exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        all_reads.extend(parse_bed_file(bed_path, sample))
    return all_reads, sample_order


def sort_reads(reads: Sequence[Read], sample_order: Sequence[str]) -> list[Read]:
    """Sort by sample order, then by descending length within a sample."""
    rank = {s: i for i, s in enumerate(sample_order)}
    return sorted(reads, key=lambda r: (rank.get(r.sample, 1_000_000), -r.length))


# --------------------------------------------------------------------------- orientation

def _flip(read: Read) -> Read:
    flipped = sorted(
        ((read.length - e, read.length - s, f) for s, e, f in read.features),
        key=lambda x: x[0],
    )
    return Read(read.sample, read.read_id, read.length, flipped)


def _avg_position(read: Read, targets: set[str]) -> float | None:
    positions = [p for s, e, f in read.features if f in targets for p in (s, e)]
    return sum(positions) / len(positions) if positions else None


def orient_reads(reads: Sequence[Read], mode: str) -> list[Read]:
    """Reorient reads so the chosen feature class sits at the top (position 0).

    ``mode`` is one of ``"telomere"``, ``"chromosome"``, or ``"satellite"``. ``"satellite"``
    orients by telomere first, falling back to satellite-dense end for reads with no telomere.
    A read is flipped when its target features average past the read midpoint.
    """
    primary = {"telomere": TELOMERE_FEATURES, "chromosome": CHROMOSOME_FEATURES,
               "satellite": TELOMERE_FEATURES}[mode]
    fallback = SATELLITE_FEATURES if mode == "satellite" else None

    oriented: list[Read] = []
    for read in reads:
        avg = _avg_position(read, primary)
        if avg is None and fallback is not None:
            avg = _avg_position(read, fallback)
        oriented.append(_flip(read) if avg is not None and avg > read.length / 2 else read)
    return oriented


# --------------------------------------------------------------------------- rasterization

def _colored_features(
    features: Sequence[Feature], colors: Mapping[str, str], min_width_exclude: Sequence[str]
) -> list[tuple[int, int, str, float, bool]]:
    """(start, end, color, opacity, skip_min) tuples for ``rasterize_features``."""
    out = []
    for start, end, feature in features:
        skip_min = any(fnmatch.fnmatch(feature, pat) for pat in min_width_exclude)
        out.append((start, end, colors.get(feature, "#ffffff"), 1.0, skip_min))
    return out


def _read_runs(read: Read, colors: Mapping[str, str], cfg: PlotConfig):
    """Rasterized pixel runs for a read, or ``None`` in raw mode (caller draws bp rects)."""
    colored = _colored_features(read.features, colors, cfg.min_width_exclude)
    bar_length = int(max((e for _s, e, _f in read.features), default=0) * cfg.ratio)
    return rasterize_features(
        colored, bar_length, cfg.ratio, cfg.feature_mode, cfg.oversample, cfg.min_feature_width
    )


# --------------------------------------------------------------------------- legend

LegendItem = namedtuple("LegendItem", ["x", "y", "label", "color", "is_header"])
LegendLayout = namedtuple(
    "LegendLayout",
    ["items", "width", "height", "col_width", "cols", "swatch_size", "font_size",
     "row_height", "item_padding"],
)


def filter_legend_features(
    features_used: set[str],
    colors: Mapping[str, str],
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
) -> list[tuple[str, str, str | None]]:
    """(display, color, section_header) for features actually used, dedup'd and cleaned."""
    used_bases = {re.sub(r"_specific$", "", f) for f in features_used}
    seen: set[str] = set()
    result: list[tuple[str, str, str | None]] = []

    def _emit(feat: str, base: str, header: str | None) -> None:
        display = re.sub(r"_multigroup\d+$", "", base).replace("_", " ")
        if display in seen:
            return
        color_hex = colors.get(feat, colors.get(base))
        if color_hex:
            seen.add(display)
            result.append((display, color_hex, header))

    if color_sections:
        for header, section_feats in color_sections:
            for feat in section_feats:
                base = re.sub(r"_specific$", "", feat)
                if base in used_bases:
                    _emit(feat, base, header)
    else:
        for feat in sorted(features_used):
            _emit(feat, re.sub(r"_specific$", "", feat), None)
    return result


def _measure_text(text: str, font_size: int) -> float:
    return len(text) * font_size * 0.48


def compute_legend_layout(
    filtered_items: Sequence[tuple[str, str, str | None]],
    *,
    max_width: int | None = None,
    max_cols: int = 6,
    swatch_size: int = 12,
    font_size: int = 12,
    padding: int = 10,
    col_gap: int = 20,
    item_padding: int = 3,
) -> LegendLayout:
    """Grid layout for legend items, filling top-to-bottom then left-to-right.

    Each section header starts a new column; ``num_rows`` grows until every section fits
    within ``max_cols`` so later sections never overwrite earlier columns.
    """
    if not filtered_items:
        return LegendLayout([], 0, 0, 0, 0, swatch_size, font_size, 0, item_padding)

    swatch_gap = int(swatch_size * 0.4)
    row_height = swatch_size + 2 * item_padding + 4
    max_label_w = max(_measure_text(d, font_size) for d, _, _ in filtered_items)
    col_width = swatch_size + swatch_gap + int(max_label_w) + 2 * item_padding + 4

    if max_width:
        available = max_width - 2 * padding
        num_cols = min(max_cols, max(1, int((available + col_gap) / (col_width + col_gap))))
    else:
        num_cols = min(max_cols, max(1, len(filtered_items)))

    sections: list[tuple[str | None, list[tuple[str, str]]]] = []
    current_header: str | None = None
    current_items: list[tuple[str, str]] = []
    for i, (display, color_hex, header) in enumerate(filtered_items):
        if header != current_header and i != 0:
            sections.append((current_header, current_items))
            current_items = []
        current_header = header
        current_items.append((display, color_hex))
    if current_items:
        sections.append((current_header, current_items))

    total_slots = sum((1 if hdr else 0) + len(its) for hdr, its in sections)
    num_rows = max(1, -(-total_slots // num_cols))

    def _col_demand(rows: int) -> int:
        return sum(max(1, -(-((1 if hdr else 0) + len(its)) // rows)) for hdr, its in sections)

    max_section_slots = max((1 if hdr else 0) + len(its) for hdr, its in sections)
    while _col_demand(num_rows) > num_cols and num_rows < max_section_slots:
        num_rows += 1
    if _col_demand(num_rows) > num_cols:
        num_cols = _col_demand(num_rows)

    positioned: list[LegendItem] = []
    col = row = 0
    for header, items in sections:
        if header and row > 0:
            col += 1
            row = 0
            num_cols = max(num_cols, col + 1)
        if header:
            positioned.append(LegendItem(
                padding + col * (col_width + col_gap), padding + row * row_height,
                header, None, True,
            ))
            row += 1
        for display, color_hex in items:
            if row >= num_rows and col + 1 < num_cols:
                col += 1
                row = 0
            positioned.append(LegendItem(
                padding + col * (col_width + col_gap), padding + row * row_height,
                display, color_hex, False,
            ))
            row += 1

    max_x = max(it.x for it in positioned) + col_width
    max_y = max(it.y for it in positioned) + row_height
    actual_cols = (max(it.x for it in positioned) - padding) // (col_width + col_gap) + 1
    return LegendLayout(positioned, int(max_x + padding), int(max_y + padding), col_width,
                        int(actual_cols), swatch_size, font_size, row_height, item_padding)


def _draw_legend(
    d: draw.Drawing, features_used: set[str], colors: Mapping[str, str], cfg: PlotConfig,
    legend_y: float, image_width: int,
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
) -> float:
    """Draw the auto-filtered feature legend below the reads; return its height."""
    filtered = filter_legend_features(features_used, colors, color_sections)
    if not filtered:
        return 0
    layout = compute_legend_layout(
        filtered, max_width=image_width, swatch_size=cfg.font_size, font_size=cfg.font_size
    )
    for item in layout.items:
        ix, iy = item.x, legend_y + item.y
        if item.is_header:
            d.append(draw.Text(
                item.label, cfg.font_size, ix, iy + cfg.font_size, fill=cfg.text_color,
                text_anchor="start", font_family=DEFAULT_FONT_FAMILY, font_weight="bold",
            ))
            continue
        sx, sy = ix + layout.item_padding, iy + layout.item_padding
        d.append(draw.Rectangle(
            sx, sy, layout.swatch_size, layout.swatch_size, fill=item.color,
            stroke=cfg.text_color, stroke_width=0.5,
        ))
        gap = int(layout.swatch_size * 0.4)
        d.append(draw.Text(
            item.label, cfg.font_size, sx + layout.swatch_size + gap, sy + layout.swatch_size - 1,
            fill=cfg.text_color, text_anchor="start", font_family=DEFAULT_FONT_FAMILY,
        ))
    return layout.height


def _legend_height(
    features_used: set[str], colors: Mapping[str, str], cfg: PlotConfig, image_width: int,
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
) -> float:
    filtered = filter_legend_features(features_used, colors, color_sections)
    if not filtered:
        return 0
    return compute_legend_layout(
        filtered, max_width=image_width, swatch_size=cfg.font_size, font_size=cfg.font_size
    ).height


# --------------------------------------------------------------------------- scale bar

def _draw_vertical_scale_bar(d: draw.Drawing, x: float, y: float, cfg: PlotConfig,
                             max_height_px: int) -> None:
    bp = cfg.scale_bar_bp or 10000
    h = int(bp * cfg.ratio)
    if not cfg.scale_bar_bp and h > max_height_px:
        bp = 5000
        h = int(bp * cfg.ratio)
    d.append(draw.Line(x, y, x, y + h, stroke=cfg.text_color, stroke_width=1))
    d.append(draw.Line(x - 2, y, x + 6, y, stroke=cfg.text_color, stroke_width=1))
    d.append(draw.Line(x - 2, y + h, x + 6, y + h, stroke=cfg.text_color, stroke_width=1))
    lx, ly = x - 5, y + h / 2
    d.append(draw.Text(
        f"{bp // 1000} Kbp", cfg.font_size, lx, ly, fill=cfg.text_color, text_anchor="middle",
        font_family=DEFAULT_FONT_FAMILY, transform=f"rotate(-90, {lx}, {ly})",
    ))


# --------------------------------------------------------------------------- renderers

def _svg_str(d: draw.Drawing) -> str:
    svg = d.as_svg()
    return svg[svg.index("<svg"):]  # drop the <?xml?> prolog (stable string API)


def render_reads_svg(reads: Sequence[Read], colors: Mapping[str, str], cfg: PlotConfig,
                     sample_order: Sequence[str]) -> str:
    """Vertical layout: each read a column of stacked feature runs; scale bar + legend."""
    max_length = max(r.length for r in reads)
    max_height_px = int(max_length * cfg.ratio)
    counts: dict[str, int] = defaultdict(int)
    for r in reads:
        counts[r.sample] += 1

    eff_left = cfg.left_margin if cfg.draw_scale_bar else 20
    image_width = eff_left + len(reads) * (cfg.bar_width + cfg.read_spacing) + \
        (len(counts) - 1) * cfg.sample_spacing + 50
    features_used = {f for r in reads for _s, _e, f in r.features}
    legend_extra = _legend_height(features_used, colors, cfg, image_width) if cfg.legend else 0
    image_height = cfg.top_margin + max_height_px + cfg.bottom_margin + int(legend_extra)

    d = draw.Drawing(image_width, image_height, id_prefix="tr")
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=cfg.background))
    if cfg.draw_scale_bar:
        _draw_vertical_scale_bar(d, cfg.left_margin - 10, cfg.top_margin, cfg, max_height_px)

    x = eff_left
    current_sample: str | None = None
    sample_x_start: dict[str, float] = {}
    sample_x_end: dict[str, float] = {}
    for read in reads:
        if current_sample is not None and read.sample != current_sample:
            sample_x_end[current_sample] = x - cfg.read_spacing
            x += cfg.sample_spacing
        sample_x_start.setdefault(read.sample, x)
        current_sample = read.sample

        if cfg.feature_mode == "raw":
            for start, end, feature in read.features:
                d.append(draw.Rectangle(
                    x, cfg.top_margin + int(start * cfg.ratio), cfg.bar_width,
                    max(1, int((end - start) * cfg.ratio)), fill=colors.get(feature, "#ffffff"),
                ))
        else:
            for run in _read_runs(read, colors, cfg) or []:
                d.append(draw.Rectangle(
                    x, cfg.top_margin + run["scaled_start"], cfg.bar_width,
                    run["scaled_stop"] - run["scaled_start"], fill=run["color"],
                    fill_opacity=run["fill_opacity"],
                ))
        if cfg.read_border:
            h = int(max((e for _s, e, _f in read.features), default=0) * cfg.ratio)
            d.append(draw.Rectangle(x, cfg.top_margin, cfg.bar_width, h, fill="none",
                                    stroke=cfg.text_color, stroke_width=0.5))
        x += cfg.bar_width + cfg.read_spacing
    if current_sample:
        sample_x_end[current_sample] = x - cfg.read_spacing

    if not cfg.no_header:
        for sample in sample_order:
            if sample in sample_x_start and sample in sample_x_end:
                d.append(draw.Line(
                    sample_x_start[sample], cfg.top_margin - 5, sample_x_end[sample],
                    cfg.top_margin - 5, stroke=cfg.text_color, stroke_width=2,
                ))
                d.append(draw.Text(
                    sample.replace("_", " "), cfg.font_size, sample_x_start[sample],
                    cfg.top_margin - 12, fill=cfg.text_color, text_anchor="start",
                    font_family=DEFAULT_FONT_FAMILY,
                ))

    if cfg.legend and features_used:
        _draw_legend(d, features_used, colors, cfg, cfg.top_margin + max_height_px + 10,
                     image_width)
    return _svg_str(d)


def render_reads_horizontal_svg(reads: Sequence[Read], colors: Mapping[str, str],
                                cfg: PlotConfig, sample_order: Sequence[str]) -> str:
    """Horizontal layout: one read per row; horizontal scale bar at the top."""
    max_length = max(r.length for r in reads)
    max_width_px = int(max_length * cfg.ratio)
    counts: dict[str, int] = defaultdict(int)
    for r in reads:
        counts[r.sample] += 1

    image_width = cfg.left_margin + max_width_px + 50
    image_height = cfg.top_margin + len(reads) * (cfg.bar_width + cfg.read_spacing) + \
        (len(counts) - 1) * cfg.sample_spacing + 50

    d = draw.Drawing(image_width, image_height, id_prefix="trh")
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=cfg.background))
    if cfg.draw_scale_bar:
        bp = cfg.scale_bar_bp or 10000
        w = int(bp * cfg.ratio)
        if not cfg.scale_bar_bp and w > max_width_px:
            bp = 5000
            w = int(bp * cfg.ratio)
        sx, sy = cfg.left_margin, cfg.top_margin - 30
        d.append(draw.Rectangle(sx, sy, w, 3, fill=cfg.text_color))
        d.append(draw.Line(sx, sy - 4, sx, sy + 4, stroke=cfg.text_color, stroke_width=1))
        d.append(draw.Line(sx + w, sy - 4, sx + w, sy + 4, stroke=cfg.text_color, stroke_width=1))
        d.append(draw.Text(
            f"{bp // 1000} Kbp", cfg.font_size, sx + w / 2, sy - 8, fill=cfg.text_color,
            text_anchor="middle", font_family=DEFAULT_FONT_FAMILY,
        ))

    y = cfg.top_margin
    current_sample: str | None = None
    sample_y_start: dict[str, float] = {}
    sample_y_end: dict[str, float] = {}
    for read in reads:
        if current_sample is not None and read.sample != current_sample:
            sample_y_end[current_sample] = y
            y += cfg.sample_spacing
        sample_y_start.setdefault(read.sample, y)
        current_sample = read.sample

        if cfg.feature_mode == "raw":
            for start, end, feature in read.features:
                d.append(draw.Rectangle(
                    cfg.left_margin + int(start * cfg.ratio), y,
                    max(1, int((end - start) * cfg.ratio)), cfg.bar_width,
                    fill=colors.get(feature, "#ffffff"),
                ))
        else:
            for run in _read_runs(read, colors, cfg) or []:
                d.append(draw.Rectangle(
                    cfg.left_margin + run["scaled_start"], y,
                    run["scaled_stop"] - run["scaled_start"], cfg.bar_width,
                    fill=run["color"], fill_opacity=run["fill_opacity"],
                ))
        if cfg.read_border:
            w = int(max((e for _s, e, _f in read.features), default=0) * cfg.ratio)
            d.append(draw.Rectangle(cfg.left_margin, y, w, cfg.bar_width, fill="none",
                                    stroke=cfg.text_color, stroke_width=0.5))
        y += cfg.bar_width + cfg.read_spacing
    if current_sample:
        sample_y_end[current_sample] = y

    if not cfg.no_header:
        for sample in sample_order:
            if sample in sample_y_start and sample in sample_y_end:
                line_x = cfg.left_margin - 5
                d.append(draw.Line(line_x, sample_y_start[sample], line_x, sample_y_end[sample],
                                   stroke=cfg.text_color, stroke_width=2))
                lx, ly = line_x - 10, sample_y_start[sample] + 50
                d.append(draw.Text(
                    sample.replace("_", " "), cfg.font_size, lx, ly, fill=cfg.text_color,
                    text_anchor="middle", font_family=DEFAULT_FONT_FAMILY,
                    transform=f"rotate(-90, {lx}, {ly})",
                ))
    return _svg_str(d)


def render(reads: Sequence[Read], colors: Mapping[str, str], cfg: PlotConfig,
           sample_order: Sequence[str], *, horizontal: bool = False) -> str:
    """Render reads to an SVG string (vertical by default, horizontal if requested)."""
    if not reads:
        raise ValueError("no reads to plot")
    fn = render_reads_horizontal_svg if horizontal else render_reads_svg
    return fn(reads, colors, cfg, sample_order)
