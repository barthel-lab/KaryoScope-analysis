"""Render per-read feature BEDs as stacked colored bars (SVG).

Migrated from the legacy ``KaryoScope_plot_reads.py`` (SVG core). Each read is a bar of
feature-colored runs; reads stack as columns (vertical, the default) or rows (horizontal).
Feature rasterization (bp -> pixel runs) is delegated to
:mod:`karyoplot.svg.reads` (``rasterize_features``); only the figure assembly — column/row
layout, orientation, scale bar, sample headers, and the auto-filtered legend — lives here.

Feature colors come from the DB palette (``colors.tsv`` collapsed to ``{feature: hex}``). The DB
is authoritative, so only the ``novel`` sentinel may be absent (it renders white); any other
feature without a color is a data/DB mismatch and raises :class:`UnknownFeatureError`. Orientation
feature classes are likewise derived from the DB hierarchy by the caller (no hardcoded biology).
This module is SVG-only; PNG/animation and the telogator preset (3c) are layered on later.
"""

from __future__ import annotations

import fnmatch
import gzip
from collections import defaultdict, namedtuple
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import drawsvg as draw
from karyoplot.core.colors import BARTHEL, qualitative_palette
from karyoplot.core.fonts import DEFAULT_FONT_FAMILY
from karyoplot.svg.reads import rasterize_features
from karyoscope.core.io.features import NOVEL_NAME

from karyoscope_analysis.core.feature_vocab import UnknownFeatureError

#: Neutral fill for heatmap cells whose value is missing (Barthel gray, not a data color).
_HEATMAP_MISSING = BARTHEL["gray"]


def resolve_feature_color(feature: str, colors: Mapping[str, str]) -> str:
    """Color for a feature from the DB palette.

    ``novel`` (the k-mer-not-in-index sentinel) always renders white. Otherwise the feature must
    be in ``colors``; the DB is authoritative (``colors.tsv`` covers every ``hierarchy.tsv``
    node), so a genuinely unknown feature is a data/DB mismatch and raises
    :class:`UnknownFeatureError` rather than rendering as a silent white bar.
    """
    if feature in colors:
        return colors[feature]
    if feature == NOVEL_NAME:
        return "#ffffff"
    raise UnknownFeatureError(feature)


def validate_feature_colors(reads: Sequence[Read], colors: Mapping[str, str]) -> None:
    """Raise :class:`UnknownFeatureError` if any read feature has no DB color (and isn't novel)."""
    unknown = sorted(
        {f for r in reads for _s, _e, f in r.features if f not in colors and f != NOVEL_NAME}
    )
    if unknown:
        raise UnknownFeatureError(
            "features have no color in the database palette (only 'novel' may be absent): "
            + ", ".join(unknown)
        )


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
    """Layout/style knobs for the SVG renderers."""

    bar_width: int = 5
    read_spacing: int = 5
    sample_spacing: int = 20
    subgroup_spacing: int | None = None  # spacing within a group; defaults to sample_spacing
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

    # --- two-tier read-list grouping (3b) ---
    label_tiers: Mapping[str, tuple[str, str | None]] = field(default_factory=dict)
    group_subgroup_order: Sequence[tuple[str, str | None]] = field(default_factory=tuple)
    tier_display_names: Mapping[int, str] = field(default_factory=dict)
    group_boundaries: frozenset[str] = field(default_factory=frozenset)

    # --- heatmap / metadata tracks (3b) ---
    metadata_columns: Sequence[str] = field(default_factory=tuple)
    read_metadata: Mapping[str, Mapping[str, str | None]] = field(default_factory=dict)
    heatmap_colors: Mapping[str, Mapping[str | None, str]] = field(default_factory=dict)
    heatmap_display_names: Mapping[str, str] = field(default_factory=dict)
    heatmap_row_gap: int = 1
    heatmap_top_gap: int = 3
    heatmap_bottom_gap: int = 10

    # --- markers (3b) ---
    markers: Mapping[str, Sequence[tuple[int, int]]] = field(default_factory=dict)
    marker_scale: float = 1.0

    #: Optional feature -> sortable key for ordering the legend (e.g. KaryoScope-style via
    #: core.legend_order.feature_sort_key); None = alphabetical.
    legend_sort_key: object | None = None

    @property
    def text_color(self) -> str:
        return "#ffffff" if self.background == "black" else "#000000"

    @property
    def has_heatmap(self) -> bool:
        return bool(self.metadata_columns and self.read_metadata)

    @property
    def heatmap_total(self) -> int:
        """Vertical space the heatmap rows occupy above the reads (0 if none)."""
        n = len(self.metadata_columns)
        if not (self.has_heatmap and n):
            return 0
        return (
            n * self.bar_width
            + (n - 1) * self.heatmap_row_gap
            + self.heatmap_top_gap
            + self.heatmap_bottom_gap
        )

    @property
    def tier_offset(self) -> int:
        return self.font_size + 10

    @property
    def effective_top_margin(self) -> int:
        """Base top margin plus room reserved for the heatmap and extra label tiers."""
        n_tiers = len(self.tier_display_names)
        extra_tiers = (n_tiers - 1) * self.tier_offset if n_tiers > 1 else 0
        return self.top_margin + self.heatmap_total + extra_tiers


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


# ----------------------------------------------------- read-list, grouping, heatmap, markers


def _composite_label(group: str, subgroup: str | None) -> str:
    """Composite sample label for a (group, subgroup) pair: ``group — subgroup`` or ``group``."""
    return f"{group} — {subgroup}" if subgroup else group


def load_read_list(path: str | Path) -> tuple[set[str], list[str], dict, list[tuple[str, dict]]]:
    """Parse a read-list TSV.

    Returns ``(read_ids, columns, read_data, read_rows)``: the set of read IDs, the column
    names from the header (columns 2+), ``{read_id: {col: value}}`` (last wins), and the rows
    in file order (preserving duplicates so grouping keeps TSV order). A header row is detected
    by a first field of ``sequence``/``read``/``read_id``.
    """
    path = str(path)
    read_ids: set[str] = set()
    columns: list[str] = []
    read_data: dict[str, dict[str, str | None]] = {}
    read_rows: list[tuple[str, dict[str, str | None]]] = []
    open_func = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with open_func(path, mode) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            rid = fields[0]
            if rid.lower() in ("sequence", "read", "read_id"):
                if len(fields) > 1:
                    columns = fields[1:]
                continue
            read_ids.add(rid)
            if columns and len(fields) > 1:
                meta = {
                    c: (fields[1 + i] if 1 + i < len(fields) and fields[1 + i] else None)
                    for i, c in enumerate(columns)
                }
                read_data[rid] = meta
                read_rows.append((rid, meta))
    return read_ids, columns, read_data, read_rows


def build_grouping(
    read_rows: Sequence[tuple[str, Mapping[str, str | None]]],
    tier_specs: Sequence[tuple[str, str]],
) -> tuple[dict[str, tuple[str, str | None]], list[tuple[str, str | None]]]:
    """Build ``{read_id: (group, subgroup)}`` and the ordered ``(group, subgroup)`` pairs.

    ``tier_specs`` is ``[(column, display_name), ...]``; the first column is the group, the
    second (if any) the subgroup. Pairs keep first-seen TSV order, then are stably re-sorted
    so all subgroups of a group stay contiguous.
    """
    if not tier_specs:
        return {}, []
    group_col = tier_specs[0][0]
    subgroup_col = tier_specs[1][0] if len(tier_specs) >= 2 else None

    read_groups: dict[str, tuple[str, str | None]] = {}
    order: list[tuple[str, str | None]] = []
    seen: set[tuple[str, str | None]] = set()
    for rid, meta in read_rows:
        group = meta.get(group_col)
        if not group:
            continue
        subgroup = meta.get(subgroup_col) if subgroup_col else None
        pair = (group, subgroup)
        read_groups[rid] = pair
        if pair not in seen:
            seen.add(pair)
            order.append(pair)

    if subgroup_col:
        group_rank: dict[str, int] = {}
        subgroup_rank: dict[str | None, int] = {}
        for g, s in order:
            group_rank.setdefault(g, len(group_rank))
            subgroup_rank.setdefault(s, len(subgroup_rank))
        order.sort(key=lambda p: (group_rank[p[0]], subgroup_rank[p[1]]))
    return read_groups, order


def apply_grouping(
    reads: Sequence[Read],
    read_groups: Mapping[str, tuple[str, str | None]],
    group_subgroup_order: Sequence[tuple[str, str | None]],
) -> tuple[list[Read], list[str], frozenset[str]]:
    """Relabel reads' ``sample`` to composite group labels; return (reads, order, boundaries).

    ``group_boundaries`` is the set of composite labels that are the last subgroup before a
    group change — drawn with a wider gap / used for the tier-1 spanning lines.
    """
    sample_order = [_composite_label(g, s) for g, s in group_subgroup_order]
    relabeled = [
        Read(_composite_label(*read_groups[r.read_id]), r.read_id, r.length, r.features)
        if r.read_id in read_groups
        else r
        for r in reads
    ]
    boundaries = {
        sample_order[i]
        for i in range(len(group_subgroup_order) - 1)
        if group_subgroup_order[i][0] != group_subgroup_order[i + 1][0]
    }
    return relabeled, sample_order, frozenset(boundaries)


def compute_group_spans(
    sample_order: Sequence[str],
    starts: Mapping[str, float],
    ends: Mapping[str, float],
    group_subgroup_order: Sequence[tuple[str, str | None]],
) -> list[tuple[str, float, float]]:
    """Top-level group spans ``(group, span_start, span_end)`` from per-subgroup extents."""
    result: list[tuple[str, float, float]] = []
    current_group: str | None = None
    span_start: float | None = None
    span_end: float | None = None
    for group, subgroup in group_subgroup_order:
        label = _composite_label(group, subgroup)
        if label not in starts:
            continue
        if group != current_group:
            if current_group is not None and span_start is not None:
                result.append((current_group, span_start, span_end))
            current_group = group
            span_start = starts[label]
            span_end = ends.get(label, span_start)
        else:
            span_end = ends.get(label, span_end)
    if current_group is not None and span_start is not None:
        result.append((current_group, span_start, span_end))
    return result


def assign_heatmap_colors(
    metadata_columns: Sequence[str],
    read_metadata: Mapping[str, Mapping[str, str | None]],
) -> dict[str, dict[str | None, str]]:
    """Per-column ``{value: hex}`` maps, coloring unique values from the shared TAB20 palette.

    Categorical metadata colors (not biological features) come from the single library palette
    (:func:`karyoplot.core.colors.qualitative_palette`), never a per-script literal list. Missing
    values map to a neutral gray.
    """
    color_map: dict[str, dict[str | None, str]] = {}
    for col in metadata_columns:
        seen: list[str] = []
        for meta in read_metadata.values():
            val = meta.get(col)
            if val is not None and val not in seen:
                seen.append(val)
        palette = qualitative_palette(len(seen)) if seen else []
        val_colors: dict[str | None, str] = {val: palette[i] for i, val in enumerate(seen)}
        val_colors[None] = _HEATMAP_MISSING
        color_map[col] = val_colors
    return color_map


def build_heatmap_legend_items(
    metadata_columns: Sequence[str],
    heatmap_colors: Mapping[str, Mapping[str | None, str]],
    heatmap_display_names: Mapping[str, str],
) -> list[tuple[str, str, str | None]]:
    """Heatmap value swatches as ``(display, color, section_header)`` legend rows."""
    items: list[tuple[str, str, str | None]] = []
    for col in metadata_columns:
        header = heatmap_display_names.get(col, col)
        for val, color in heatmap_colors.get(col, {}).items():
            if val is not None:
                items.append((str(val).replace("_", " "), color, header))
    return items


def parse_markers(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Parse a markers TSV (``read_id  start  end``) into ``{read_id: [(start, end), ...]}``."""
    markers: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with Path(path).open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0].lower() not in ("read_id", "read", "sequence"):
                try:
                    markers[parts[0]].append((int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    return dict(markers)


@dataclass
class ReadListData:
    """Grouping + heatmap pieces derived from a read-list TSV (see :func:`process_read_list`)."""

    reads: list[Read]
    read_groups: dict[str, tuple[str, str | None]]
    group_subgroup_order: list[tuple[str, str | None]]
    metadata_columns: list[str]
    read_metadata: dict[str, dict[str, str | None]]
    tier_specs: list[tuple[str, str]]
    heatmap_specs: list[tuple[str, str]]


def _parse_specs(raw: Sequence[str]) -> list[tuple[str, str]]:
    """Parse ``COLUMN[:DISPLAY]`` specs into ``(column, display_name)`` pairs."""
    out = []
    for spec in raw:
        col, _, name = spec.partition(":")
        out.append((col, name or col))
    return out


def process_read_list(
    reads: Sequence[Read],
    path: str | Path,
    *,
    label_tier_specs: Sequence[str] = (),
    heatmap_track_specs: Sequence[str] = (),
    heatmap: bool = False,
    filter_groups: Sequence[str] = (),
) -> ReadListData:
    """Filter reads to a read-list and derive its grouping + heatmap structure.

    Tiers default to the first two TSV data columns; heatmap tracks default to the remaining
    columns when ``heatmap`` is set. ``filter_groups`` keeps only reads whose group (the first
    tier column) is in the set. Mirrors the legacy ``--read-list`` semantics.
    """
    read_ids, columns, read_data, read_rows = load_read_list(path)
    kept = [r for r in reads if r.read_id in read_ids]

    tier_specs = _parse_specs(label_tier_specs)
    if not tier_specs and len(columns) >= 2:
        tier_specs = [(columns[0], columns[0]), (columns[1], columns[1])]

    heatmap_specs = _parse_specs(heatmap_track_specs)
    if not heatmap_specs and heatmap and columns:
        tier_cols = {c for c, _ in tier_specs}
        heatmap_specs = [(c, c) for c in columns if c not in tier_cols]

    read_groups, group_subgroup_order = build_grouping(read_rows, tier_specs)

    if filter_groups and tier_specs:
        allowed = set(filter_groups)
        group_col = tier_specs[0][0]
        sub_col = tier_specs[1][0] if len(tier_specs) >= 2 else None
        keep_ids: set[str] = set()
        read_groups, group_subgroup_order, seen = {}, [], set()
        for rid, meta in read_rows:
            g = meta.get(group_col)
            if g not in allowed:
                continue
            keep_ids.add(rid)
            pair = (g, meta.get(sub_col) if sub_col else None)
            read_groups[rid] = pair
            if pair not in seen:
                seen.add(pair)
                group_subgroup_order.append(pair)
        kept = [r for r in kept if r.read_id in keep_ids]

    metadata_columns = [c for c, _ in heatmap_specs]
    read_metadata = {
        rid: {c: read_data.get(rid, {}).get(c) for c in metadata_columns} for rid in read_ids
    }
    return ReadListData(
        kept,
        read_groups,
        group_subgroup_order,
        metadata_columns,
        read_metadata,
        tier_specs,
        heatmap_specs,
    )


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


def orient_reads(
    reads: Sequence[Read],
    primary_features: set[str],
    fallback_features: set[str] | None = None,
) -> list[Read]:
    """Reorient reads so the given feature class sits at the top (position 0).

    A read is flipped when its ``primary_features`` average past the read midpoint; reads with
    no primary feature fall back to ``fallback_features`` if given. The feature *sets* are
    supplied by the caller (derived from the DB hierarchy, e.g.
    ``FeatureHierarchy.telomere_features`` / ``.chromosomes`` / ``.satellite_features``), so this
    stays free of hardcoded biology.
    """
    oriented: list[Read] = []
    for read in reads:
        avg = _avg_position(read, primary_features)
        if avg is None and fallback_features is not None:
            avg = _avg_position(read, fallback_features)
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
        out.append((start, end, resolve_feature_color(feature, colors), 1.0, skip_min))
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
    [
        "items",
        "width",
        "height",
        "col_width",
        "cols",
        "swatch_size",
        "font_size",
        "row_height",
        "item_padding",
    ],
)


def filter_legend_features(
    features_used: set[str],
    colors: Mapping[str, str],
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
    sort_key=None,
) -> list[tuple[str, str, str | None]]:
    """(display, color, section_header) for features actually used, dedup'd; ``_`` -> space.

    Without ``color_sections``, features are ordered by ``sort_key`` (on the feature name) if
    given, else alphabetically.
    """
    seen: set[str] = set()
    result: list[tuple[str, str, str | None]] = []

    def _emit(feat: str, header: str | None) -> None:
        display = feat.replace("_", " ")
        if display in seen:
            return
        color_hex = colors.get(feat)
        if color_hex:
            seen.add(display)
            result.append((display, color_hex, header))

    if color_sections:
        for header, section_feats in color_sections:
            for feat in section_feats:
                if feat in features_used:
                    _emit(feat, header)
    else:
        ordered = sorted(features_used, key=sort_key) if sort_key else sorted(features_used)
        for feat in ordered:
            _emit(feat, None)
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
            positioned.append(
                LegendItem(
                    padding + col * (col_width + col_gap),
                    padding + row * row_height,
                    header,
                    None,
                    True,
                )
            )
            row += 1
        for display, color_hex in items:
            if row >= num_rows and col + 1 < num_cols:
                col += 1
                row = 0
            positioned.append(
                LegendItem(
                    padding + col * (col_width + col_gap),
                    padding + row * row_height,
                    display,
                    color_hex,
                    False,
                )
            )
            row += 1

    max_x = max(it.x for it in positioned) + col_width
    max_y = max(it.y for it in positioned) + row_height
    actual_cols = (max(it.x for it in positioned) - padding) // (col_width + col_gap) + 1
    return LegendLayout(
        positioned,
        int(max_x + padding),
        int(max_y + padding),
        col_width,
        int(actual_cols),
        swatch_size,
        font_size,
        row_height,
        item_padding,
    )


def _draw_legend(
    d: draw.Drawing,
    features_used: set[str],
    colors: Mapping[str, str],
    cfg: PlotConfig,
    legend_y: float,
    image_width: int,
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
    extra_items: Sequence[tuple[str, str, str | None]] | None = None,
) -> float:
    """Draw the legend (optional heatmap items prepended, then features); return its height."""
    filtered = list(extra_items or []) + filter_legend_features(
        features_used, colors, color_sections, cfg.legend_sort_key
    )
    if not filtered:
        return 0
    layout = compute_legend_layout(
        filtered, max_width=image_width, swatch_size=cfg.font_size, font_size=cfg.font_size
    )
    for item in layout.items:
        ix, iy = item.x, legend_y + item.y
        if item.is_header:
            d.append(
                draw.Text(
                    item.label,
                    cfg.font_size,
                    ix,
                    iy + cfg.font_size,
                    fill=cfg.text_color,
                    text_anchor="start",
                    font_family=DEFAULT_FONT_FAMILY,
                    font_weight="bold",
                )
            )
            continue
        sx, sy = ix + layout.item_padding, iy + layout.item_padding
        d.append(
            draw.Rectangle(
                sx,
                sy,
                layout.swatch_size,
                layout.swatch_size,
                fill=item.color,
                stroke=cfg.text_color,
                stroke_width=0.5,
            )
        )
        gap = int(layout.swatch_size * 0.4)
        d.append(
            draw.Text(
                item.label,
                cfg.font_size,
                sx + layout.swatch_size + gap,
                sy + layout.swatch_size - 1,
                fill=cfg.text_color,
                text_anchor="start",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )
    return layout.height


def _legend_height(
    features_used: set[str],
    colors: Mapping[str, str],
    cfg: PlotConfig,
    image_width: int,
    color_sections: Sequence[tuple[str | None, Sequence[str]]] | None = None,
    extra_items: Sequence[tuple[str, str, str | None]] | None = None,
) -> float:
    filtered = list(extra_items or []) + filter_legend_features(
        features_used, colors, color_sections, cfg.legend_sort_key
    )
    if not filtered:
        return 0
    return compute_legend_layout(
        filtered, max_width=image_width, swatch_size=cfg.font_size, font_size=cfg.font_size
    ).height


# --------------------------------------------------------------------------- scale bar


def _draw_vertical_scale_bar(
    d: draw.Drawing, x: float, y: float, cfg: PlotConfig, max_height_px: int
) -> None:
    bp = cfg.scale_bar_bp or 10000
    h = int(bp * cfg.ratio)
    if not cfg.scale_bar_bp and h > max_height_px:
        bp = 5000
        h = int(bp * cfg.ratio)
    d.append(draw.Line(x, y, x, y + h, stroke=cfg.text_color, stroke_width=1))
    d.append(draw.Line(x - 2, y, x + 6, y, stroke=cfg.text_color, stroke_width=1))
    d.append(draw.Line(x - 2, y + h, x + 6, y + h, stroke=cfg.text_color, stroke_width=1))
    lx, ly = x - 5, y + h / 2
    d.append(
        draw.Text(
            f"{bp // 1000} Kbp",
            cfg.font_size,
            lx,
            ly,
            fill=cfg.text_color,
            text_anchor="middle",
            font_family=DEFAULT_FONT_FAMILY,
            transform=f"rotate(-90, {lx}, {ly})",
        )
    )


# --------------------------------------------------------------------------- heatmap / markers / headers


def _effective_left_margin(cfg: PlotConfig, base_left: int) -> int:
    """Grow the left margin to fit tier / heatmap row labels drawn in the margin."""
    labels = [dn for dn in cfg.tier_display_names.values()]
    labels += [cfg.heatmap_display_names.get(c, c) for c in cfg.metadata_columns]
    if not labels:
        return base_left
    needed = max(60, int(max(len(s) for s in labels) * 8.5) + 10)
    return max(base_left, needed)


def _draw_heatmap_grid(
    d: draw.Drawing,
    read_positions: Sequence[tuple[str, float]],
    cfg: PlotConfig,
    top: int,
    left: int,
) -> None:
    """Draw the metadata heatmap: one row of colored cells per column, stacked above the reads."""
    stroke = cfg.text_color
    for row_idx, col in enumerate(reversed(list(cfg.metadata_columns))):
        row_y = (
            top
            - cfg.heatmap_bottom_gap
            - (row_idx + 1) * cfg.bar_width
            - row_idx * cfg.heatmap_row_gap
        )
        val_colors = cfg.heatmap_colors.get(col, {})
        for read_id, rx in read_positions:
            val = cfg.read_metadata.get(read_id, {}).get(col)
            color = val_colors.get(val, val_colors.get(None, _HEATMAP_MISSING))
            d.append(
                draw.Rectangle(
                    rx,
                    row_y,
                    cfg.bar_width,
                    cfg.bar_width,
                    fill=color,
                    stroke=stroke,
                    stroke_width=0.5,
                )
            )
        d.append(
            draw.Text(
                cfg.heatmap_display_names.get(col, col),
                cfg.font_size,
                left - 5,
                row_y + cfg.bar_width / 2 + cfg.font_size * 0.35,
                fill=cfg.text_color,
                text_anchor="end",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )


def _draw_markers(
    d: draw.Drawing, positions: Sequence[tuple[float, float]], cfg: PlotConfig, arrow_size: int
) -> None:
    """Draw left-pointing arrowheads at the given (x, mid_y) read positions."""
    for x, mid_y in positions:
        d.append(
            draw.Lines(
                x - arrow_size - 1,
                mid_y - arrow_size,
                x - 1,
                mid_y,
                x - arrow_size - 1,
                mid_y + arrow_size,
                fill=cfg.text_color,
                close=True,
            )
        )


def _draw_vertical_header(
    d: draw.Drawing,
    cfg: PlotConfig,
    top: int,
    left: int,
    sample_order: Sequence[str],
    x_start: Mapping[str, float],
    x_end: Mapping[str, float],
) -> None:
    """Draw sample separator lines + labels above the reads (single- or two-tier)."""
    hm = cfg.heatmap_total
    if not cfg.label_tiers:
        for sample in sample_order:
            if sample in x_start and sample in x_end:
                d.append(
                    draw.Line(
                        x_start[sample],
                        top - 5 - hm,
                        x_end[sample],
                        top - 5 - hm,
                        stroke=cfg.text_color,
                        stroke_width=2,
                    )
                )
                d.append(
                    draw.Text(
                        sample.replace("_", " "),
                        cfg.font_size,
                        x_start[sample],
                        top - 12 - hm,
                        fill=cfg.text_color,
                        text_anchor="start",
                        font_family=DEFAULT_FONT_FAMILY,
                    )
                )
        return

    # Tier 2 (inner): subgroup lines + labels.
    tier2_y = top - 5 - hm
    for sample in sample_order:
        if sample in x_start and sample in x_end:
            d.append(
                draw.Line(
                    x_start[sample],
                    tier2_y,
                    x_end[sample],
                    tier2_y,
                    stroke=cfg.text_color,
                    stroke_width=1.5,
                )
            )
            _g, subgroup = cfg.label_tiers.get(sample, (sample, None))
            d.append(
                draw.Text(
                    (subgroup or sample).replace("_", " "),
                    cfg.font_size,
                    x_start[sample],
                    top - 10 - hm,
                    fill=cfg.text_color,
                    text_anchor="start",
                    font_family=DEFAULT_FONT_FAMILY,
                )
            )
    if 1 in cfg.tier_display_names:
        d.append(
            draw.Text(
                cfg.tier_display_names[1],
                cfg.font_size,
                left - 5,
                top - 10 - hm,
                fill=cfg.text_color,
                text_anchor="end",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )

    # Tier 1 (outer): group spanning lines + labels.
    tier1_y = top - 5 - cfg.tier_offset - hm
    for group, gx_start, gx_end in compute_group_spans(
        sample_order, x_start, x_end, cfg.group_subgroup_order
    ):
        d.append(
            draw.Line(gx_start, tier1_y, gx_end, tier1_y, stroke=cfg.text_color, stroke_width=1.5)
        )
        d.append(
            draw.Text(
                group.replace("_", " "),
                cfg.font_size,
                gx_start,
                tier1_y - 5,
                fill=cfg.text_color,
                text_anchor="start",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )
    if 0 in cfg.tier_display_names:
        d.append(
            draw.Text(
                cfg.tier_display_names[0],
                cfg.font_size,
                left - 5,
                tier1_y - 5,
                fill=cfg.text_color,
                text_anchor="end",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )


# --------------------------------------------------------------------------- renderers


def _svg_str(d: draw.Drawing) -> str:
    svg = d.as_svg()
    return svg[svg.index("<svg") :]  # drop the <?xml?> prolog (stable string API)


def render_reads_svg(
    reads: Sequence[Read], colors: Mapping[str, str], cfg: PlotConfig, sample_order: Sequence[str]
) -> str:
    """Vertical layout: each read a column of stacked feature runs.

    Reserves room above the reads for the heatmap rows and extra label tiers; draws the scale
    bar, optional heatmap grid, sample headers (single- or two-tier), markers, and legend.
    """
    top = cfg.effective_top_margin
    subgroup_spacing = (
        cfg.subgroup_spacing if cfg.subgroup_spacing is not None else cfg.sample_spacing
    )
    arrow_size = max(2, int(cfg.bar_width // 3 * cfg.marker_scale))
    max_length = max(r.length for r in reads)
    max_height_px = int(max_length * cfg.ratio)
    counts: dict[str, int] = defaultdict(int)
    for r in reads:
        counts[r.sample] += 1

    base_left = cfg.left_margin if cfg.draw_scale_bar else 20
    left = _effective_left_margin(cfg, base_left)
    if cfg.markers:
        left += arrow_size + 2  # room for arrowheads
    image_width = (
        left
        + len(reads) * (cfg.bar_width + cfg.read_spacing)
        + (len(counts) - 1) * max(cfg.sample_spacing, subgroup_spacing)
        + 50
    )

    features_used = {f for r in reads for _s, _e, f in r.features}
    hm_items = (
        build_heatmap_legend_items(
            cfg.metadata_columns, cfg.heatmap_colors, cfg.heatmap_display_names
        )
        if cfg.has_heatmap
        else None
    )
    legend_extra = (
        _legend_height(features_used, colors, cfg, image_width, extra_items=hm_items)
        if cfg.legend
        else 0
    )
    image_height = top + max_height_px + cfg.bottom_margin + int(legend_extra)

    d = draw.Drawing(image_width, image_height, id_prefix="tr")
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=cfg.background))
    if cfg.draw_scale_bar:
        _draw_vertical_scale_bar(d, base_left - 10, top, cfg, max_height_px)

    x = left
    current_sample: str | None = None
    x_start: dict[str, float] = {}
    x_end: dict[str, float] = {}
    read_positions: list[tuple[str, float]] = []
    marker_positions: list[tuple[float, float]] = []
    for read in reads:
        if current_sample is not None and read.sample != current_sample:
            x_end[current_sample] = x - cfg.read_spacing
            x += cfg.sample_spacing if current_sample in cfg.group_boundaries else subgroup_spacing
        x_start.setdefault(read.sample, x)
        current_sample = read.sample
        read_positions.append((read.read_id, x))

        if cfg.feature_mode == "raw":
            for start, end, feature in read.features:
                d.append(
                    draw.Rectangle(
                        x,
                        top + int(start * cfg.ratio),
                        cfg.bar_width,
                        max(1, int((end - start) * cfg.ratio)),
                        fill=resolve_feature_color(feature, colors),
                    )
                )
        else:
            for run in _read_runs(read, colors, cfg) or []:
                d.append(
                    draw.Rectangle(
                        x,
                        top + run["scaled_start"],
                        cfg.bar_width,
                        run["scaled_stop"] - run["scaled_start"],
                        fill=run["color"],
                        fill_opacity=run["fill_opacity"],
                    )
                )
        if cfg.read_border:
            h = int(max((e for _s, e, _f in read.features), default=0) * cfg.ratio)
            d.append(
                draw.Rectangle(
                    x, top, cfg.bar_width, h, fill="none", stroke=cfg.text_color, stroke_width=0.5
                )
            )
        for m_start, m_end in cfg.markers.get(read.read_id, []):
            marker_positions.append((x, top + int((m_start + m_end) / 2 * cfg.ratio)))
        x += cfg.bar_width + cfg.read_spacing
    if current_sample:
        x_end[current_sample] = x - cfg.read_spacing

    if cfg.has_heatmap:
        _draw_heatmap_grid(d, read_positions, cfg, top, left)
    _draw_markers(d, marker_positions, cfg, arrow_size)
    if not cfg.no_header:
        _draw_vertical_header(d, cfg, top, left, sample_order, x_start, x_end)

    if cfg.legend and (features_used or hm_items):
        _draw_legend(
            d,
            features_used,
            colors,
            cfg,
            top + max_height_px + 10,
            image_width,
            extra_items=hm_items,
        )
    return _svg_str(d)


def render_reads_horizontal_svg(
    reads: Sequence[Read], colors: Mapping[str, str], cfg: PlotConfig, sample_order: Sequence[str]
) -> str:
    """Horizontal layout: one read per row; horizontal scale bar at the top."""
    max_length = max(r.length for r in reads)
    max_width_px = int(max_length * cfg.ratio)
    counts: dict[str, int] = defaultdict(int)
    for r in reads:
        counts[r.sample] += 1

    image_width = cfg.left_margin + max_width_px + 50
    image_height = (
        cfg.top_margin
        + len(reads) * (cfg.bar_width + cfg.read_spacing)
        + (len(counts) - 1) * cfg.sample_spacing
        + 50
    )

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
        d.append(
            draw.Text(
                f"{bp // 1000} Kbp",
                cfg.font_size,
                sx + w / 2,
                sy - 8,
                fill=cfg.text_color,
                text_anchor="middle",
                font_family=DEFAULT_FONT_FAMILY,
            )
        )

    subgroup_spacing = (
        cfg.subgroup_spacing if cfg.subgroup_spacing is not None else cfg.sample_spacing
    )
    y = cfg.top_margin
    current_sample: str | None = None
    sample_y_start: dict[str, float] = {}
    sample_y_end: dict[str, float] = {}
    for read in reads:
        if current_sample is not None and read.sample != current_sample:
            sample_y_end[current_sample] = y
            y += cfg.sample_spacing if current_sample in cfg.group_boundaries else subgroup_spacing
        sample_y_start.setdefault(read.sample, y)
        current_sample = read.sample

        if cfg.feature_mode == "raw":
            for start, end, feature in read.features:
                d.append(
                    draw.Rectangle(
                        cfg.left_margin + int(start * cfg.ratio),
                        y,
                        max(1, int((end - start) * cfg.ratio)),
                        cfg.bar_width,
                        fill=resolve_feature_color(feature, colors),
                    )
                )
        else:
            for run in _read_runs(read, colors, cfg) or []:
                d.append(
                    draw.Rectangle(
                        cfg.left_margin + run["scaled_start"],
                        y,
                        run["scaled_stop"] - run["scaled_start"],
                        cfg.bar_width,
                        fill=run["color"],
                        fill_opacity=run["fill_opacity"],
                    )
                )
        if cfg.read_border:
            w = int(max((e for _s, e, _f in read.features), default=0) * cfg.ratio)
            d.append(
                draw.Rectangle(
                    cfg.left_margin,
                    y,
                    w,
                    cfg.bar_width,
                    fill="none",
                    stroke=cfg.text_color,
                    stroke_width=0.5,
                )
            )
        y += cfg.bar_width + cfg.read_spacing
    if current_sample:
        sample_y_end[current_sample] = y

    if not cfg.no_header:
        for sample in sample_order:
            if sample in sample_y_start and sample in sample_y_end:
                line_x = cfg.left_margin - 5
                d.append(
                    draw.Line(
                        line_x,
                        sample_y_start[sample],
                        line_x,
                        sample_y_end[sample],
                        stroke=cfg.text_color,
                        stroke_width=2,
                    )
                )
                lx, ly = line_x - 10, sample_y_start[sample] + 50
                d.append(
                    draw.Text(
                        sample.replace("_", " "),
                        cfg.font_size,
                        lx,
                        ly,
                        fill=cfg.text_color,
                        text_anchor="middle",
                        font_family=DEFAULT_FONT_FAMILY,
                        transform=f"rotate(-90, {lx}, {ly})",
                    )
                )
    return _svg_str(d)


def render(
    reads: Sequence[Read],
    colors: Mapping[str, str],
    cfg: PlotConfig,
    sample_order: Sequence[str],
    *,
    horizontal: bool = False,
) -> str:
    """Render reads to an SVG string (vertical by default, horizontal if requested)."""
    if not reads:
        raise ValueError("no reads to plot")
    validate_feature_colors(reads, colors)
    fn = render_reads_horizontal_svg if horizontal else render_reads_svg
    return fn(reads, colors, cfg, sample_order)
