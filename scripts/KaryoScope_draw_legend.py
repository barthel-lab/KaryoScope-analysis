#!/usr/bin/env python3
"""
KaryoScope_draw_legend.py

Generate SVG legends from KaryoScope color mapping files.

Usage:
    # Basic: auto-layout all features in 4 columns
    python KaryoScope_draw_legend.py \
        --colors KS_human_CHM13.chromosome_acrocentric.colors.txt \
        --output legend.svg --columns 4

    # Merge same-color features with custom labels
    python KaryoScope_draw_legend.py \
        --colors KS_human_CHM13.chromosome_acrocentric.colors.txt \
        --output legend.svg --columns 4 \
        --merge-same-color \
        --merge-label "chr13_specific=acrocentric,autosome_multigroup1=categorized"

    # With category headers and filtering
    python KaryoScope_draw_legend.py \
        --colors KS_human_CHM13.chromosome_acrocentric.colors.txt \
        --output legend.svg --columns 3 \
        --exclude "categorized" \
        --groups "Chromosomes:chr1,chr2,chr3;Features:DJ,PJ,rDNA"
"""

import argparse
import math
import os
import sys
from collections import OrderedDict

import drawsvg as draw


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate SVG legends from KaryoScope color mapping files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--colors",
        required=True,
        help="Path to TSV colors file (feature, color columns)",
    )
    parser.add_argument(
        "--output",
        default="legend.svg",
        help="Output SVG path (default: legend.svg)",
    )

    # Layout
    parser.add_argument(
        "--columns",
        type=int,
        default=None,
        help="Number of columns (auto-calculated if not set)",
    )
    parser.add_argument(
        "--rows",
        type=int,
        default=None,
        help="Number of rows (auto-calculated if not set)",
    )

    # Filtering
    parser.add_argument(
        "--include",
        default=None,
        help="Comma-separated feature names to include (only these are shown)",
    )
    parser.add_argument(
        "--exclude",
        default=None,
        help="Comma-separated feature names to exclude",
    )

    # Merging
    parser.add_argument(
        "--merge-same-color",
        action="store_true",
        help="Collapse features sharing the same hex color into one entry",
    )
    parser.add_argument(
        "--merge-label",
        default=None,
        help="Override merged labels: feature=label,feature=label,...\n"
             "e.g. 'chr13_specific=acrocentric,autosome_multigroup1=categorized'",
    )

    # Grouping
    parser.add_argument(
        "--groups",
        default=None,
        help="Category headers: 'Header1:feat1,feat2;Header2:feat3,feat4'\n"
             "Features matched by prefix (e.g. 'chr' matches chr1, chr2, ...)",
    )

    # Styling
    parser.add_argument(
        "--swatch-size", type=int, default=8,
        help="Swatch width/height in px (default: 8)",
    )
    parser.add_argument(
        "--font-size", type=int, default=12,
        help="Label font size in px (default: 12)",
    )
    parser.add_argument(
        "--background", default="#000000",
        help="Background color (default: #000000)",
    )
    parser.add_argument(
        "--text-color", default="#FFFFFF",
        help="Text color (default: #FFFFFF)",
    )
    parser.add_argument(
        "--stroke-color", default="#FFFFFF",
        help="Swatch stroke color (default: #FFFFFF)",
    )
    parser.add_argument(
        "--row-spacing", type=int, default=14,
        help="Vertical spacing between rows in px (default: 14)",
    )
    parser.add_argument(
        "--col-spacing", type=int, default=None,
        help="Horizontal spacing between columns in px (default: auto)",
    )
    parser.add_argument(
        "--padding", type=int, default=15,
        help="Edge padding in px (default: 15)",
    )

    return parser.parse_args()


def load_colors(colors_file):
    """Load ordered list of (feature, color) from TSV colors file.

    Returns:
        list of (feature_name, hex_color) tuples in file order
    """
    items = []
    with open(colors_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2 or parts[0].lower() == "feature":
                continue
            items.append((parts[0], parts[1]))
    return items


def strip_suffix(name):
    """Clean feature name into a display label.

    Strips _specific and _multigroup1 suffixes, replaces underscores with spaces.
    """
    for suffix in ("_specific", "_multigroup1"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.replace("_", " ")


def filter_items(items, include=None, exclude=None):
    """Filter items by include/exclude lists.

    Matching is done on the raw feature name (before suffix stripping).
    """
    if include is not None:
        include_set = {s.strip() for s in include.split(",")}
        items = [(f, c) for f, c in items if f in include_set]
    if exclude is not None:
        exclude_set = {s.strip() for s in exclude.split(",")}
        items = [(f, c) for f, c in items if f not in exclude_set]
    return items


def merge_by_color(items, merge_label_str=None):
    """Collapse features sharing the same color into one entry.

    Uses the shortest cleaned label by default, or a merge-label override.

    Args:
        items: list of (feature, color) tuples
        merge_label_str: "feature=label,feature=label,..." overrides

    Returns:
        list of (label, color) tuples (deduplicated)
    """
    # Parse merge label overrides (map raw feature name -> desired label)
    overrides = {}
    if merge_label_str:
        for pair in merge_label_str.split(","):
            if "=" in pair:
                feat, label = pair.split("=", 1)
                overrides[feat.strip()] = label.strip()

    # Group by color, preserving first-seen order
    color_groups = OrderedDict()
    for feature, color in items:
        color_upper = color.upper()
        if color_upper not in color_groups:
            color_groups[color_upper] = []
        color_groups[color_upper].append(feature)

    merged = []
    for color, features in color_groups.items():
        # Check for override label
        label = None
        for feat in features:
            if feat in overrides:
                label = overrides[feat]
                break
        if label is None:
            # Use shortest cleaned label
            cleaned = [(strip_suffix(f), f) for f in features]
            cleaned.sort(key=lambda x: len(x[0]))
            label = cleaned[0][0]
        merged.append((label, color))

    return merged


def parse_groups(groups_str):
    """Parse groups specification string.

    Format: "Header1:feat1,feat2;Header2:feat3,feat4"
    Feature names are treated as prefixes for matching.

    Returns:
        list of (header, [prefix1, prefix2, ...]) tuples
    """
    groups = []
    for group in groups_str.split(";"):
        group = group.strip()
        if ":" not in group:
            continue
        header, feats = group.split(":", 1)
        prefixes = [p.strip() for p in feats.split(",")]
        groups.append((header.strip(), prefixes))
    return groups


def group_items(items, groups_str):
    """Organize items by group with headers.

    Items matching a group's prefixes are placed under that group header.
    Unmatched items go at the end without a header.

    Args:
        items: list of (label, color) tuples
        groups_str: groups specification string

    Returns:
        list of (label, color, is_header) tuples
    """
    groups = parse_groups(groups_str)
    assigned = set()
    result = []

    for header, prefixes in groups:
        group_items_list = []
        for i, (label, color) in enumerate(items):
            if i in assigned:
                continue
            for prefix in prefixes:
                if label.lower().startswith(prefix.lower()):
                    group_items_list.append((label, color, False))
                    assigned.add(i)
                    break
        if group_items_list:
            result.append((header, None, True))
            result.extend(group_items_list)

    # Add unmatched items
    for i, (label, color) in enumerate(items):
        if i not in assigned:
            result.append((label, color, False))

    return result


def estimate_text_width(text, font_size):
    """Approximate text width in pixels.

    Uses a simple heuristic: average character width ~0.6 * font_size.
    """
    return len(text) * font_size * 0.6


def calculate_layout(items, rows=None, cols=None):
    """Calculate grid dimensions.

    Args:
        items: list of items to lay out (headers count as items in their column)
        rows: desired rows (None for auto)
        cols: desired columns (None for auto)

    Returns:
        (rows, cols) tuple
    """
    n = len(items)
    if n == 0:
        return (0, 0)

    if rows and cols:
        return (rows, cols)
    elif cols:
        rows = math.ceil(n / cols)
        return (rows, cols)
    elif rows:
        cols = math.ceil(n / rows)
        return (rows, cols)
    else:
        # Default: aim for roughly square, prefer wider
        cols = max(1, math.ceil(math.sqrt(n * 1.5)))
        rows = math.ceil(n / cols)
        return (rows, cols)


def draw_legend(items, config):
    """Generate SVG legend.

    Args:
        items: list of (label, color, is_header) tuples
        config: dict with swatch_size, font_size, background, text_color,
                stroke_color, row_spacing, col_spacing, padding, rows, cols

    Returns:
        drawsvg.Drawing
    """
    swatch = config["swatch_size"]
    font_size = config["font_size"]
    bg = config["background"]
    text_color = config["text_color"]
    stroke_color = config["stroke_color"]
    row_sp = config["row_spacing"]
    padding = config["padding"]
    rows = config["rows"]
    cols = config["cols"]

    # Separate headers from regular items, tracking which items each header covers
    has_headers = any(is_header for _, _, is_header in items)
    regular_items = []
    # Map: index in regular_items where a header's group starts -> header label
    header_at = {}
    for label, color, is_header in items:
        if is_header:
            header_at[len(regular_items)] = label
        else:
            regular_items.append((label, color, False))

    # Distribute regular items into columns (fill top-to-bottom, then next column)
    columns = [[] for _ in range(cols)]
    col_idx = 0
    row_idx = 0
    # Track which column each header lands in
    col_headers = {}  # col_index -> header label
    for i, item in enumerate(regular_items):
        if col_idx >= cols:
            break
        if i in header_at:
            # If this item starts a new group mid-column, advance to next column
            if row_idx > 0:
                col_idx += 1
                row_idx = 0
                if col_idx >= cols:
                    break
            col_headers[col_idx] = header_at[i]
        columns[col_idx].append(item)
        row_idx += 1
        if row_idx >= rows:
            row_idx = 0
            col_idx += 1

    # Calculate per-column widths based on longest label (including header)
    swatch_gap = swatch + 3  # gap between swatch and text
    col_widths = []
    for ci, col in enumerate(columns):
        if not col:
            col_widths.append(0)
            continue
        labels = [label for label, _, _ in col]
        if ci in col_headers:
            labels.append(col_headers[ci])
        max_label_w = max(estimate_text_width(l, font_size) for l in labels)
        col_widths.append(swatch_gap + max_label_w)

    col_spacing = config["col_spacing"]
    if col_spacing is None:
        col_spacing = 20  # default gap between columns

    # Canvas dimensions
    total_width = round(padding * 2 + sum(col_widths) + col_spacing * max(0, len(col_widths) - 1), 1)
    max_rows_in_col = max((len(c) for c in columns), default=0)
    header_offset = row_sp if has_headers else 0
    total_height = round(padding * 2 + max_rows_in_col * row_sp + header_offset, 1)

    d = draw.Drawing(total_width, total_height, id_prefix="legend")
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=bg))

    # Draw items
    x_offset = padding
    for ci, col in enumerate(columns):
        # Draw header if this column has one
        if ci in col_headers:
            d.append(
                draw.Text(
                    col_headers[ci],
                    font_size,
                    x_offset,
                    padding + swatch - 1,
                    fill=text_color,
                    font_family="Basic Sans",
                    font_weight="bold",
                )
            )

        # Draw regular items below the header row
        y_offset = padding + header_offset
        for label, color, is_header in col:
            # Colored swatch
            d.append(
                draw.Rectangle(
                    x_offset,
                    y_offset,
                    swatch,
                    swatch,
                    fill=color,
                    stroke=stroke_color,
                    stroke_width=0.5,
                )
            )
            # Label text
            d.append(
                draw.Text(
                    label,
                    font_size,
                    x_offset + swatch_gap,
                    y_offset + swatch - 1,
                    fill=text_color,
                    font_family="Basic Sans",
                )
            )
            y_offset += row_sp
        x_offset += col_widths[ci] + col_spacing

    return d


def main():
    args = parse_args()

    # Load colors
    print(f"Loading colors from: {args.colors}")
    items = load_colors(args.colors)
    print(f"  Loaded {len(items)} entries")

    # Filter
    items = filter_items(items, include=args.include, exclude=args.exclude)
    if args.include or args.exclude:
        print(f"  After filtering: {len(items)} entries")

    # Merge
    if args.merge_same_color:
        merge_labels = args.merge_label
        items = merge_by_color(items, merge_labels)
        print(f"  After merging same-color: {len(items)} entries")
        # Items are now (label, color) tuples — add is_header=False
        items = [(label, color, False) for label, color in items]
    else:
        # Clean labels
        items = [(strip_suffix(f), color, False) for f, color in items]

    # Group
    if args.groups:
        # Strip the is_header flag for grouping input
        flat_items = [(label, color) for label, color, _ in items]
        items = group_items(flat_items, args.groups)
        print(f"  Organized into groups: {sum(1 for _, _, h in items if h)} headers")

    # Calculate layout
    n_items = len(items)
    rows, cols = calculate_layout(items, args.rows, args.columns)
    print(f"  Layout: {rows} rows x {cols} columns ({n_items} items)")

    # Build config
    config = {
        "swatch_size": args.swatch_size,
        "font_size": args.font_size,
        "background": args.background,
        "text_color": args.text_color,
        "stroke_color": args.stroke_color,
        "row_spacing": args.row_spacing,
        "col_spacing": args.col_spacing,
        "padding": args.padding,
        "rows": rows,
        "cols": cols,
    }

    # Draw and save
    d = draw_legend(items, config)
    d.save_svg(args.output)
    print(f"Saved legend: {args.output}")
    print(f"  Dimensions: {d.width:.0f} x {d.height:.0f} px")


if __name__ == "__main__":
    main()
