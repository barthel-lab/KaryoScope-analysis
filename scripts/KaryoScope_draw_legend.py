#!/usr/bin/env python3
"""KaryoScope_draw_legend.py — generate SVG legends from KaryoScope color files.

Thin CLI wrapper around :func:`karyoplot.svg.legend.make_legend_drawing`.

Usage:
    # Auto-layout in 4 columns (defaults to dark theme)
    python KaryoScope_draw_legend.py \\
        --colors KS_human_CHM13.chromosome_acrocentric.colors.txt \\
        --output legend.svg --columns 4

    # Merge features sharing a color, with custom labels
    python KaryoScope_draw_legend.py \\
        --colors ... --output legend.svg --columns 4 \\
        --merge-same-color \\
        --merge-label "chr13_specific=acrocentric,autosome_multigroup1=categorized"

    # Category headers + filtering
    python KaryoScope_draw_legend.py \\
        --colors ... --output legend.svg --columns 3 \\
        --exclude "categorized" \\
        --groups "Chromosomes:chr1,chr2,chr3;Features:DJ,PJ,rDNA"
"""

import argparse

from karyoplot.core.theme import Theme
from karyoplot.svg.legend import (
    make_legend_drawing,
    merge_by_color,
    strip_label_suffixes,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate SVG legends from KaryoScope color mapping files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--colors", required=True,
                        help="Path to TSV colors file (feature, color columns)")
    parser.add_argument("--output", default="legend.svg",
                        help="Output SVG path (default: legend.svg)")

    # Layout
    parser.add_argument("--columns", type=int, default=None,
                        help="Number of columns (auto-calculated if not set)")
    parser.add_argument("--rows", type=int, default=None,
                        help="Number of rows (auto-calculated if not set)")

    # Filtering
    parser.add_argument("--include", default=None,
                        help="Comma-separated feature names to include")
    parser.add_argument("--exclude", default=None,
                        help="Comma-separated feature names to exclude")

    # Merging
    parser.add_argument("--merge-same-color", action="store_true",
                        help="Collapse features sharing the same hex color")
    parser.add_argument("--merge-label", default=None,
                        help="Override merged labels: feature=label,feature=label,...")

    # Grouping
    parser.add_argument("--groups", default=None,
                        help="Category headers: 'Header1:feat1,feat2;Header2:feat3,feat4'\n"
                             "Features matched by prefix")

    # Styling
    parser.add_argument("--swatch-size", type=int, default=8)
    parser.add_argument("--font-size", type=int, default=12)
    parser.add_argument("--background", default="#000000")
    parser.add_argument("--text-color", default="#FFFFFF")
    parser.add_argument("--stroke-color", default="#FFFFFF")
    parser.add_argument("--row-spacing", type=int, default=14)
    parser.add_argument("--col-spacing", type=int, default=None)
    parser.add_argument("--padding", type=int, default=15)

    return parser.parse_args()


def load_colors_file(filepath: str) -> list[tuple[str, str]]:
    """Load ordered list of ``(feature, color)`` from TSV colors file."""
    items = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2 or parts[0].lower() == "feature":
                continue
            items.append((parts[0], parts[1]))
    return items


def filter_items(items, include=None, exclude=None):
    if include is not None:
        include_set = {s.strip() for s in include.split(",")}
        items = [(f, c) for f, c in items if f in include_set]
    if exclude is not None:
        exclude_set = {s.strip() for s in exclude.split(",")}
        items = [(f, c) for f, c in items if f not in exclude_set]
    return items


def parse_groups(groups_str: str):
    """Parse 'Header1:feat1,feat2;Header2:feat3,feat4' into [(header, [prefixes])]."""
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
    """Organize items by group, inserting headers as ``(label, color, True)`` rows."""
    groups = parse_groups(groups_str)
    assigned = set()
    result = []
    for header, prefixes in groups:
        bucket = []
        for i, (label, color) in enumerate(items):
            if i in assigned:
                continue
            if any(label.lower().startswith(p.lower()) for p in prefixes):
                bucket.append((label, color, False))
                assigned.add(i)
        if bucket:
            result.append((header, "", True))
            result.extend(bucket)
    for i, (label, color) in enumerate(items):
        if i not in assigned:
            result.append((label, color, False))
    return result


def main():
    args = parse_args()

    print(f"Loading colors from: {args.colors}")
    items = load_colors_file(args.colors)
    print(f"  Loaded {len(items)} entries")

    items = filter_items(items, include=args.include, exclude=args.exclude)
    if args.include or args.exclude:
        print(f"  After filtering: {len(items)} entries")

    if args.merge_same_color:
        overrides = {}
        if args.merge_label:
            for pair in args.merge_label.split(","):
                if "=" in pair:
                    feat, label = pair.split("=", 1)
                    overrides[feat.strip()] = label.strip()
        merged = merge_by_color(items, label_overrides=overrides)
        items = [(label, color, False) for label, color in merged]
        print(f"  After merging same-color: {len(items)} entries")
    else:
        items = [(strip_label_suffixes(f), color, False) for f, color in items]

    if args.groups:
        flat = [(label, color) for label, color, _ in items]
        items = group_items(flat, args.groups)
        print(f"  Organized into groups: {sum(1 for _, _, h in items if h)} headers")

    # Build a Theme that matches the legacy CLI args (background/text/font)
    theme = Theme(
        name="custom",
        background=args.background,
        text=args.text_color,
        line=args.text_color,
        muted_line=args.text_color,
    )

    d = make_legend_drawing(
        items,
        theme=theme,
        rows=args.rows,
        cols=args.columns,
        swatch_size=args.swatch_size,
        font_size=args.font_size,
        row_spacing=args.row_spacing,
        col_spacing=args.col_spacing,
        padding=args.padding,
        stroke_color=args.stroke_color,
    )
    d.save_svg(args.output)
    print(f"Saved legend: {args.output}")
    print(f"  Dimensions: {d.width:.0f} x {d.height:.0f} px")


if __name__ == "__main__":
    main()
