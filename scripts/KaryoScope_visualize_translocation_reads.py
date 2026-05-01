#!/usr/bin/env python3
"""
KaryoScope Visualize Translocation Reads

Visualize translocation reads with multi-featureset tracks. Supports batch mode
from find_translocation_reads output TSV or direct read specification.

Usage (TSV batch mode):
    python KaryoScope_visualize_translocation_reads.py \\
        --input-tsv translocation_reads.tsv \\
        --results-dir results \\
        --output-dir outputs/translocation_reads \\
        --colors-dir /path/to/KS_human_CHM13

Usage (direct read specification):
    python KaryoScope_visualize_translocation_reads.py \\
        --reads READ_ID:SAMPLE:DATA_TYPE:REPLICATE:TRANS_TYPE \\
        --results-dir results \\
        --output-dir outputs \\
        --colors-dir /path/to/KS_human_CHM13
"""

import argparse
import gzip
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import drawsvg as draw

# Capture original command line for logging
_original_command = ' '.join(sys.argv)


class TeeLogger:
    """Write to both stdout and a log file."""
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log = open(log_path, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()


# ---------------------------------------------------------------------------
# Display constants
# ---------------------------------------------------------------------------
FS_HEIGHT = 12
ROW_SPACING = 20
MARGIN_X = 40
MARGIN_Y = 80
SAMPLE_LABEL_WIDTH = 90
TRACK_LABEL_WIDTH = 70
LABEL_WIDTH = SAMPLE_LABEL_WIDTH + TRACK_LABEL_WIDTH
BAR_GAP = 5
RIGHT_MARGIN = 40
LEGEND_WIDTH = 170
LEGEND_INTERNAL_OFFSET = 25


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# abbreviate_read_name moved to karyoplot.core.text (Phase 13.2)
from karyoplot.core.text import abbreviate_read_name  # noqa: F401, E402


def load_color_files(colors_dir, database, featuresets):
    """Load color mappings for featuresets via karyoplot's library."""
    from karyoplot.core.colors import load_featureset_palettes
    out = load_featureset_palettes(
        str(colors_dir), database, featuresets,
        on_missing="warn", value_format="tuple",
    )
    for fs in featuresets:
        print(f"  {fs}: {len(out.get(fs, {}))} colors")
    return out


def get_translocation_bed_path(results_dir, database, sample, data_type,
                               replicate, trans_type, featureset):
    """Get path to a translocation featureset BED file.

    Chromosome smoothed files use a different naming convention:
      {sample}.{db}.chromosome.smoothed.{trans_type}.translocations.bed.gz
    Other featuresets use:
      {sample}.{db}.{trans_type}.{featureset}.smoothed.translocations.bed.gz
    """
    base = (results_dir / sample / data_type / replicate
            / "KaryoScope" / database)
    if featureset == "chromosome":
        return base / (f"{sample}.{data_type}.{replicate}.{database}"
                       f".{featureset}.smoothed.{trans_type}"
                       f".translocations.bed.gz")
    return base / (f"{sample}.{data_type}.{replicate}.{database}"
                   f".{trans_type}.{featureset}.smoothed"
                   f".translocations.bed.gz")


def get_standard_bed_path(results_dir, database, sample, data_type,
                          replicate, featureset):
    """Get path to a non-translocation (standard) featureset BED file.

    Pattern: {sample}.{data_type}.{replicate}.{database}.{featureset}.smoothed.features.bed.gz
    """
    base = (results_dir / sample / data_type / replicate
            / "KaryoScope" / database)
    return base / (f"{sample}.{data_type}.{replicate}.{database}"
                   f".{featureset}.smoothed.features.bed.gz")


def load_all_features_batch(results_dir, database, sample, data_type,
                            replicate, trans_type, featuresets):
    """Batch-load ALL features from translocation BED files.

    Returns: {read_id: {featureset: [{start, stop, feature}, ...]}}
    """
    all_features = defaultdict(lambda: defaultdict(list))

    for fs in featuresets:
        bed_path = get_translocation_bed_path(
            results_dir, database, sample, data_type, replicate,
            trans_type, fs)

        if not bed_path.exists():
            print(f"    Warning: BED file not found: {bed_path}")
            continue

        with gzip.open(bed_path, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                read_id = parts[0]
                all_features[read_id][fs].append({
                    "start": int(parts[1]),
                    "stop": int(parts[2]),
                    "feature": parts[3],
                })

    return dict(all_features)


def load_features_for_reads(results_dir, database, sample, data_type,
                            replicate, read_ids, featuresets):
    """Load features for specific reads from standard (non-translocation) BED files.

    Scans each featureset BED, keeping only rows whose read_id is in read_ids.
    """
    all_features = defaultdict(lambda: defaultdict(list))

    for fs in featuresets:
        bed_path = get_standard_bed_path(
            results_dir, database, sample, data_type, replicate, fs)

        if not bed_path.exists():
            print(f"    Warning: BED file not found: {bed_path}")
            continue

        with gzip.open(bed_path, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                rid = parts[0]
                if rid not in read_ids:
                    continue
                all_features[rid][fs].append({
                    "start": int(parts[1]),
                    "stop": int(parts[2]),
                    "feature": parts[3],
                })

    return dict(all_features)


def load_translocation_reads(tsv_path):
    """Load translocation reads from TSV file."""
    reads = []
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            reads.append({
                "read_id": row["read_id"],
                "sample": row["sample"],
                "data_type": row["data_type"],
                "replicate": row["replicate"],
                "length": int(row["read_length"]),
                "type": row["translocation_type"],
            })
    return reads


def natural_sort_key(name):
    """Natural sort key for chromosome-style names (chr1, chr2, ..., chr10, chrX)."""
    c = name.lstrip('chr')
    if c.isdigit():
        return (0, int(c))
    order = {'X': 100, 'Y': 101, 'M': 102}
    return (1, order.get(c.upper(), 999))


def compute_scale_bar_bp(min_read_length):
    """Pick a round scale-bar value ~25% of the smallest read."""
    target = min_read_length * 0.25
    nice = [500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
    chosen = nice[0]
    for n in nice:
        if n <= target:
            chosen = n
        else:
            break
    return chosen


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def draw_scale_bar(d, x_start, y_pos, ratio, scale_bp, text_color="white"):
    """Draw a horizontal scale bar."""
    scale_width = scale_bp * ratio
    label = f"{scale_bp // 1000} kb" if scale_bp >= 1000 else f"{scale_bp} bp"

    d.append(draw.Line(x_start, y_pos, x_start + scale_width, y_pos,
                       stroke=text_color, stroke_width=2))
    d.append(draw.Line(x_start, y_pos - 5, x_start, y_pos + 5,
                       stroke=text_color, stroke_width=2))
    d.append(draw.Line(x_start + scale_width, y_pos - 5,
                       x_start + scale_width, y_pos + 5,
                       stroke=text_color, stroke_width=2))
    d.append(draw.Text(label, 12, x_start + scale_width / 2, y_pos - 10,
                       fill=text_color, text_anchor="middle",
                       font_family="Basic Sans"))


def draw_legend(d, featureset_colors, displayed_features, featuresets,
                legend_x, legend_y_start, text_color="white"):
    """Draw vertical color legend on the right side, grouped by featureset."""
    featureset_display_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomere',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'acrocentric': 'Acrocentric',
    }

    item_height = 14
    swatch_size = 10
    section_gap = 18
    current_y = legend_y_start
    first_section = True

    for fs in featuresets:
        fs_colors = featureset_colors.get(fs, {})
        fs_displayed = displayed_features.get(fs, set())
        if not fs_colors or not fs_displayed:
            continue

        if not first_section:
            current_y += section_gap
        first_section = False

        track_name = featureset_display_names.get(fs, fs.title())
        d.append(draw.Text(track_name, 9, legend_x, current_y,
                           fill=text_color, font_family="Basic Sans",
                           text_anchor="start", font_weight="bold"))
        current_y += 4

        # Group features by color hex
        color_to_features = {}
        for feature_name in fs_displayed:
            color_info = fs_colors.get(feature_name)
            if color_info is None:
                color_hex = '#808080'
            elif isinstance(color_info, tuple):
                color_hex = color_info[0]
            else:
                color_hex = color_info
            if color_hex not in color_to_features:
                color_to_features[color_hex] = []
            color_to_features[color_hex].append(feature_name)

        if fs == "chromosome":
            sorted_items = sorted(
                color_to_features.items(),
                key=lambda x: natural_sort_key(x[1][0]) if x[1] else (99,))
        else:
            sorted_items = sorted(
                color_to_features.items(),
                key=lambda x: x[1][0] if x[1] else '')

        for color_hex, feature_names in sorted_items:
            display_name = (feature_names[0]
                            .replace('_', ' ')
                            .replace(' specific', ''))
            d.append(draw.Rectangle(legend_x + 3, current_y,
                                    swatch_size, swatch_size, fill=color_hex))
            d.append(draw.Text(
                display_name, 9,
                legend_x + swatch_size + 6,
                current_y + swatch_size / 2 + 1,
                fill=text_color, font_family="Basic Sans",
                text_anchor="start", dominant_baseline="middle"))
            current_y += item_height


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_group(reads, features_dict, featureset_colors, featuresets,
                 output_dir, sample, data_type, replicate, trans_type,
                 title=None, part=None, total_in_group=None, pad_width=1,
                 total_parts=None, read_start=None):
    """Render a single SVG/PNG for one (sample, data_type, replicate, trans_type) group.

    part: optional part number (1-based) when a group is split across files.
    total_in_group: total reads in the full group (used in title when split).
    pad_width: zero-pad width for part numbers in filenames and titles.
    """
    n_reads = len(reads)
    if n_reads == 0:
        return

    # Sort by read length descending
    reads.sort(key=lambda x: -x["length"])

    # Pre-scan to collect displayed features for dynamic legend
    displayed_features = defaultdict(set)
    for read_info in reads:
        read_data = features_dict.get(read_info["read_id"], {})
        for fs in featuresets:
            for feat in read_data.get(fs, []):
                displayed_features[fs].add(feat["feature"])

    max_length = max(r["length"] for r in reads)
    panel_width = 1100
    ratio = (panel_width - LABEL_WIDTH - RIGHT_MARGIN) / max_length

    plot_height = (MARGIN_Y
                   + n_reads * (len(featuresets) * FS_HEIGHT + ROW_SPACING)
                   + 20)
    legend_item_count = (sum(len(v) for v in displayed_features.values())
                         + len(displayed_features) * 2)
    min_legend_height = (MARGIN_Y + legend_item_count * 14
                         + len(displayed_features) * 18 + 20)
    canvas_height = max(plot_height, min_legend_height)
    canvas_width = panel_width + 2 * MARGIN_X + LEGEND_WIDTH

    bg_color = "black"
    text_color = "white"

    d = draw.Drawing(canvas_width, canvas_height, displayInline=False)
    d.append(draw.Rectangle(0, 0, canvas_width, canvas_height, fill=bg_color))

    # Title
    if title is not None:
        display_title = title
    else:
        trans_display = trans_type.replace("_", "+")
        if part is not None:
            part_str = str(part).zfill(pad_width)
            total_str = str(total_parts).zfill(pad_width)
            read_end = read_start + n_reads - 1
            display_title = (
                f"{sample} | {data_type}.{replicate} | "
                f"{trans_display} translocation reads "
                f"(part {part_str} of {total_str}, "
                f"reads {read_start}\u2013{read_end} of {total_in_group})")
        else:
            display_title = (
                f"{sample} | {data_type}.{replicate} | "
                f"{trans_display} translocation reads (n={n_reads})")

    d.append(draw.Text(display_title, 16, canvas_width / 2, 35,
                       fill=text_color, font_weight="bold",
                       text_anchor="middle", font_family="Basic Sans"))

    bars_x_start = MARGIN_X + LABEL_WIDTH + BAR_GAP

    for j, read_info in enumerate(reads):
        read_id = read_info["read_id"]
        read_length = read_info["length"]
        read_data = features_dict.get(read_id, {})

        row_y = MARGIN_Y + j * (len(featuresets) * FS_HEIGHT + ROW_SPACING)

        # Label
        label_color = "#FF4444" if "chr2_chr13" in trans_type else "#60A5FA"
        short_id = abbreviate_read_name(read_id)

        d.append(draw.Text(short_id, 9, MARGIN_X, row_y + FS_HEIGHT - 2,
                           fill=label_color, font_family="Basic Sans",
                           font_weight="bold"))
        d.append(draw.Text(f"{read_length:,} bp", 9, MARGIN_X,
                           row_y + FS_HEIGHT + 10,
                           fill="#888888", font_family="Basic Sans"))

        # Feature tracks
        track_label_x = bars_x_start - 5
        for fs_idx, fs in enumerate(featuresets):
            fs_y = row_y + fs_idx * FS_HEIGHT

            # Track labels
            fs_label = fs.capitalize()
            if fs == "subtelomeric":
                fs_label = "Subtel"
            elif fs == "acrocentric":
                fs_label = "Acro"
            d.append(draw.Text(fs_label, 9, track_label_x,
                               fs_y + FS_HEIGHT - 2,
                               fill=text_color, font_family="Basic Sans",
                               text_anchor="end"))

            # Background
            d.append(draw.Rectangle(bars_x_start, fs_y,
                                    read_length * ratio, FS_HEIGHT,
                                    fill="#1a1a1a"))

            # Features
            for feat in read_data.get(fs, []):
                x = bars_x_start + feat["start"] * ratio
                w = max((feat["stop"] - feat["start"]) * ratio, 1)
                color, opacity = featureset_colors[fs].get(
                    feat["feature"], ("#ffffff", 1.0))
                d.append(draw.Rectangle(x, fs_y, w, FS_HEIGHT,
                                        fill=color, fill_opacity=opacity))

            # Border
            d.append(draw.Rectangle(bars_x_start, fs_y,
                                    read_length * ratio, FS_HEIGHT,
                                    fill="none", stroke="#333333",
                                    stroke_width=0.5))

        # Outer border
        total_height = len(featuresets) * FS_HEIGHT
        d.append(draw.Rectangle(bars_x_start, row_y,
                                read_length * ratio, total_height,
                                fill="none", stroke="#555555",
                                stroke_width=1))

    # Scale bar (top, below title)
    scale_y = 68
    min_length = min(r["length"] for r in reads)
    scale_bp = compute_scale_bar_bp(min_length)
    draw_scale_bar(d, bars_x_start, scale_y, ratio, scale_bp, text_color)

    # Vertical legend (right side)
    legend_x = panel_width + 2 * MARGIN_X + LEGEND_INTERNAL_OFFSET
    draw_legend(d, featureset_colors, displayed_features, featuresets,
                legend_x, MARGIN_Y, text_color)

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    if part is not None:
        part_str = str(part).zfill(pad_width)
        total_str = str(total_parts).zfill(pad_width)
        stem = (f"{sample}.{data_type}.{replicate}.{trans_type}"
                f".n{total_in_group}"
                f".part{part_str}of{total_str}.translocation_reads")
    else:
        stem = (f"{sample}.{data_type}.{replicate}.{trans_type}"
                f".n{n_reads}.translocation_reads")
    svg_path = output_dir / f"{stem}.svg"
    png_path = output_dir / f"{stem}.png"

    d.save_svg(str(svg_path))

    try:
        zoom = min(2.0 + max(0, (n_reads - 5)) // 5 * 0.5, 4.0)
        subprocess.run(["rsvg-convert", f"--zoom={zoom}",
                        "-o", str(png_path), str(svg_path)],
                       check=True, capture_output=True)
        print(f"    Saved: {png_path.name} ({zoom}x zoom)")
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    finally:
        svg_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_read_spec(spec):
    """Parse a single read specification string.

    Format: READ_ID:SAMPLE:DATA_TYPE:REPLICATE:TRANS_TYPE
    """
    parts = spec.split(":")
    if len(parts) != 5:
        print(f"Error: invalid read spec '{spec}'. "
              f"Expected READ_ID:SAMPLE:DATA_TYPE:REPLICATE:TRANS_TYPE",
              file=sys.stderr)
        sys.exit(1)
    return {
        "read_id": parts[0],
        "sample": parts[1],
        "data_type": parts[2],
        "replicate": parts[3],
        "type": parts[4],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Visualize translocation reads with multi-featureset "
                    "tracks. Supports batch mode from find_translocation_reads "
                    "output TSV or direct read specification.",
        formatter_class=argparse.RawTextHelpFormatter)

    # Input mode (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input-tsv", type=Path,
        help="Batch mode: path to find_translocation_reads output TSV.")
    input_group.add_argument(
        "--reads", nargs="+", metavar="SPEC",
        help="Direct read specification.\n"
             "Each spec: READ_ID:SAMPLE:DATA_TYPE:REPLICATE:TRANS_TYPE")

    # Required directories
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Root results directory")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for generated images")
    parser.add_argument("--colors-dir", type=Path, required=True,
                        help="Directory containing {database}.{featureset}.colors.txt files")

    # Optional
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--title",
                        help="Custom title for the visualization (auto-generated if unset)")
    parser.add_argument("--featuresets", nargs="+",
                        default=["chromosome", "subtelomeric", "acrocentric", "repeat", "region"],
                        help="Featuresets to display (default: chromosome subtelomeric acrocentric repeat region)")
    parser.add_argument("--max-reads-per-file", type=int, default=50, dest="max_reads_per_file",
                        help="Maximum reads per output file before splitting (default: 50)")
    parser.add_argument("--filter",
                        help="Only render groups whose key contains this substring (TSV mode)")
    parser.add_argument("--parts", type=int,
                        help="Only render the first N parts of each group (TSV mode)")
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to log file (default: True)")

    args = parser.parse_args()

    # Validate directories
    for name, path in [("--results-dir", args.results_dir),
                       ("--colors-dir", args.colors_dir)]:
        if not path.is_dir():
            print(f"Error: {name} does not exist: {path}", file=sys.stderr)
            sys.exit(1)

    # Set up logging
    if args.log_file:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        log_path = args.output_dir / "visualize_translocation_reads.log"
        sys.stdout = TeeLogger(str(log_path))

    print("=" * 60)
    print("KaryoScope Visualize Translocation Reads")
    print("=" * 60)

    # Parameters table
    mode = "TSV batch" if args.input_tsv else "direct reads"
    print(f"\n{'Parameter':<22} {'Value':<35}")
    print(f"{'-' * 22} {'-' * 35}")
    print(f"{'mode':<22} {mode}")
    if args.input_tsv:
        print(f"{'input-tsv':<22} {args.input_tsv}")
    else:
        print(f"{'reads':<22} {len(args.reads)} read spec(s)")
    print(f"{'results-dir':<22} {args.results_dir}")
    print(f"{'output-dir':<22} {args.output_dir}")
    print(f"{'colors-dir':<22} {args.colors_dir}")
    print(f"{'database':<22} {args.database}")
    print(f"{'featuresets':<22} {' '.join(args.featuresets)}")
    print(f"{'max-reads-per-file':<22} {args.max_reads_per_file}")
    if args.title:
        print(f"{'title':<22} {args.title}")
    if args.filter:
        print(f"{'filter':<22} {args.filter}")
    if args.parts:
        print(f"{'parts':<22} {args.parts}")
    print(f"{'log-file':<22} {args.log_file}")

    print(f"\n{'=' * 60}")
    print("Command")
    print(f"{'=' * 60}")
    print(_original_command)

    # Load colors
    print("\nLoading color files...")
    featureset_colors = load_color_files(
        args.colors_dir, args.database, args.featuresets)

    if args.input_tsv:
        _run_tsv_mode(args, featureset_colors)
    else:
        _run_reads_mode(args, featureset_colors)

    print("\nDone.")


def _run_tsv_mode(args, featureset_colors):
    """Batch mode: read TSV, group, split, and render."""
    print(f"\nLoading reads from: {args.input_tsv}")
    all_reads = load_translocation_reads(args.input_tsv)
    print(f"  Found {len(all_reads)} total reads\n")

    # Group reads by (sample, data_type, replicate, trans_type)
    groups = defaultdict(list)
    for r in all_reads:
        key = (r["sample"], r["data_type"], r["replicate"], r["type"])
        groups[key].append(r)

    print(f"{len(groups)} groups to render\n")

    for (sample, data_type, replicate, trans_type) in sorted(groups):
        group_key = f"{sample}.{data_type}.{replicate}.{trans_type}"
        if args.filter and args.filter not in group_key:
            continue

        group_reads = groups[(sample, data_type, replicate, trans_type)]
        n = len(group_reads)
        print(f"  {sample}.{data_type}.{replicate} / {trans_type} "
              f"({n} reads)...")

        # Batch load features
        features_dict = load_all_features_batch(
            args.results_dir, args.database,
            sample, data_type, replicate, trans_type, args.featuresets)

        # Sort once before chunking so parts are consistently ordered
        group_reads.sort(key=lambda x: -x["length"])

        if n > args.max_reads_per_file:
            n_parts = ((n + args.max_reads_per_file - 1)
                       // args.max_reads_per_file)
            pad_width = len(str(n_parts))
            max_parts = args.parts if args.parts else n_parts
            print(f"    Splitting into {n_parts} parts "
                  f"({args.max_reads_per_file} reads each), "
                  f"rendering {max_parts}")
            for i in range(min(max_parts, n_parts)):
                start_idx = i * args.max_reads_per_file
                end_idx = (i + 1) * args.max_reads_per_file
                chunk = group_reads[start_idx:end_idx]
                render_group(
                    chunk, features_dict, featureset_colors,
                    args.featuresets, args.output_dir,
                    sample, data_type, replicate, trans_type,
                    title=args.title,
                    part=i + 1, total_in_group=n, pad_width=pad_width,
                    total_parts=n_parts,
                    read_start=start_idx + 1)
        else:
            render_group(
                group_reads, features_dict, featureset_colors,
                args.featuresets, args.output_dir,
                sample, data_type, replicate, trans_type,
                title=args.title)


def _run_reads_mode(args, featureset_colors):
    """Direct read specification mode."""
    read_specs = [_parse_read_spec(s) for s in args.reads]

    print(f"\nProcessing {len(read_specs)} directly specified read(s)...\n")

    # Group by (sample, data_type, replicate) to batch-load standard BED files
    load_groups = defaultdict(list)
    for spec in read_specs:
        key = (spec["sample"], spec["data_type"], spec["replicate"])
        load_groups[key].append(spec)

    # Load features and compute read lengths
    features_dict = {}
    for (sample, data_type, replicate), specs in load_groups.items():
        read_ids = {s["read_id"] for s in specs}
        print(f"  Loading features for {len(read_ids)} read(s) from "
              f"{sample}.{data_type}.{replicate}...")
        group_features = load_features_for_reads(
            args.results_dir, args.database,
            sample, data_type, replicate, read_ids, args.featuresets)
        features_dict.update(group_features)

    # Compute read lengths from loaded features (max stop across all featuresets)
    for spec in read_specs:
        rid = spec["read_id"]
        rd = features_dict.get(rid, {})
        max_stop = 0
        for fs_features in rd.values():
            for feat in fs_features:
                if feat["stop"] > max_stop:
                    max_stop = feat["stop"]
        spec["length"] = max_stop
        if spec["length"] == 0:
            print(f"    Warning: read {rid} has no features / zero length")

    # Filter out zero-length reads
    read_specs = [s for s in read_specs if s["length"] > 0]
    if not read_specs:
        print("Error: no reads with features found.", file=sys.stderr)
        sys.exit(1)

    # Use first read's group info for the filename; all go into one file
    first = read_specs[0]
    render_group(
        read_specs, features_dict, featureset_colors,
        args.featuresets, args.output_dir,
        first["sample"], first["data_type"], first["replicate"], first["type"],
        title=args.title)


if __name__ == "__main__":
    main()
