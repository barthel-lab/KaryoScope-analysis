#!/usr/bin/env python3
"""
KaryoScope_telogator_reads_viz.py

Visualize telomeric reads as vertical bars with region (satellite) features.
Reads are displayed side by side, sorted by sample then by read length.

Usage:
    python KaryoScope_telogator_reads_viz.py \
        --samples NHA_p1 NHA_E6E7_3_PDL48 ... \
        --results-dir results \
        --colors KS_human_CHM13.region.colors.txt \
        --output output.svg
"""

import argparse
import gzip
import os
import subprocess
from collections import defaultdict

import drawsvg as draw
from PIL import Image, ImageDraw, ImageFont


def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def draw_reads_png(reads, colors, output_path, config):
    """Draw reads directly to PNG using Pillow (handles large images better than SVG)."""
    bar_width = config["bar_width"]
    ratio = config["ratio"]
    top_margin = config["top_margin"]
    left_margin = config["left_margin"]
    bottom_margin = config["bottom_margin"]
    read_spacing = config["read_spacing"]
    sample_spacing = config["sample_spacing"]
    background = config["background"]
    sample_order = config["sample_order"]

    # Calculate dimensions
    max_length = max(r[2] for r in reads)
    max_height_px = int(max_length * ratio)

    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1

    total_reads = len(reads)
    num_samples = len(sample_read_counts)

    effective_left_margin = 20  # Reduced since scale bar is separate
    image_width = (
        effective_left_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + 50
    )
    image_height = top_margin + max_height_px + bottom_margin

    # Create image
    bg_color = (0, 0, 0) if background == "black" else (255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    img = Image.new('RGB', (image_width, image_height), bg_color)
    draw_ctx = ImageDraw.Draw(img)

    # Try to load font
    try:
        font = ImageFont.truetype("Basic Sans", 16)
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("/Users/fbarthel/Documents/Barthel-Custom-Powerpoint-Theme/fonts/BasicSans-Regular.otf", 16)
        except (OSError, IOError):
            font = ImageFont.load_default()

    # Track positions
    current_x = effective_left_margin
    current_sample = None
    sample_x_start = {}
    sample_x_end = {}

    # Draw reads
    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_x_end[current_sample] = current_x
            current_x += sample_spacing

        if sample not in sample_x_start:
            sample_x_start[sample] = current_x

        current_sample = sample

        for start, end, feature in features:
            y_start = top_margin + int(start * ratio)
            height = max(1, int((end - start) * ratio))
            color_hex = colors.get(feature, "#ffffff")
            color_rgb = hex_to_rgb(color_hex)

            draw_ctx.rectangle(
                [current_x, y_start, current_x + bar_width, y_start + height],
                fill=color_rgb
            )

        current_x += bar_width + read_spacing

    if current_sample:
        sample_x_end[current_sample] = current_x

    # Draw labels and lines
    label_interval = 300
    read_width = bar_width + read_spacing

    for sample in sample_order:
        if sample in sample_x_start and sample in sample_x_end:
            # White line
            draw_ctx.line(
                [(sample_x_start[sample], top_margin - 5),
                 (sample_x_end[sample], top_margin - 5)],
                fill=text_color, width=2
            )

            # Labels
            label_text = sample.replace("_", " ")
            num_reads = sample_read_counts[sample]

            for i in range(0, num_reads, label_interval):
                x_pos = sample_x_start[sample] + (i * read_width)
                draw_ctx.text((x_pos, top_margin - 25), label_text, fill=text_color, font=font)

    # Draw separator lines
    for sample in sample_order[:-1]:
        if sample in sample_x_end:
            sep_x = sample_x_end[sample] + sample_spacing // 2
            draw_ctx.line(
                [(sep_x, top_margin), (sep_x, top_margin + max_height_px)],
                fill=(*text_color, 77), width=1  # 0.3 opacity approximation
            )

    img.save(output_path)
    print(f"Saved PNG: {output_path}")
    print(f"  Dimensions: {image_width} x {image_height} pixels")
    return image_height


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize telomeric reads as vertical bars with region features"
    )

    # Input options (either --samples + --results-dir OR --bed)
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Sample names to include (in display order). Use with --results-dir.",
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Base results directory containing sample subdirectories. Use with --samples.",
    )
    parser.add_argument(
        "--bed",
        nargs="+",
        default=None,
        help="BED files to visualize directly (alternative to --samples). "
             "Format: SAMPLE:PATH or just PATH (sample name derived from filename).",
    )
    parser.add_argument(
        "--colors",
        required=True,
        help="Path to region colors file (e.g., KS_human_CHM13.region.colors.txt)",
    )
    parser.add_argument(
        "--output", "-o", required=True, help="Output SVG file path"
    )
    parser.add_argument(
        "--scale-bar-output",
        default=None,
        help="Output SVG file path for separate scale bar (optional)",
    )

    # Data options
    parser.add_argument(
        "--database",
        default="KS_human_CHM13",
        help="Database name (default: KS_human_CHM13)",
    )
    parser.add_argument(
        "--featureset",
        default="region",
        help="Featureset to use (default: region)",
    )
    parser.add_argument(
        "--smoothness",
        default="smoothed",
        help="Feature smoothness level (default: smoothed)",
    )

    # Display options
    parser.add_argument(
        "--bar-width",
        type=int,
        default=3,
        help="Width of each read bar in pixels (default: 3)",
    )
    parser.add_argument(
        "--read-spacing",
        type=int,
        default=3,
        help="Horizontal spacing between reads (default: 3)",
    )
    parser.add_argument(
        "--sample-spacing",
        type=int,
        default=20,
        help="Horizontal spacing between sample groups (default: 20)",
    )
    parser.add_argument(
        "--ratio",
        type=float,
        default=1 / 300,
        help="bp to pixel ratio for read height (default: 1/300)",
    )
    parser.add_argument(
        "--background",
        default="black",
        choices=["white", "black"],
        help="Background color (default: black)",
    )
    parser.add_argument(
        "--top-margin",
        type=int,
        default=80,
        help="Top margin for labels (default: 80)",
    )
    parser.add_argument(
        "--left-margin",
        type=int,
        default=60,
        help="Left margin for scale bar (default: 60)",
    )
    parser.add_argument(
        "--bottom-margin",
        type=int,
        default=40,
        help="Bottom margin (default: 40)",
    )
    parser.add_argument(
        "--orient-telomere-top",
        action="store_true",
        help="Reorient reads so telomere features are always at the top",
    )

    return parser.parse_args()


def load_color_mapping(colors_file):
    """Load feature -> color mapping from colors file.

    Handles both '_specific' suffix variants and bare feature names.

    Returns:
        dict: feature_name -> hex_color
    """
    colors = {"novel": "#ffffff"}  # Default for unknown features

    with open(colors_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2 or parts[0].lower() == "feature":
                continue
            feature, color = parts[0], parts[1]
            colors[feature] = color
            # Also map without _specific suffix for smoothed BED files
            if feature.endswith("_specific"):
                colors[feature[:-9]] = color
            # Also map with _specific suffix
            if not feature.endswith("_specific") and not feature.endswith("_multigroup1"):
                colors[feature + "_specific"] = color

    return colors


def load_sample_bed_data(samples, results_dir, database, featureset, smoothness):
    """Load BED data for all samples.

    Returns:
        list of tuples: [(sample, read_id, read_length, features), ...]
        where features = [(start, end, feature_name), ...]
    """
    all_reads = []

    for sample in samples:
        # Build path to BED file
        bed_path = os.path.join(
            results_dir,
            sample,
            "telogator/1/KaryoScope",
            database,
            f"{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed",
        )

        # Try gzipped version if uncompressed doesn't exist
        if not os.path.exists(bed_path):
            bed_path_gz = bed_path + ".gz"
            if os.path.exists(bed_path_gz):
                bed_path = bed_path_gz
            else:
                print(f"Warning: BED file not found for {sample}: {bed_path}")
                continue

        # Parse BED file, group by read
        read_features = defaultdict(list)
        open_func = gzip.open if bed_path.endswith(".gz") else open
        mode = "rt" if bed_path.endswith(".gz") else "r"

        with open_func(bed_path, mode) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    read_id = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    feature = parts[3]
                    read_features[read_id].append((start, end, feature))

        # Calculate read lengths and add to list
        for read_id, features in read_features.items():
            read_length = max(end for _, end, _ in features)
            all_reads.append((sample, read_id, read_length, features))

        print(f"  {sample}: {len(read_features)} reads loaded")

    return all_reads


def load_bed_files_direct(bed_specs):
    """Load BED data directly from file paths.

    Args:
        bed_specs: List of "SAMPLE:PATH" or "PATH" strings

    Returns:
        tuple: (all_reads, sample_order)
        where all_reads = [(sample, read_id, read_length, features), ...]
    """
    all_reads = []
    sample_order = []

    for spec in bed_specs:
        # Parse SAMPLE:PATH or just PATH
        if ":" in spec and not spec.startswith("/"):
            sample, bed_path = spec.split(":", 1)
        else:
            bed_path = spec
            # Derive sample name from filename
            basename = os.path.basename(bed_path)
            sample = basename.split(".")[0]

        sample_order.append(sample)

        if not os.path.exists(bed_path):
            print(f"Warning: BED file not found: {bed_path}")
            continue

        # Parse BED file, group by read
        read_features = defaultdict(list)
        open_func = gzip.open if bed_path.endswith(".gz") else open
        mode = "rt" if bed_path.endswith(".gz") else "r"

        with open_func(bed_path, mode) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    read_id = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    feature = parts[3]
                    read_features[read_id].append((start, end, feature))

        # Calculate read lengths and add to list
        for read_id, features in read_features.items():
            read_length = max(end for _, end, _ in features)
            all_reads.append((sample, read_id, read_length, features))

        print(f"  {sample}: {len(read_features)} reads loaded")

    return all_reads, sample_order


def sort_reads(reads, sample_order):
    """Sort reads by sample order, then by length (descending) within sample.

    Args:
        reads: List of (sample, read_id, read_length, features) tuples
        sample_order: List of sample names in desired order

    Returns:
        Sorted list of reads
    """
    sample_rank = {s: i for i, s in enumerate(sample_order)}

    def sort_key(read):
        sample, read_id, length, features = read
        return (sample_rank.get(sample, 999), -length)  # Descending length

    return sorted(reads, key=sort_key)


# Telomere features used for orientation detection
TELOMERE_FEATURES = {'canonical_telomere', 'noncanonical_telomere'}


def orient_telomere_top(reads):
    """Reorient reads so telomere features are at the top (position 0).

    For each read, checks if telomere features are closer to start or end.
    If closer to end, flips the read coordinates.

    Args:
        reads: List of (sample, read_id, read_length, features) tuples

    Returns:
        List of reoriented reads
    """
    oriented_reads = []
    flipped_count = 0

    for sample, read_id, read_length, features in reads:
        # Find telomere positions
        telomere_positions = []
        for start, end, feature in features:
            if feature in TELOMERE_FEATURES:
                telomere_positions.extend([start, end])

        if telomere_positions:
            # Calculate average telomere position
            avg_telomere_pos = sum(telomere_positions) / len(telomere_positions)
            midpoint = read_length / 2

            # If telomere is in second half, flip the read
            if avg_telomere_pos > midpoint:
                # Flip coordinates: new_start = length - old_end, new_end = length - old_start
                flipped_features = [
                    (read_length - end, read_length - start, feature)
                    for start, end, feature in features
                ]
                # Sort by start position
                flipped_features.sort(key=lambda x: x[0])
                oriented_reads.append((sample, read_id, read_length, flipped_features))
                flipped_count += 1
            else:
                oriented_reads.append((sample, read_id, read_length, features))
        else:
            # No telomere features, keep as-is
            oriented_reads.append((sample, read_id, read_length, features))

    print(f"  Reoriented {flipped_count} of {len(reads)} reads (telomere now at top)")
    return oriented_reads


def draw_scale_bar(d, x, y, ratio, text_color, max_height_px):
    """Draw a vertical scale bar showing read length scale."""
    scale_bar_bp = 10000  # 10 kbp
    scale_bar_height = int(scale_bar_bp * ratio)

    # Don't draw if scale bar is too tall
    if scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    # Draw vertical bar
    d.append(draw.Rectangle(x, y, 3, scale_bar_height, fill=text_color))

    # Draw tick marks at top and bottom
    d.append(draw.Line(x - 3, y, x + 6, y, stroke=text_color, stroke_width=1))
    d.append(
        draw.Line(
            x - 3,
            y + scale_bar_height,
            x + 6,
            y + scale_bar_height,
            stroke=text_color,
            stroke_width=1,
        )
    )

    # Draw label (rotated)
    label = f"{scale_bar_bp // 1000} Kbp"
    label_x = x - 10
    label_y = y + scale_bar_height / 2
    d.append(
        draw.Text(
            label,
            12,
            label_x,
            label_y,
            fill=text_color,
            text_anchor="middle",
            transform=f"rotate(-90, {label_x}, {label_y})",
        )
    )


def draw_scale_bar_svg(output_path, image_height, top_margin, ratio, background):
    """Draw a separate SVG file containing just the scale bar."""
    text_color = "#ffffff" if background == "black" else "#000000"
    max_height_px = image_height - top_margin - 40  # Approximate bottom margin

    # Calculate scale bar dimensions
    scale_bar_bp = 10000  # 10 kbp
    scale_bar_height = int(scale_bar_bp * ratio)
    if scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    # Create a narrow SVG for the scale bar
    width = 60
    d = draw.Drawing(width, image_height, id_prefix="sb")
    d.append(draw.Rectangle(0, 0, width, image_height, fill=background))

    x = 45
    y = top_margin

    # Draw vertical bar
    d.append(draw.Rectangle(x, y, 3, scale_bar_height, fill=text_color))

    # Draw tick marks at top and bottom
    d.append(draw.Line(x - 3, y, x + 6, y, stroke=text_color, stroke_width=1))
    d.append(
        draw.Line(
            x - 3,
            y + scale_bar_height,
            x + 6,
            y + scale_bar_height,
            stroke=text_color,
            stroke_width=1,
        )
    )

    # Draw label (rotated)
    label = f"{scale_bar_bp // 1000} Kbp"
    label_x = x - 15
    label_y = y + scale_bar_height / 2
    d.append(
        draw.Text(
            label,
            16,
            label_x,
            label_y,
            fill=text_color,
            text_anchor="middle",
            font_family="Basic Sans",
            transform=f"rotate(-90, {label_x}, {label_y})",
        )
    )

    d.save_svg(output_path)
    print(f"Saved scale bar: {output_path}")
    print(f"  Dimensions: {width} x {image_height} pixels")

    # Render scale bar PNG (small enough for rsvg-convert)
    png_path = output_path.rsplit('.', 1)[0] + '.png'
    try:
        subprocess.run(
            ['rsvg-convert', '-o', png_path, output_path],
            check=True,
            capture_output=True,
        )
        print(f"Saved scale bar PNG: {png_path}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass  # Scale bar PNG optional


def draw_reads_svg(reads, colors, output_path, config):
    """Draw all reads as vertical bars in an SVG.

    Args:
        reads: Sorted list of (sample, read_id, read_length, features)
        colors: Feature color mapping
        output_path: Output SVG file path
        config: Dict with bar_width, ratio, background, margins, etc.

    Returns:
        image_height: The height of the generated image (for scale bar alignment)
    """
    bar_width = config["bar_width"]
    ratio = config["ratio"]
    top_margin = config["top_margin"]
    left_margin = config["left_margin"]
    bottom_margin = config["bottom_margin"]
    read_spacing = config["read_spacing"]
    sample_spacing = config["sample_spacing"]
    background = config["background"]
    sample_order = config["sample_order"]
    draw_scale = config.get("draw_scale_bar", True)

    # Calculate max read length for scaling
    max_length = max(r[2] for r in reads)
    max_height_px = int(max_length * ratio)

    # Count reads per sample for width calculation
    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1

    total_reads = len(reads)
    num_samples = len(sample_read_counts)

    # Calculate image dimensions
    # Reduce left margin if scale bar is separate
    effective_left_margin = left_margin if draw_scale else 20
    image_width = (
        effective_left_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + 50  # Right margin
    )
    image_height = top_margin + max_height_px + bottom_margin

    # Create drawing
    d = draw.Drawing(image_width, image_height, id_prefix="tr")
    text_color = "#ffffff" if background == "black" else "#000000"
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background))

    # Draw scale bar (left side) only if not creating separate file
    if draw_scale:
        draw_scale_bar(d, left_margin - 40, top_margin, ratio, text_color, max_height_px)

    # Use effective left margin for positioning reads
    left_margin = effective_left_margin

    # Track x position and sample boundaries
    current_x = left_margin
    current_sample = None
    sample_x_start = {}
    sample_x_end = {}

    for sample, read_id, read_length, features in reads:
        # Add sample spacing when sample changes
        if current_sample is not None and sample != current_sample:
            sample_x_end[current_sample] = current_x
            current_x += sample_spacing

        if sample not in sample_x_start:
            sample_x_start[sample] = current_x

        current_sample = sample

        # Draw this read's features as colored rectangles
        for start, end, feature in features:
            y_start = top_margin + int(start * ratio)
            height = max(1, int((end - start) * ratio))
            color = colors.get(feature, "#ffffff")

            d.append(
                draw.Rectangle(current_x, y_start, bar_width, height, fill=color)
            )

        current_x += bar_width + read_spacing

    # Record final sample boundary
    if current_sample:
        sample_x_end[current_sample] = current_x

    # Draw horizontal line and sample labels (repeated every ~300 reads)
    label_interval = 300  # Repeat label every N reads
    font_family = "Basic Sans"

    for sample in sample_order:
        if sample in sample_x_start and sample in sample_x_end:
            # Draw horizontal line at top spanning this sample
            d.append(
                draw.Line(
                    sample_x_start[sample],
                    top_margin - 5,
                    sample_x_end[sample],
                    top_margin - 5,
                    stroke=text_color,
                    stroke_width=2,
                )
            )

            # Parse label: replace underscores with spaces
            label_text = sample.replace("_", " ")

            # Calculate label positions (starting at first read, then every N reads)
            num_reads_in_sample = sample_read_counts[sample]
            read_width = bar_width + read_spacing

            # Determine label positions - first label at start, then every interval
            label_positions = []
            for i in range(0, num_reads_in_sample, label_interval):
                x_pos = sample_x_start[sample] + (i * read_width)
                label_positions.append(x_pos)

            # Draw labels above the white line
            for x_pos in label_positions:
                d.append(
                    draw.Text(
                        label_text,
                        16,
                        x_pos,
                        top_margin - 12,
                        fill=text_color,
                        text_anchor="start",
                        font_family=font_family,
                    )
                )

    # Draw separator lines between samples
    for sample in sample_order[:-1]:
        if sample in sample_x_end:
            sep_x = sample_x_end[sample] + sample_spacing / 2
            d.append(
                draw.Line(
                    sep_x,
                    top_margin,
                    sep_x,
                    top_margin + max_height_px,
                    stroke=text_color,
                    stroke_width=0.5,
                    stroke_opacity=0.3,
                )
            )

    # Save SVG
    d.save_svg(output_path)
    print(f"\nSaved: {output_path}")
    print(f"  Dimensions: {image_width} x {image_height} pixels")

    return image_height


def main():
    args = parse_args()

    print("NHA Telomeric Reads Region Feature Visualization")
    print("=" * 50)

    # Validate input arguments
    if args.bed:
        # Direct BED file mode
        use_direct_bed = True
    elif args.samples and args.results_dir:
        # Sample + results-dir mode
        use_direct_bed = False
    else:
        print("Error: Must provide either --bed OR (--samples AND --results-dir)")
        return

    # Load color mapping
    print(f"\nLoading colors from: {args.colors}")
    colors = load_color_mapping(args.colors)
    print(f"  Loaded {len(colors)} color mappings")

    # Load BED data
    if use_direct_bed:
        print(f"\nLoading feature data from BED files...")
        reads, sample_order = load_bed_files_direct(args.bed)
    else:
        print(f"\nLoading feature data from: {args.results_dir}")
        reads = load_sample_bed_data(
            args.samples,
            args.results_dir,
            args.database,
            args.featureset,
            args.smoothness,
        )
        sample_order = args.samples

    print(f"  Total reads: {len(reads)}")

    if not reads:
        print("Error: No reads loaded. Check sample names and paths.")
        return

    # Orient reads so telomere is at top (if requested)
    if args.orient_telomere_top:
        print(f"\nOrienting reads (telomere at top)...")
        reads = orient_telomere_top(reads)

    # Sort reads
    print(f"\nSorting reads by sample order, then by length (descending)")
    sorted_reads = sort_reads(reads, sample_order)

    # Generate SVG
    print(f"\nGenerating SVG: {args.output}")
    config = {
        "bar_width": args.bar_width,
        "ratio": args.ratio,
        "background": args.background,
        "top_margin": args.top_margin,
        "left_margin": args.left_margin,
        "bottom_margin": args.bottom_margin,
        "read_spacing": args.read_spacing,
        "sample_spacing": args.sample_spacing,
        "sample_order": sample_order,
        "draw_scale_bar": args.scale_bar_output is None,  # Don't draw if separate file
    }
    image_height = draw_reads_svg(sorted_reads, colors, args.output, config)

    # Generate PNG directly (bypasses SVG element limits)
    png_output = args.output.rsplit('.', 1)[0] + '.png'
    print(f"\nGenerating PNG: {png_output}")
    draw_reads_png(sorted_reads, colors, png_output, config)

    # Generate separate scale bar if requested
    if args.scale_bar_output:
        draw_scale_bar_svg(
            args.scale_bar_output,
            image_height,
            args.top_margin,
            args.ratio,
            args.background,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
