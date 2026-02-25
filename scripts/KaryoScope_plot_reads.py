#!/usr/bin/env python3
"""
KaryoScope_plot_reads.py

Visualize telomeric reads with region (satellite) features.
Supports vertical bars (default) or horizontal bars (--horizontal).

Input modes:
  --samples + --results-dir    Load from KaryoScope results directory
  --bed                        Load from BED files directly
  --clusters + --cluster-prefix  Load all reads from specific clusters

Features:
  --legend         Auto-filtered color legend
  --horizontal     Reads as horizontal bars (1 read per row)
  --animate        Generate panning MP4 animation

Usage:
    python KaryoScope_plot_reads.py \\
        --samples NHA_p1 NHA_E6E7_3_PDL48 ... \\
        --results-dir results \\
        --colors KS_human_CHM13.region.colors.txt \\
        --output output.svg
"""

import argparse
import gzip
import os
import subprocess
from collections import defaultdict

import drawsvg as draw
from PIL import Image, ImageDraw, ImageFont
import numpy as np


def hex_to_rgba(hex_color):
    """Convert hex color to RGBA tuple.

    Supports both 6-character (#RRGGBB) and 8-character (#RRGGBBAA) hex codes.
    For 6-character codes, alpha defaults to 255 (fully opaque).
    """
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 6:
        r, g, b = (int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        return (r, g, b, 255)
    elif len(hex_color) == 8:
        r, g, b, a = (int(hex_color[i:i+2], 16) for i in (0, 2, 4, 6))
        return (r, g, b, a)
    else:
        # Fallback for invalid formats
        return (255, 255, 255, 255)


def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple (for backward compatibility)."""
    rgba = hex_to_rgba(hex_color)
    return rgba[:3]


MAX_PNG_DIMENSION = 32000  # Stay under common image library limits


def calculate_png_params(reads, config):
    """Calculate optimal PNG parameters within dimension limits.

    For panning animations, we prioritize HEIGHT resolution (for zoom quality)
    over width. Width dimensions (bar_width, spacing) only need to scale
    enough to maintain visual clarity, not for zoom.

    Returns:
        tuple: (height_scale, width_scale)
        - height_scale: Scale factor for vertical dimensions (ratio, margins)
        - width_scale: Scale factor for horizontal dimensions (bar, spacing)
    """
    base_ratio = config["ratio"]
    requested_scale = config.get("png_scale", 8.0)
    bar_width = config["bar_width"]
    read_spacing = config["read_spacing"]
    sample_spacing = config["sample_spacing"]
    top_margin = config["top_margin"]
    bottom_margin = config["bottom_margin"]

    # Calculate base dimensions (at scale=1)
    max_length = max(r[2] for r in reads)
    base_height = int(max_length * base_ratio) + top_margin + bottom_margin

    # Count samples for width calculation
    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1
    num_samples = len(sample_read_counts)
    total_reads = len(reads)

    base_width = (
        20  # Left margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + 50  # Right margin
    )

    # Height scale: allow full requested scale if within limits
    max_height_scale = MAX_PNG_DIMENSION / base_height
    height_scale = min(requested_scale, max_height_scale)

    # Width scale: scale less aggressively, just enough for clarity
    # For panning, width mainly needs to show reads clearly, not zoom
    # Use sqrt of height_scale as a reasonable width multiplier
    desired_width_scale = min(height_scale, max(1.0, height_scale ** 0.5))
    max_width_scale = MAX_PNG_DIMENSION / base_width
    width_scale = min(desired_width_scale, max_width_scale)

    if height_scale < requested_scale:
        print(f"  Note: Capping height scale from {requested_scale}x to {height_scale:.2f}x "
              f"(dimension limit: {MAX_PNG_DIMENSION}px)")

    return height_scale, width_scale


def draw_reads_png(reads, colors, output_path, config):
    """Draw reads directly to PNG using Pillow at high resolution.

    This bypasses SVG element limits and produces images suitable for
    high-zoom panning animations. Height is scaled more aggressively than
    width to prioritize zoom quality on features.
    """
    base_bar_width = config["bar_width"]
    base_ratio = config["ratio"]
    base_top_margin = config["top_margin"]
    base_left_margin = config["left_margin"]
    base_bottom_margin = config["bottom_margin"]
    base_read_spacing = config["read_spacing"]
    base_sample_spacing = config["sample_spacing"]
    background = config["background"]
    sample_order = config["sample_order"]

    # Calculate optimal scales (height may be higher than width)
    height_scale, width_scale = calculate_png_params(reads, config)

    # Scale dimensions - height uses full scale, width uses reduced scale
    bar_width = max(1, int(base_bar_width * width_scale))
    ratio = base_ratio * height_scale  # Full scale for feature detail
    top_margin = int(base_top_margin * height_scale)
    left_margin = int(base_left_margin * width_scale) if config.get("draw_scale_bar", True) else int(20 * width_scale)
    bottom_margin = int(base_bottom_margin * height_scale)
    read_spacing = max(1, int(base_read_spacing * width_scale))
    sample_spacing = int(base_sample_spacing * width_scale)

    # Calculate dimensions
    max_length = max(r[2] for r in reads)
    max_height_px = int(max_length * ratio)

    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1

    total_reads = len(reads)
    num_samples = len(sample_read_counts)

    effective_left_margin = int(20 * width_scale)  # Reduced since scale bar is separate
    image_width = (
        effective_left_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + int(50 * width_scale)
    )
    image_height = top_margin + max_height_px + bottom_margin

    print(f"  PNG scale: height={height_scale:.2f}x, width={width_scale:.2f}x")
    print(f"  Target dimensions: {image_width} x {image_height} pixels")

    # Create image in RGBA mode to support transparency
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    img = Image.new('RGBA', (image_width, image_height), bg_color)
    draw_ctx = ImageDraw.Draw(img)

    # Try to load font at scaled size
    font_size = int(16 * height_scale)
    try:
        font = ImageFont.truetype(
            "/Users/fbarthel/Documents/Barthel-Custom-Powerpoint-Theme/fonts/BasicSans-Regular.otf",
            font_size
        )
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("Arial", font_size)
        except (OSError, IOError):
            font = ImageFont.load_default()

    # Track positions
    current_x = effective_left_margin
    current_sample = None
    sample_x_start = {}
    sample_x_end = {}

    # Draw reads (using numpy for faster rectangle drawing on very large images)
    print(f"  Drawing {total_reads} reads...")
    img_array = np.array(img)

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_x_end[current_sample] = current_x
            current_x += sample_spacing

        if sample not in sample_x_start:
            sample_x_start[sample] = current_x

        current_sample = sample

        for start, end, feature in features:
            y_start = top_margin + int(start * ratio)
            y_end = top_margin + int(end * ratio)
            height = max(1, y_end - y_start)
            color_hex = colors.get(feature, "#ffffff")
            color_rgba = hex_to_rgba(color_hex)

            # Draw directly to numpy array for speed
            x1 = max(0, current_x)
            x2 = min(image_width, current_x + bar_width)
            y1 = max(0, y_start)
            y2 = min(image_height, y_start + height)
            if x1 < x2 and y1 < y2:
                alpha = color_rgba[3] / 255.0
                if alpha >= 1.0:
                    # Fully opaque - direct assignment
                    img_array[y1:y2, x1:x2] = color_rgba
                else:
                    # Alpha blending: new = alpha * fg + (1-alpha) * bg
                    fg = np.array(color_rgba[:3], dtype=np.float32)
                    bg = img_array[y1:y2, x1:x2, :3].astype(np.float32)
                    blended = alpha * fg + (1 - alpha) * bg
                    img_array[y1:y2, x1:x2, :3] = blended.astype(np.uint8)
                    img_array[y1:y2, x1:x2, 3] = 255  # Keep alpha channel at full

        current_x += bar_width + read_spacing

    if current_sample:
        sample_x_end[current_sample] = current_x

    # Convert back to PIL Image for text rendering
    img = Image.fromarray(img_array)
    draw_ctx = ImageDraw.Draw(img)

    # Draw labels and lines (unless --no-header)
    if not config.get("no_header", False):
        label_interval = 300
        read_width = bar_width + read_spacing
        line_y = top_margin - int(5 * height_scale)
        label_y = top_margin - int(25 * height_scale)

        for sample in sample_order:
            if sample in sample_x_start and sample in sample_x_end:
                # White line
                draw_ctx.line(
                    [(sample_x_start[sample], line_y),
                     (sample_x_end[sample], line_y)],
                    fill=text_color, width=max(1, int(2 * height_scale))
                )

                # Labels
                label_text = sample.replace("_", " ")
                num_reads = sample_read_counts[sample]

                for i in range(0, num_reads, label_interval):
                    x_pos = sample_x_start[sample] + (i * read_width)
                    draw_ctx.text((x_pos, label_y), label_text, fill=text_color, font=font)

        # Draw separator lines
        for sample in sample_order[:-1]:
            if sample in sample_x_end:
                sep_x = sample_x_end[sample] + sample_spacing // 2
                # Create semi-transparent line by drawing in gray
                sep_color = (77, 77, 77) if background == "black" else (200, 200, 200)
                draw_ctx.line(
                    [(sep_x, top_margin), (sep_x, top_margin + max_height_px)],
                    fill=sep_color, width=1
                )

    # Save PNG - blend alpha with background to produce RGB (faster for animation)
    print(f"  Saving PNG...")
    if img.mode == 'RGBA':
        # Alpha composite onto solid background
        bg_img = Image.new('RGB', img.size, bg_color[:3])
        bg_img.paste(img, mask=img.split()[3])  # Use alpha channel as mask
        bg_img.save(output_path, optimize=True)
    else:
        img.save(output_path, optimize=True)
    file_size = os.path.getsize(output_path)
    print(f"Saved PNG: {output_path}")
    print(f"  Dimensions: {image_width} x {image_height} pixels")
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")

    return image_height


def calculate_horizontal_png_params(reads, config):
    """Calculate optimal PNG parameters for horizontal mode within dimension limits.

    Returns:
        tuple: (height_scale, width_scale)
    """
    base_ratio = config["ratio"]
    requested_scale = config.get("png_scale", 8.0)
    bar_width = config["bar_width"]
    read_spacing = config["read_spacing"]
    sample_spacing = config["sample_spacing"]
    left_margin = config["left_margin"]
    top_margin = config["top_margin"]

    # Base dimensions (at scale=1)
    max_length = max(r[2] for r in reads)
    base_width = int(max_length * base_ratio) + left_margin + 50  # right margin

    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1
    num_samples = len(sample_read_counts)
    total_reads = len(reads)

    base_height = (
        top_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + 50  # bottom margin
    )

    # Width scale: allow full requested scale for feature detail
    max_width_scale = MAX_PNG_DIMENSION / base_width
    width_scale = min(requested_scale, max_width_scale)

    # Height scale: scale less aggressively
    desired_height_scale = min(width_scale, max(1.0, width_scale ** 0.5))
    max_height_scale = MAX_PNG_DIMENSION / base_height
    height_scale = min(desired_height_scale, max_height_scale)

    if width_scale < requested_scale:
        print(f"  Note: Capping width scale from {requested_scale}x to {width_scale:.2f}x "
              f"(dimension limit: {MAX_PNG_DIMENSION}px)")

    return height_scale, width_scale


def draw_reads_horizontal_png(reads, colors, output_path, config):
    """Draw reads as horizontal bars (1 read per row) directly to PNG.

    Each read is a horizontal bar with features laid out left-to-right.
    Reads are stacked vertically, grouped by sample/cluster.
    """
    base_bar_width = config["bar_width"]
    base_ratio = config["ratio"]
    base_top_margin = config["top_margin"]
    base_left_margin = config["left_margin"]
    base_read_spacing = config["read_spacing"]
    base_sample_spacing = config["sample_spacing"]
    background = config["background"]
    sample_order = config["sample_order"]

    # Calculate optimal scales
    height_scale, width_scale = calculate_horizontal_png_params(reads, config)

    # Scale dimensions - width uses full scale for feature detail, height uses reduced
    bar_width = max(1, int(base_bar_width * height_scale))
    ratio = base_ratio * width_scale  # Full scale for feature detail along x-axis
    top_margin = int(base_top_margin * height_scale)
    left_margin = int(base_left_margin * height_scale)
    read_spacing = max(1, int(base_read_spacing * height_scale))
    sample_spacing = int(base_sample_spacing * height_scale)

    # Calculate dimensions
    max_length = max(r[2] for r in reads)
    max_width_px = int(max_length * ratio)

    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1

    total_reads = len(reads)
    num_samples = len(sample_read_counts)

    image_width = left_margin + max_width_px + int(50 * width_scale)
    image_height = (
        top_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + int(50 * height_scale)
    )

    print(f"  PNG scale: height={height_scale:.2f}x, width={width_scale:.2f}x")
    print(f"  Target dimensions: {image_width} x {image_height} pixels")

    # Create image
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    img = Image.new('RGBA', (image_width, image_height), bg_color)

    # Try to load font
    font_size = int(16 * height_scale)
    try:
        font = ImageFont.truetype(
            "/Users/fbarthel/Documents/Barthel-Custom-Powerpoint-Theme/fonts/BasicSans-Regular.otf",
            font_size
        )
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("Arial", font_size)
        except (OSError, IOError):
            font = ImageFont.load_default()

    # Track positions
    current_y = top_margin
    current_sample = None
    sample_y_start = {}
    sample_y_end = {}

    # Draw reads using numpy for speed
    print(f"  Drawing {total_reads} reads...")
    img_array = np.array(img)

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_y_end[current_sample] = current_y
            current_y += sample_spacing

        if sample not in sample_y_start:
            sample_y_start[sample] = current_y

        current_sample = sample

        for start, end, feature in features:
            x_start = left_margin + int(start * ratio)
            x_end = left_margin + int(end * ratio)
            width = max(1, x_end - x_start)
            color_hex = colors.get(feature, "#ffffff")
            color_rgba = hex_to_rgba(color_hex)

            # Draw directly to numpy array
            x1 = max(0, x_start)
            x2 = min(image_width, x_start + width)
            y1 = max(0, current_y)
            y2 = min(image_height, current_y + bar_width)
            if x1 < x2 and y1 < y2:
                alpha = color_rgba[3] / 255.0
                if alpha >= 1.0:
                    img_array[y1:y2, x1:x2] = color_rgba
                else:
                    fg = np.array(color_rgba[:3], dtype=np.float32)
                    bg = img_array[y1:y2, x1:x2, :3].astype(np.float32)
                    blended = alpha * fg + (1 - alpha) * bg
                    img_array[y1:y2, x1:x2, :3] = blended.astype(np.uint8)
                    img_array[y1:y2, x1:x2, 3] = 255

        current_y += bar_width + read_spacing

    if current_sample:
        sample_y_end[current_sample] = current_y

    # Convert back to PIL for text rendering
    img = Image.fromarray(img_array)
    draw_ctx = ImageDraw.Draw(img)

    # Draw labels and separator lines (unless --no-header)
    if not config.get("no_header", False):
        label_interval = 300
        read_height = bar_width + read_spacing

        for sample in sample_order:
            if sample in sample_y_start and sample in sample_y_end:
                # Vertical line at left spanning this sample
                line_x = left_margin - int(5 * height_scale)
                draw_ctx.line(
                    [(line_x, sample_y_start[sample]),
                     (line_x, sample_y_end[sample])],
                    fill=text_color, width=max(1, int(2 * height_scale))
                )

                # Labels (rotated vertically on left margin)
                label_text = sample.replace("_", " ")
                num_reads_in_sample = sample_read_counts[sample]

                for i in range(0, num_reads_in_sample, label_interval):
                    y_pos = sample_y_start[sample] + (i * read_height)
                    # Draw rotated text using a temporary image
                    bbox = font.getbbox(label_text)
                    tw = bbox[2] - bbox[0]
                    th = bbox[3] - bbox[1]
                    txt_img = Image.new('RGBA', (tw + 4, th + 4), (0, 0, 0, 0))
                    txt_draw = ImageDraw.Draw(txt_img)
                    txt_draw.text((2, 2), label_text, fill=text_color + (255,), font=font)
                    txt_img = txt_img.rotate(90, expand=True)
                    paste_x = max(0, line_x - txt_img.width - int(5 * height_scale))
                    paste_y = y_pos
                    if paste_y + txt_img.height <= image_height:
                        img.paste(txt_img, (paste_x, paste_y), txt_img)

        # Draw separator lines between samples
        for sample in sample_order[:-1]:
            if sample in sample_y_end:
                sep_y = sample_y_end[sample] + sample_spacing // 2
                sep_color = (77, 77, 77) if background == "black" else (200, 200, 200)
                draw_ctx = ImageDraw.Draw(img)  # Refresh after paste
                draw_ctx.line(
                    [(left_margin, sep_y), (left_margin + max_width_px, sep_y)],
                    fill=sep_color, width=1
                )

    # Draw horizontal scale bar at top
    if config.get("draw_scale_bar", True):
        draw_ctx = ImageDraw.Draw(img)
        scale_bar_bp = 10000
        scale_bar_width_px = int(scale_bar_bp * ratio)
        if scale_bar_width_px > max_width_px:
            scale_bar_bp = 5000
            scale_bar_width_px = int(scale_bar_bp * ratio)

        sb_x = left_margin
        sb_y = top_margin - int(30 * height_scale)
        if sb_y < 5:
            sb_y = 5
        # Horizontal bar
        draw_ctx.rectangle(
            [sb_x, sb_y, sb_x + scale_bar_width_px, sb_y + max(2, int(3 * height_scale))],
            fill=text_color
        )
        # Tick marks
        tick_h = max(2, int(4 * height_scale))
        draw_ctx.line([(sb_x, sb_y - tick_h), (sb_x, sb_y + tick_h)],
                      fill=text_color, width=1)
        draw_ctx.line([(sb_x + scale_bar_width_px, sb_y - tick_h),
                       (sb_x + scale_bar_width_px, sb_y + tick_h)],
                      fill=text_color, width=1)
        # Label
        label = f"{scale_bar_bp // 1000} Kbp"
        label_x = sb_x + scale_bar_width_px // 2
        label_y = sb_y - int(15 * height_scale)
        if label_y < 2:
            label_y = 2
        draw_ctx.text((label_x, label_y), label, fill=text_color, font=font, anchor="mt")

    # Save PNG
    print(f"  Saving PNG...")
    if img.mode == 'RGBA':
        bg_img = Image.new('RGB', img.size, bg_color[:3])
        bg_img.paste(img, mask=img.split()[3])
        bg_img.save(output_path, optimize=True)
    else:
        img.save(output_path, optimize=True)
    file_size = os.path.getsize(output_path)
    print(f"Saved PNG: {output_path}")
    print(f"  Dimensions: {image_width} x {image_height} pixels")
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")

    return image_width


def draw_reads_horizontal_svg(reads, colors, output_path, config):
    """Draw reads as horizontal bars (1 read per row) in SVG."""
    bar_width = config["bar_width"]
    ratio = config["ratio"]
    top_margin = config["top_margin"]
    left_margin = config["left_margin"]
    read_spacing = config["read_spacing"]
    sample_spacing = config["sample_spacing"]
    background = config["background"]
    sample_order = config["sample_order"]
    draw_scale = config.get("draw_scale_bar", True)

    max_length = max(r[2] for r in reads)
    max_width_px = int(max_length * ratio)

    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1

    total_reads = len(reads)
    num_samples = len(sample_read_counts)

    image_width = left_margin + max_width_px + 50
    image_height = (
        top_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + 50
    )

    d = draw.Drawing(image_width, image_height, id_prefix="trh")
    text_color = "#ffffff" if background == "black" else "#000000"
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background))

    # Draw horizontal scale bar at top
    if draw_scale:
        scale_bar_bp = 10000
        scale_bar_width_px = int(scale_bar_bp * ratio)
        if scale_bar_width_px > max_width_px:
            scale_bar_bp = 5000
            scale_bar_width_px = int(scale_bar_bp * ratio)

        sb_x = left_margin
        sb_y = top_margin - 30
        d.append(draw.Rectangle(sb_x, sb_y, scale_bar_width_px, 3, fill=text_color))
        d.append(draw.Line(sb_x, sb_y - 4, sb_x, sb_y + 4,
                           stroke=text_color, stroke_width=1))
        d.append(draw.Line(sb_x + scale_bar_width_px, sb_y - 4,
                           sb_x + scale_bar_width_px, sb_y + 4,
                           stroke=text_color, stroke_width=1))
        label = f"{scale_bar_bp // 1000} Kbp"
        d.append(draw.Text(label, 12, sb_x + scale_bar_width_px / 2, sb_y - 8,
                           fill=text_color, text_anchor="middle",
                           font_family="Basic Sans"))

    current_y = top_margin
    current_sample = None
    sample_y_start = {}
    sample_y_end = {}

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_y_end[current_sample] = current_y
            current_y += sample_spacing

        if sample not in sample_y_start:
            sample_y_start[sample] = current_y

        current_sample = sample

        for start, end, feature in features:
            x_start = left_margin + int(start * ratio)
            width = max(1, int((end - start) * ratio))
            color = colors.get(feature, "#ffffff")
            d.append(draw.Rectangle(x_start, current_y, width, bar_width, fill=color))

        current_y += bar_width + read_spacing

    if current_sample:
        sample_y_end[current_sample] = current_y

    # Draw labels and separators (unless --no-header)
    if not config.get("no_header", False):
        font_family = "Basic Sans"
        label_interval = 300
        read_height = bar_width + read_spacing

        for sample in sample_order:
            if sample in sample_y_start and sample in sample_y_end:
                # Vertical line at left
                line_x = left_margin - 5
                d.append(draw.Line(line_x, sample_y_start[sample],
                                   line_x, sample_y_end[sample],
                                   stroke=text_color, stroke_width=2))

                label_text = sample.replace("_", " ")
                num_reads_in_sample = sample_read_counts[sample]

                for i in range(0, num_reads_in_sample, label_interval):
                    y_pos = sample_y_start[sample] + (i * read_height)
                    lx = line_x - 10
                    ly = y_pos + 50
                    d.append(draw.Text(
                        label_text, 16, lx, ly,
                        fill=text_color, text_anchor="middle",
                        font_family=font_family,
                        transform=f"rotate(-90, {lx}, {ly})",
                    ))

        for sample in sample_order[:-1]:
            if sample in sample_y_end:
                sep_y = sample_y_end[sample] + sample_spacing / 2
                d.append(draw.Line(
                    left_margin, sep_y, left_margin + max_width_px, sep_y,
                    stroke=text_color, stroke_width=0.5, stroke_opacity=0.3,
                ))

    d.save_svg(output_path)
    print(f"\nSaved: {output_path}")
    print(f"  Dimensions: {image_width} x {image_height} pixels")

    return image_width


def draw_legend_png(features_used, colors, config, horizontal=False):
    """Render an auto-filtered color legend as a PIL Image.

    Args:
        features_used: set of feature names that appeared in the plotted reads
        colors: {feature_name: hex_color} mapping
        config: dict with background, etc.
        horizontal: if True, render single-column layout (for horizontal mode)

    Returns:
        PIL.Image: legend image ready for compositing
    """
    background = config["background"]
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)

    swatch_size = 14
    font_size = 14
    padding = 10
    row_height = swatch_size + 6

    # Load font
    try:
        font = ImageFont.truetype(
            "/Users/fbarthel/Documents/Barthel-Custom-Powerpoint-Theme/fonts/BasicSans-Regular.otf",
            font_size
        )
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("Arial", font_size)
        except (OSError, IOError):
            font = ImageFont.load_default()

    # Filter to only features in use, skip _specific duplicates
    filtered = {}
    for feat in sorted(features_used):
        base = feat.replace("_specific", "")
        if base in filtered:
            continue
        color_hex = colors.get(feat, colors.get(base, None))
        if color_hex:
            filtered[base] = color_hex

    if not filtered:
        return None

    items = list(filtered.items())

    # Measure max label width
    max_label_w = 0
    for name, _ in items:
        bbox = font.getbbox(name)
        max_label_w = max(max_label_w, bbox[2] - bbox[0])

    col_width = swatch_size + 8 + max_label_w + padding * 2

    if horizontal:
        # Single-column layout (placed to the right of horizontal reads)
        num_cols = 1
        num_rows = len(items)
    else:
        # Multi-column layout (placed below vertical reads)
        max_cols = 6
        num_cols = min(max_cols, max(1, len(items)))
        num_rows = (len(items) + num_cols - 1) // num_cols

    legend_width = num_cols * col_width + padding * 2
    legend_height = num_rows * row_height + padding * 2

    legend_img = Image.new('RGBA', (legend_width, legend_height), bg_color)
    ctx = ImageDraw.Draw(legend_img)

    for idx, (name, color_hex) in enumerate(items):
        if horizontal:
            col = 0
            row = idx
        else:
            col = idx // num_rows
            row = idx % num_rows

        x = padding + col * col_width
        y = padding + row * row_height

        rgba = hex_to_rgba(color_hex)
        ctx.rectangle([x, y, x + swatch_size, y + swatch_size], fill=rgba)
        ctx.text((x + swatch_size + 8, y), name, fill=text_color, font=font)

    return legend_img


def composite_legend(png_path, legend_img, position="below"):
    """Composite a legend image onto a rendered PNG.

    Args:
        png_path: path to existing PNG to modify in-place
        legend_img: PIL.Image of the legend
        position: "below" or "right"
    """
    if legend_img is None:
        return

    main_img = Image.open(png_path)
    if main_img.mode != 'RGB':
        main_img = main_img.convert('RGB')

    lw, lh = legend_img.size
    mw, mh = main_img.size

    if position == "right":
        new_width = mw + lw
        new_height = max(mh, lh)
        bg_pixel = main_img.getpixel((0, 0))
        combined = Image.new('RGB', (new_width, new_height), bg_pixel)
        combined.paste(main_img, (0, 0))
        if legend_img.mode == 'RGBA':
            combined.paste(legend_img, (mw, 0), legend_img)
        else:
            combined.paste(legend_img, (mw, 0))
    else:  # below
        new_width = max(mw, lw)
        new_height = mh + lh
        bg_pixel = main_img.getpixel((0, 0))
        combined = Image.new('RGB', (new_width, new_height), bg_pixel)
        combined.paste(main_img, (0, 0))
        if legend_img.mode == 'RGBA':
            combined.paste(legend_img, (0, mh), legend_img)
        else:
            combined.paste(legend_img, (0, mh))

    combined.save(png_path, optimize=True)
    print(f"  Legend composited ({position}): {png_path}")


def load_cluster_reads(cluster_prefix, cluster_ids, results_dir,
                       database, featureset, smoothness, analysis):
    """Load all reads from specified clusters as plot_reads-format tuples.

    Returns:
        tuple: (reads_list, sample_order)
    """
    import pandas as pd

    assignments_file = f"{cluster_prefix}.sequence_assignments.tsv"
    if not os.path.exists(assignments_file):
        print(f"Error: Cluster assignments file not found: {assignments_file}")
        return [], []

    print(f"  Loading cluster assignments: {assignments_file}")
    df = pd.read_csv(assignments_file, sep='\t')
    if 'read' in df.columns:
        df = df.rename(columns={'read': 'sequence'})

    df = df[df['cluster'].isin(cluster_ids)]
    print(f"  Found {len(df)} reads in clusters: {cluster_ids}")

    # Load BED data for all samples represented in these clusters
    samples = df['sample'].unique().tolist()
    all_bed_reads = load_sample_bed_data(
        samples, results_dir, database, featureset, smoothness, analysis)

    # Build lookup by read ID
    bed_lookup = {read_id: feats for (_, read_id, _, feats) in all_bed_reads}

    # Group by cluster
    result = []
    sample_order = []
    for cid in cluster_ids:
        label = f"Cluster {cid}"
        sample_order.append(label)
        cdata = df[df['cluster'] == cid]
        found = 0
        for _, row in cdata.iterrows():
            feats = bed_lookup.get(row['sequence'])
            if not feats:
                continue
            read_length = max(end for _, end, _ in feats)
            result.append((label, row['sequence'], read_length, feats))
            found += 1
        print(f"  {label}: {found} reads with BED data")

    return result, sample_order


def draw_scale_bar_png(image_height, top_margin, ratio, background, output_path,
                       height_scale=1.0):
    """Render a standalone vertical scale bar as PNG for animation overlay.

    Generated at source PNG dimensions with scaled ratio. The animation function
    will crop+resize this to match the viewport, keeping the bar proportional
    to the visible reads.
    """
    width = 60
    bg_color = (0, 0, 0) if background == "black" else (255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)

    scale_bar_bp = 10000
    scale_bar_height = int(scale_bar_bp * ratio)
    max_height_px = image_height - top_margin - 40
    if scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    img = Image.new('RGB', (width, image_height), bg_color)
    draw_ctx = ImageDraw.Draw(img)

    # Load font — scale with height_scale so text survives downscale to viewport
    font_size = max(14, int(14 * height_scale))
    try:
        font = ImageFont.truetype(
            "/Users/fbarthel/Documents/Barthel-Custom-Powerpoint-Theme/fonts/BasicSans-Regular.otf",
            font_size
        )
    except (OSError, IOError):
        try:
            font = ImageFont.truetype("Arial", font_size)
        except (OSError, IOError):
            font = ImageFont.load_default()

    x = 45
    y = top_margin
    bar_w = max(2, int(3 * height_scale))
    tick_ext = max(3, int(5 * height_scale))

    # Vertical bar
    draw_ctx.rectangle([x, y, x + bar_w, y + scale_bar_height], fill=text_color)
    # Tick marks
    draw_ctx.line([(x - tick_ext, y), (x + tick_ext, y)], fill=text_color, width=1)
    draw_ctx.line([(x - tick_ext, y + scale_bar_height),
                   (x + tick_ext, y + scale_bar_height)], fill=text_color, width=1)

    # Rotated label
    label = f"{scale_bar_bp // 1000} Kbp"
    bbox = font.getbbox(label)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    txt_img = Image.new('RGBA', (tw + 4, th + 4), (0, 0, 0, 0))
    txt_draw = ImageDraw.Draw(txt_img)
    txt_draw.text((2, 2), label, fill=text_color + (255,), font=font)
    txt_img = txt_img.rotate(90, expand=True)
    paste_x = max(0, x - tick_ext - txt_img.width - 2)
    paste_y = y + (scale_bar_height - txt_img.height) // 2
    img.paste(txt_img, (paste_x, paste_y), txt_img)

    img.save(output_path)
    print(f"  Scale bar saved for animation: {output_path}")
    print(f"  Dimensions: {width} x {image_height} pixels")


def _run_animation(png_path, num_reads, args, config=None,
                    legend_path=None, scale_bar_path=None):
    """Generate panning animation MP4 from a rendered PNG."""
    import sys as _sys

    script_dir = os.path.dirname(os.path.abspath(__file__))
    _sys.path.insert(0, script_dir)
    from create_panning_animation import (
        create_horizontal_panning, create_vertical_panning,
        create_adaptive_horizontal_panning,
    )

    base = png_path.rsplit('.', 1)[0]
    mp4_path = f"{base}.mp4"

    # Auto-detect pan direction from read orientation if not specified
    direction = args.animate_direction
    if direction is None:
        direction = "vertical" if args.horizontal else "horizontal"

    # Auto-calculate duration
    duration = args.animate_duration or (num_reads / args.animate_reads_per_second)

    # Viewport defaults
    if args.animate_viewport:
        vw, vh = [int(x) for x in args.animate_viewport.split('x')]
    elif direction == "vertical":
        vw, vh = 640, 864
    else:
        vw, vh = 1920, 1080

    bg = args.background
    # Use built-in legend if available, otherwise fall back to --animate-legend
    legend = legend_path or args.animate_legend
    scale_bar = scale_bar_path

    # Compute the scaled ratio (px/bp in the source PNG)
    # The source PNG uses ratio * height_scale for vertical dimensions
    height_scale = config.get("png_scale", 8.0) if config else 1.0
    scaled_ratio = args.ratio * height_scale

    print(f"\nGenerating animation: {mp4_path}")
    print(f"  Direction: {direction}, Duration: {duration:.1f}s, FPS: {args.animate_fps}")
    print(f"  Viewport: {vw}x{vh}, Zoom: {args.animate_zoom}")
    if legend:
        print(f"  Legend: {legend}")
    if scale_bar:
        print(f"  Scale bar: {scale_bar}")

    if direction == "horizontal":
        # Pass scaled margins/ratio matching the source PNG's coordinate space
        scaled_top_margin = int(args.top_margin * height_scale)
        scaled_left_margin = int(args.left_margin * max(1.0, height_scale ** 0.5))
        create_adaptive_horizontal_panning(
            png_path, mp4_path, duration, args.animate_fps, vw, vh,
            legend_path=legend, background=bg,
            top_margin=scaled_top_margin,
            left_margin=scaled_left_margin,
            ratio=scaled_ratio,
            max_zoom=1.0, scale_bar_padding=10)
    else:
        create_vertical_panning(
            png_path, mp4_path, duration, args.animate_fps, vw, vh,
            scale_bar_path=scale_bar, legend_path=legend,
            background=bg, scale_bar_padding=10,
            ratio=scaled_ratio)

    print(f"  Animation: {mp4_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize telomeric reads with region features"
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
        "--analysis",
        default="telogator",
        help="Analysis type subdirectory name (default: telogator). "
             "Used to construct path: {sample}/{analysis}/1/KaryoScope/{database}/",
    )
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
    parser.add_argument(
        "--orient-chromosome-top",
        action="store_true",
        help="Reorient reads so chromosome-specific features are always at the top",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="If set, split reads into batches of this size and output separate SVGs",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Skip drawing sample labels, horizontal lines, and separator lines",
    )
    parser.add_argument(
        "--no-scale-bar",
        action="store_true",
        help="Skip drawing the scale bar in the left margin",
    )
    parser.add_argument(
        "--format",
        choices=["svg", "png", "both"],
        default="svg",
        help="Output format: svg, png, or both (default: svg)",
    )
    parser.add_argument(
        "--png-scale",
        type=float,
        default=8.0,
        help="Resolution multiplier for PNG output (default: 8.0 for 8x zoom support). "
             "Higher values = sharper at high zoom but larger files.",
    )

    # ── Orientation ──
    parser.add_argument(
        "--horizontal",
        action="store_true",
        help="Draw reads as horizontal bars (1 read per row). "
             "Default is vertical (1 read per column).",
    )

    # ── Legend ──
    parser.add_argument(
        "--legend",
        action="store_true",
        help="Draw a color legend. Auto-filtered to features present in plotted reads. "
             "Vertical mode: legend below in multi-column layout. "
             "Horizontal mode: legend to the right in single-column layout.",
    )

    # ── Cluster input mode ──
    parser.add_argument(
        "--cluster-prefix",
        default=None,
        help="Prefix for cluster analysis outputs. Auto-discovers "
             "{prefix}.sequence_assignments.tsv. Use with --clusters.",
    )
    parser.add_argument(
        "--clusters",
        default=None,
        help="Comma-separated cluster IDs to plot (e.g., '22,69,112'). "
             "Requires --cluster-prefix and --results-dir.",
    )

    # ── Animation ──
    parser.add_argument(
        "--animate",
        action="store_true",
        help="Generate panning animation MP4 from the rendered PNG. "
             "Implies --format png.",
    )
    parser.add_argument(
        "--animate-direction",
        default=None,
        choices=["horizontal", "vertical"],
        help="Pan direction (default: auto from read orientation — "
             "horizontal pan for vertical reads, vertical pan for horizontal reads)",
    )
    parser.add_argument(
        "--animate-duration",
        type=float,
        default=None,
        help="Duration in seconds (auto-calculated from read count if omitted)",
    )
    parser.add_argument(
        "--animate-reads-per-second",
        type=float,
        default=20.625,
        help="Scroll rate for auto-duration (default: 20.625)",
    )
    parser.add_argument(
        "--animate-zoom",
        default="fixed",
        choices=["fixed", "adaptive"],
        help="Zoom mode (default: fixed)",
    )
    parser.add_argument(
        "--animate-crop-ratio",
        type=float,
        default=0.5,
        help="Height crop ratio for fixed zoom (default: 0.5)",
    )
    parser.add_argument(
        "--animate-fps",
        type=int,
        default=30,
        help="Frames per second (default: 30)",
    )
    parser.add_argument(
        "--animate-viewport",
        default=None,
        help="Viewport as WxH (default: auto based on direction)",
    )
    parser.add_argument(
        "--animate-legend",
        default=None,
        help="Separate legend image (PNG/SVG) to composite into animation. "
             "If --legend is used, the built-in legend is composited automatically.",
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


def load_sample_bed_data(samples, results_dir, database, featureset, smoothness, analysis="telogator"):
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
            f"{analysis}/1/KaryoScope",
            database,
            f"{sample}.{analysis}.1.{database}.{featureset}.{smoothness}.features.bed",
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

# Chromosome features used for orientation detection
CHROMOSOME_FEATURES = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',
    'chr1_specific', 'chr2_specific', 'chr3_specific', 'chr4_specific',
    'chr5_specific', 'chr6_specific', 'chr7_specific', 'chr8_specific',
    'chr9_specific', 'chr10_specific', 'chr11_specific', 'chr12_specific',
    'chr13_specific', 'chr14_specific', 'chr15_specific', 'chr16_specific',
    'chr17_specific', 'chr18_specific', 'chr19_specific', 'chr20_specific',
    'chr21_specific', 'chr22_specific', 'chrX_specific', 'chrY_specific',
}


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


def orient_chromosome_top(reads):
    """Reorient reads so chromosome features are at the top (position 0).

    For each read, checks if chromosome-specific features are closer to start or end.
    If closer to end, flips the read coordinates.

    Args:
        reads: List of (sample, read_id, read_length, features) tuples

    Returns:
        List of reoriented reads
    """
    oriented_reads = []
    flipped_count = 0

    for sample, read_id, read_length, features in reads:
        # Find chromosome feature positions
        chromosome_positions = []
        for start, end, feature in features:
            if feature in CHROMOSOME_FEATURES:
                chromosome_positions.extend([start, end])

        if chromosome_positions:
            # Calculate average chromosome position
            avg_chromosome_pos = sum(chromosome_positions) / len(chromosome_positions)
            midpoint = read_length / 2

            # If chromosome is in second half, flip the read
            if avg_chromosome_pos > midpoint:
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
            # No chromosome features, keep as-is
            oriented_reads.append((sample, read_id, read_length, features))

    print(f"  Reoriented {flipped_count} of {len(reads)} reads (chromosome now at top)")
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

    # Draw horizontal line, sample labels, and separators (unless --no-header)
    if not config.get("no_header", False):
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

    print("KaryoScope Read Visualization")
    print("=" * 50)

    # If --animate is set, ensure PNG format
    fmt = args.format
    if args.animate and fmt == "svg":
        print("  --animate requires PNG output, switching format to 'both'")
        fmt = "both"

    # Validate input arguments
    if args.bed:
        input_mode = "bed"
    elif args.clusters and args.cluster_prefix:
        input_mode = "cluster"
        if not args.results_dir:
            print("Error: --clusters requires --results-dir")
            return
    elif args.samples and args.results_dir:
        input_mode = "samples"
    else:
        print("Error: Must provide --bed, --samples+--results-dir, "
              "or --clusters+--cluster-prefix+--results-dir")
        return

    # Load color mapping
    print(f"\nLoading colors from: {args.colors}")
    colors = load_color_mapping(args.colors)
    print(f"  Loaded {len(colors)} color mappings")

    # Load BED data
    if input_mode == "bed":
        print(f"\nLoading feature data from BED files...")
        reads, sample_order = load_bed_files_direct(args.bed)
    elif input_mode == "cluster":
        cluster_ids = [int(c.strip()) for c in args.clusters.split(',')]
        print(f"\nLoading reads for clusters: {cluster_ids}")
        reads, sample_order = load_cluster_reads(
            args.cluster_prefix, cluster_ids, args.results_dir,
            args.database, args.featureset, args.smoothness, args.analysis)
    else:
        print(f"\nLoading feature data from: {args.results_dir}")
        reads = load_sample_bed_data(
            args.samples,
            args.results_dir,
            args.database,
            args.featureset,
            args.smoothness,
            args.analysis,
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

    # Orient reads so chromosome is at top (if requested)
    if args.orient_chromosome_top:
        print(f"\nOrienting reads (chromosome at top)...")
        reads = orient_chromosome_top(reads)

    # Sort reads
    print(f"\nSorting reads by sample order, then by length (descending)")
    sorted_reads = sort_reads(reads, sample_order)

    # Build config
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
        "draw_scale_bar": args.scale_bar_output is None and not args.no_scale_bar,
        "no_header": args.no_header,
        "png_scale": args.png_scale,
    }

    # Determine output paths
    base_output = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
    svg_output = args.output if args.output.endswith('.svg') else f"{base_output}.svg"
    png_output = f"{base_output}.png"

    image_size = None  # height for vertical, width for horizontal

    if args.batch_size:
        num_batches = (len(sorted_reads) + args.batch_size - 1) // args.batch_size
        print(f"\nGenerating {num_batches} batched outputs ({args.batch_size} reads each)")

        for i in range(0, len(sorted_reads), args.batch_size):
            batch = sorted_reads[i:i + args.batch_size]
            batch_num = (i // args.batch_size) + 1
            batch_base = f"{base_output}.batch{batch_num:03d}"

            if args.horizontal:
                if fmt in ("svg", "both"):
                    batch_svg = f"{batch_base}.svg"
                    print(f"\nGenerating SVG (horizontal): {batch_svg}")
                    image_size = draw_reads_horizontal_svg(batch, colors, batch_svg, config)
                if fmt in ("png", "both"):
                    batch_png = f"{batch_base}.png"
                    print(f"\nGenerating PNG (horizontal): {batch_png}")
                    image_size = draw_reads_horizontal_png(batch, colors, batch_png, config)
            else:
                if fmt in ("svg", "both"):
                    batch_svg = f"{batch_base}.svg"
                    print(f"\nGenerating SVG: {batch_svg}")
                    image_size = draw_reads_svg(batch, colors, batch_svg, config)
                if fmt in ("png", "both"):
                    batch_png = f"{batch_base}.png"
                    print(f"\nGenerating PNG: {batch_png}")
                    image_size = draw_reads_png(batch, colors, batch_png, config)
    else:
        if args.horizontal:
            if fmt in ("svg", "both"):
                print(f"\nGenerating SVG (horizontal): {svg_output}")
                image_size = draw_reads_horizontal_svg(sorted_reads, colors, svg_output, config)
            if fmt in ("png", "both"):
                print(f"\nGenerating PNG (horizontal): {png_output}")
                image_size = draw_reads_horizontal_png(sorted_reads, colors, png_output, config)
        else:
            if fmt in ("svg", "both"):
                print(f"\nGenerating SVG: {svg_output}")
                image_size = draw_reads_svg(sorted_reads, colors, svg_output, config)
            if fmt in ("png", "both"):
                print(f"\nGenerating PNG: {png_output}")
                image_size = draw_reads_png(sorted_reads, colors, png_output, config)

    # Handle legend: save separately for animation, or composite onto PNG
    legend_png_path = None
    if args.legend and fmt in ("png", "both"):
        features_used = {feat for _, _, _, feats in sorted_reads for _, _, feat in feats}
        legend_img = draw_legend_png(features_used, colors, config,
                                     horizontal=args.horizontal)
        if legend_img is not None:
            if args.animate:
                # Save as separate file — animation function will overlay it
                legend_png_path = f"{base_output}.legend.png"
                if legend_img.mode == 'RGBA':
                    bg_c = (0, 0, 0) if args.background == "black" else (255, 255, 255)
                    bg_img = Image.new('RGB', legend_img.size, bg_c)
                    bg_img.paste(legend_img, mask=legend_img.split()[3])
                    bg_img.save(legend_png_path)
                else:
                    legend_img.save(legend_png_path)
                print(f"  Legend saved for animation: {legend_png_path}")
            else:
                # Composite directly onto the reads image
                position = "right" if args.horizontal else "below"
                composite_legend(png_output, legend_img, position=position)

    # Scale bar: adaptive horizontal panning draws its own; vertical still uses overlay
    scale_bar_png_path = None

    # Generate animation if requested
    if args.animate:
        if fmt not in ("png", "both"):
            # Should not reach here due to earlier correction, but be safe
            print(f"\nGenerating PNG for animation: {png_output}")
            if args.horizontal:
                draw_reads_horizontal_png(sorted_reads, colors, png_output, config)
            else:
                draw_reads_png(sorted_reads, colors, png_output, config)
        _run_animation(png_output, len(sorted_reads), args, config,
                       legend_path=legend_png_path, scale_bar_path=scale_bar_png_path)

    # Generate separate scale bar if requested
    if args.scale_bar_output and image_size and not args.horizontal:
        draw_scale_bar_svg(
            args.scale_bar_output,
            image_size,
            args.top_margin,
            args.ratio,
            args.background,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
