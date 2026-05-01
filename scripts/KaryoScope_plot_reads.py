#!/usr/bin/env python3
"""
KaryoScope_plot_reads.py

Visualize telomeric reads with region (satellite) features as colored bars.
Supports vertical bars (default) or horizontal bars (--horizontal).

Input modes:
  --samples + --results-dir          Load from KaryoScope results directory
  --bed                              Load from BED files directly
  --clusters + --cluster-prefix      Load reads from specific clusters

Output:
  --format svg|png|both              Output format (default: svg)
  --output PATH                      Output file path (extension auto-adjusted)
  --png-scale FACTOR                 Resolution multiplier for PNG (default: 8.0)
  --scale-bar-output PATH            Separate scale bar SVG
  --read-list PATH                   Filter to read IDs from file

Display controls:
  --horizontal                       Reads as horizontal bars (1 read per row)
  --legend                           Auto-filtered color legend (composited or separate)
  --no-header                        Skip sample labels and separator lines
  --no-scale-bar                     Skip scale bar in left margin
  --orient-telomere-top              Reorient reads so telomere is at top
  --orient-chromosome-top            Reorient reads so chromosome is at top
  --batch-size N                     Split output into batches of N reads

Animation:
  --animate                          Generate panning MP4 (implies PNG output)
  --animate-direction                Pan direction (auto-detected by default)
  --animate-viewport WxH             Viewport size
  --animate-duration / --animate-fps Duration and framerate controls

Usage:
    python KaryoScope_plot_reads.py \\
        --samples NHA_p1 NHA_E6E7_3_PDL48 \\
        --results-dir results \\
        --colors KS_human_CHM13.region.colors.txt \\
        --output output.svg
"""

import argparse
import fnmatch
import gzip
import logging
import os
import subprocess
import sys
import time
import re
from collections import defaultdict, namedtuple

import drawsvg as draw
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None
import numpy as np

logger = logging.getLogger(__name__)


def setup_logging(output_path):
    """Configure logging with console (INFO) and file (DEBUG) handlers."""
    base = output_path.rsplit('.', 1)[0] if '.' in output_path else output_path
    log_path = f"{base}.log"

    root = logging.getLogger(__name__)
    root.setLevel(logging.DEBUG)
    root.handlers.clear()

    # Console: INFO, simple format
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(ch)

    # File: DEBUG, detailed format
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    root.addHandler(fh)

    return log_path


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


def _save_png(img, output_path, bg_color):
    """Alpha-composite RGBA image onto solid background and save as PNG."""
    logger.info("  Saving PNG...")
    if img.mode == 'RGBA':
        bg_img = Image.new('RGB', img.size, bg_color[:3])
        bg_img.paste(img, mask=img.split()[3])
        bg_img.save(output_path, optimize=True)
    else:
        img.save(output_path, optimize=True)
    file_size = os.path.getsize(output_path)
    w, h = img.size
    logger.info("Saved PNG: %s", output_path)
    logger.info("  Dimensions: %d x %d pixels", w, h)
    logger.info("  File size: %.1f MB", file_size / 1024 / 1024)


from karyoplot.core.fonts import pil_font as _load_font  # noqa: F401


def _draw_rect_rgba(img_array, x1, y1, x2, y2, color_rgba):
    """Draw an RGBA rectangle onto a numpy image array with alpha blending."""
    if x1 >= x2 or y1 >= y2:
        return
    alpha = color_rgba[3] / 255.0
    if alpha >= 1.0:
        img_array[y1:y2, x1:x2] = color_rgba
    else:
        fg = np.array(color_rgba[:3], dtype=np.float32)
        bg = img_array[y1:y2, x1:x2, :3].astype(np.float32)
        blended = alpha * fg + (1 - alpha) * bg
        img_array[y1:y2, x1:x2, :3] = blended.astype(np.uint8)
        img_array[y1:y2, x1:x2, 3] = 255


from karyoplot.svg.reads import (  # noqa: E402  (placed after argparse cache)
    features_to_pixels_direct,
    rasterize_features,
    smooth_features_to_pixels,
)



def _build_colored_features(features, colors, min_width_exclude):
    """Build colored feature tuples from raw BED features for rasterization.

    Returns list of (bp_start, bp_stop, color_hex, fill_opacity, skip_min).
    """
    result = []
    for start, end, feature in features:
        color_hex = colors.get(feature, "#ffffff")
        skip_min = any(fnmatch.fnmatch(feature, pat) for pat in min_width_exclude)
        result.append((start, end, color_hex, 1.0, skip_min))
    return result


def calculate_png_params(reads, config):
    """Calculate optimal PNG parameters within dimension limits.

    Uses uniform scaling so the PNG is a proportionally faithful copy
    of the SVG, capped by MAX_PNG_DIMENSION on each axis.

    Returns:
        tuple: (height_scale, width_scale)
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

    # Width scale: match height_scale for uniform scaling (PNG matches SVG proportions)
    desired_width_scale = height_scale
    max_width_scale = MAX_PNG_DIMENSION / base_width
    width_scale = min(desired_width_scale, max_width_scale)

    if height_scale < requested_scale:
        logger.info("  Note: Capping height scale from %.1fx to %.2fx (dimension limit: %dpx)",
                    requested_scale, height_scale, MAX_PNG_DIMENSION)

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

    # Use full left margin when display-name labels exist, otherwise reduced
    has_left_labels = (config.get("tier_display_names") or
                       (config.get("heatmap") and config.get("heatmap_display_names")))
    effective_left_margin = left_margin if has_left_labels else int(20 * width_scale)
    image_width = (
        effective_left_margin
        + (total_reads * (bar_width + read_spacing))
        + ((num_samples - 1) * sample_spacing)
        + int(50 * width_scale)
    )
    image_height = top_margin + max_height_px + bottom_margin

    logger.info("  PNG scale: height=%.2fx, width=%.2fx", height_scale, width_scale)
    logger.info("  Target dimensions: %d x %d pixels", image_width, image_height)

    # Create image in RGBA mode to support transparency
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    img = Image.new('RGBA', (image_width, image_height), bg_color)
    draw_ctx = ImageDraw.Draw(img)

    # Try to load font at scaled size
    font_size = int(config.get("font_size", 11) * height_scale)
    font = _load_font(font_size)

    # Track positions
    current_x = effective_left_margin
    current_sample = None
    sample_x_start = {}
    sample_x_end = {}

    # Draw reads (using numpy for faster rectangle drawing on very large images)
    logger.info("  Drawing %d reads...", total_reads)
    img_array = np.array(img)

    feature_mode = config.get("feature_mode", "raw")
    min_feature_width = config.get("min_feature_width", 0.5)
    min_width_exclude = config.get("min_width_exclude", [])
    oversample = config.get("oversample", 1)
    read_border = config.get("read_border", False)
    read_positions = []
    markers = config.get("markers", {})
    marker_scale = config.get("marker_scale", 1.0)
    arrow_size = max(2, int(bar_width // 3 * marker_scale))
    pending_markers = []

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_x_end[current_sample] = current_x - read_spacing
            current_x += sample_spacing

        if sample not in sample_x_start:
            sample_x_start[sample] = current_x

        current_sample = sample
        read_positions.append((read_id, current_x))

        if feature_mode == "raw":
            for start, end, feature in features:
                y_start = top_margin + int(start * ratio)
                y_end = top_margin + int(end * ratio)
                height = max(1, y_end - y_start)
                color_hex = colors.get(feature, "#ffffff")
                color_rgba = hex_to_rgba(color_hex)

                x1 = max(0, current_x)
                x2 = min(image_width, current_x + bar_width)
                y1 = max(0, y_start)
                y2 = min(image_height, y_start + height)
                _draw_rect_rgba(img_array, x1, y1, x2, y2, color_rgba)
        else:
            colored_feats = _build_colored_features(features, colors, min_width_exclude)
            max_feat_end = max((e for _, e, _ in features), default=0)
            bar_length = int(max_feat_end * ratio)
            rasterized = rasterize_features(colored_feats, bar_length, ratio,
                                            feature_mode, oversample, min_feature_width)
            if rasterized:
                for run in rasterized:
                    y1 = top_margin + int(run['scaled_start'])
                    y2 = top_margin + int(run['scaled_stop'] + 0.5)
                    color_rgba = hex_to_rgba(run['color'])
                    if run['fill_opacity'] < 1.0:
                        color_rgba = (color_rgba[0], color_rgba[1], color_rgba[2],
                                      int(color_rgba[3] * run['fill_opacity']))
                    x1 = max(0, current_x)
                    x2 = min(image_width, current_x + bar_width)
                    y1 = max(0, y1)
                    y2 = min(image_height, y2)
                    _draw_rect_rgba(img_array, x1, y1, x2, y2, color_rgba)

        if read_border:
            bx1 = max(0, current_x)
            bx2 = min(image_width, current_x + bar_width)
            read_top = top_margin
            max_feat_end = max((e for _, e, _ in features), default=0)
            read_bottom = min(image_height, top_margin + int(max_feat_end * ratio))
            if bx2 > bx1 and read_bottom > read_top:
                border_color = np.array(text_color + (255,), dtype=np.uint8)
                bt = max(1, int(width_scale))
                img_array[read_top:read_top + bt, bx1:bx2] = border_color
                rb = min(read_bottom, image_height - 1)
                img_array[rb - bt + 1:rb + 1, bx1:bx2] = border_color
                img_array[read_top:read_bottom + 1, bx1:bx1 + bt] = border_color
                rx = min(bx2, image_width - 1)
                img_array[read_top:read_bottom + 1, rx - bt + 1:rx + 1] = border_color

        if read_id in markers:
            for m_start, m_end in markers[read_id]:
                mid_y = top_margin + int((m_start + m_end) / 2 * ratio)
                pending_markers.append((current_x, mid_y))

        current_x += bar_width + read_spacing

    if current_sample:
        sample_x_end[current_sample] = current_x - read_spacing

    # Draw heatmap grid if enabled
    if config.get("heatmap") and config.get("metadata_columns"):
        draw_heatmap_grid_png(img_array, read_positions, config, width_scale, height_scale)

    # Convert back to PIL Image for text rendering
    img = Image.fromarray(img_array)
    draw_ctx = ImageDraw.Draw(img)

    # Draw accumulated marker arrowheads
    if pending_markers:
        text_color_fill = text_color  # tuple (R, G, B)
        for mx, my, in pending_markers:
            draw_ctx.polygon(
                [(mx - arrow_size - 1, my - arrow_size),
                 (mx - 1, my),
                 (mx - arrow_size - 1, my + arrow_size)],
                fill=text_color_fill)

    # Draw inline scale bar (matching SVG draw_scale_bar at line 2786)
    if config.get("draw_scale_bar", True):
        scale_bar_bp = config.get("scale_bar_bp") or 10000
        scale_bar_height = int(scale_bar_bp * ratio)
        max_scale_height = max_height_px
        if not config.get("scale_bar_bp") and scale_bar_height > max_scale_height:
            scale_bar_bp = 5000
            scale_bar_height = int(scale_bar_bp * ratio)

        sb_bar_w = max(1, int(1.5 * height_scale))
        sb_tick_ext = max(2, int(3 * height_scale))
        sb_x = effective_left_margin - int(10 * width_scale)
        sb_y = top_margin

        # Vertical bar
        draw_ctx.rectangle([sb_x, sb_y, sb_x + sb_bar_w, sb_y + scale_bar_height],
                           fill=text_color)
        # Tick marks (symmetric around bar center)
        bar_center_x = sb_x + sb_bar_w // 2
        tick_w = max(1, int(1 * height_scale))
        draw_ctx.line([(bar_center_x - sb_tick_ext, sb_y),
                       (bar_center_x + sb_tick_ext, sb_y)],
                      fill=text_color, width=tick_w)
        draw_ctx.line([(bar_center_x - sb_tick_ext, sb_y + scale_bar_height),
                       (bar_center_x + sb_tick_ext, sb_y + scale_bar_height)],
                      fill=text_color, width=tick_w)

        # Rotated label
        sb_label = f"{scale_bar_bp // 1000} Kbp"
        sb_font = _load_font(max(font_size, int(config.get("font_size", 14) * height_scale)))
        bbox = sb_font.getbbox(sb_label)
        tw = bbox[2] - bbox[0]
        th = bbox[3] - bbox[1]
        pad = max(4, int(2 * height_scale))
        txt_img = Image.new('RGBA', (tw + 2 * pad, th + 2 * pad), (0, 0, 0, 0))
        txt_draw = ImageDraw.Draw(txt_img)
        txt_draw.text((pad, pad - bbox[1]), sb_label, fill=text_color + (255,), font=sb_font)
        txt_img = txt_img.rotate(90, expand=True)
        bar_center_x = sb_x + sb_bar_w // 2
        paste_x = max(0, bar_center_x - sb_tick_ext - txt_img.width - int(3 * width_scale))
        paste_y = sb_y + (scale_bar_height - txt_img.height) // 2
        img.paste(txt_img, (paste_x, paste_y), txt_img)

    # Draw heatmap row labels (matching SVG draw_heatmap_grid_svg lines 2206-2216)
    if config.get("heatmap") and config.get("metadata_columns"):
        heatmap_display_names = config.get("heatmap_display_names", {})
        metadata_columns = config["metadata_columns"]
        base_bar_width = config["bar_width"]
        box_h = max(1, int(base_bar_width * width_scale))
        scaled_bottom_gap = int(config["heatmap_bottom_gap"] * height_scale)
        scaled_row_gap = max(1, int(config["heatmap_row_gap"] * height_scale))
        for row_idx, col_name in enumerate(reversed(metadata_columns)):
            row_y = top_margin - scaled_bottom_gap - (row_idx + 1) * box_h - row_idx * scaled_row_gap
            display_name = heatmap_display_names.get(col_name, col_name)
            bbox = font.getbbox(display_name)
            text_w = bbox[2] - bbox[0]
            text_h = bbox[3] - bbox[1]
            label_y = row_y + box_h // 2 - text_h // 2 - bbox[1]
            draw_ctx.text(
                (left_margin - int(5 * width_scale) - text_w, label_y),
                display_name, fill=text_color, font=font)

    # Draw labels and lines (unless --no-header)
    if not config.get("no_header", False):
        label_interval = 300
        read_width = bar_width + read_spacing
        hm_offset = int(config.get("heatmap_total", 0) * height_scale)
        label_tiers = config.get("label_tiers", {})
        group_subgroup_order = config.get("group_subgroup_order", [])
        has_grouping = bool(label_tiers)

        if has_grouping:
            base_font_size = config.get("font_size", 11)
            tier_offset = int((base_font_size + 10) * height_scale)
            tier_display_names = config.get("tier_display_names", {})

            # Tier 2 (bottom): subgroup lines + labels
            # SVG: line at top_margin - 5 - hm_offset, label baseline at top_margin - 10 - hm_offset
            # PIL y = SVG baseline y - font_size (PIL draws from top-left, SVG from baseline)
            tier2_line_y = top_margin - int(5 * height_scale) - hm_offset
            tier2_label_y = top_margin - int(10 * height_scale) - hm_offset - font_size
            for sample in sample_order:
                if sample in sample_x_start and sample in sample_x_end:
                    draw_ctx.line(
                        [(sample_x_start[sample], tier2_line_y),
                         (sample_x_end[sample], tier2_line_y)],
                        fill=text_color, width=max(1, int(1.5 * height_scale))
                    )
                    _, subgroup = label_tiers.get(sample, (sample, None))
                    label_text = (subgroup or sample).replace("_", " ")
                    draw_ctx.text(
                        (sample_x_start[sample], tier2_label_y),
                        label_text, fill=text_color, font=font)

            # Tier 2 display name at left margin (e.g., "Condition")
            if 1 in tier_display_names:
                label_text = tier_display_names[1]
                bbox = font.getbbox(label_text)
                text_w = bbox[2] - bbox[0]
                draw_ctx.text(
                    (left_margin - int(5 * width_scale) - text_w, tier2_label_y),
                    label_text, fill=text_color, font=font)

            # Tier 1 (top): group lines + labels
            # SVG: line at top_margin - 5 - tier_offset - hm_offset, label baseline at tier1_line_y - 5
            group_spans = compute_group_spans(
                sample_order, sample_x_start, sample_x_end, group_subgroup_order)
            tier1_line_y = top_margin - int(5 * height_scale) - tier_offset - hm_offset
            tier1_label_y = tier1_line_y - int(5 * height_scale) - font_size
            for group_name, gx_start, gx_end in group_spans:
                draw_ctx.line(
                    [(gx_start, tier1_line_y), (gx_end, tier1_line_y)],
                    fill=text_color, width=max(1, int(1.5 * height_scale))
                )
                draw_ctx.text(
                    (gx_start, tier1_label_y),
                    group_name.replace("_", " "), fill=text_color, font=font)

            # Tier 1 display name at left margin (e.g., "Satellite")
            if 0 in tier_display_names:
                label_text = tier_display_names[0]
                bbox = font.getbbox(label_text)
                text_w = bbox[2] - bbox[0]
                draw_ctx.text(
                    (left_margin - int(5 * width_scale) - text_w, tier1_label_y),
                    label_text, fill=text_color, font=font)
        else:
            # Single-tier labels
            # SVG: line at top_margin - 5 - hm_offset, label baseline at top_margin - 12 - hm_offset
            line_y = top_margin - int(5 * height_scale) - hm_offset
            label_y = top_margin - int(12 * height_scale) - hm_offset - font_size
            for sample in sample_order:
                if sample in sample_x_start and sample in sample_x_end:
                    draw_ctx.line(
                        [(sample_x_start[sample], line_y),
                         (sample_x_end[sample], line_y)],
                        fill=text_color, width=max(1, int(1.5 * height_scale))
                    )
                    label_text = sample.replace("_", " ")
                    num_reads = sample_read_counts[sample]
                    for i in range(0, num_reads, label_interval):
                        x_pos = sample_x_start[sample] + (i * read_width)
                        draw_ctx.text((x_pos, label_y), label_text, fill=text_color, font=font)

    _save_png(img, output_path, bg_color)

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

    # Height scale: match width_scale for uniform scaling (PNG matches SVG proportions)
    desired_height_scale = width_scale
    max_height_scale = MAX_PNG_DIMENSION / base_height
    height_scale = min(desired_height_scale, max_height_scale)

    if width_scale < requested_scale:
        logger.info(f"  Note: Capping width scale from {requested_scale}x to {width_scale:.2f}x "
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

    logger.info(f"  PNG scale: height={height_scale:.2f}x, width={width_scale:.2f}x")
    logger.info(f"  Target dimensions: {image_width} x {image_height} pixels")

    # Create image
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    img = Image.new('RGBA', (image_width, image_height), bg_color)

    # Try to load font
    font_size = int(config.get("font_size", 11) * height_scale)
    font = _load_font(font_size)

    # Track positions
    current_y = top_margin
    current_sample = None
    sample_y_start = {}
    sample_y_end = {}

    # Draw reads using numpy for speed
    logger.info(f"  Drawing {total_reads} reads...")
    img_array = np.array(img)

    feature_mode = config.get("feature_mode", "raw")
    min_feature_width = config.get("min_feature_width", 0.5)
    min_width_exclude = config.get("min_width_exclude", [])
    oversample = config.get("oversample", 1)
    read_border = config.get("read_border", False)

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_y_end[current_sample] = current_y
            current_y += sample_spacing

        if sample not in sample_y_start:
            sample_y_start[sample] = current_y

        current_sample = sample

        if feature_mode == "raw":
            for start, end, feature in features:
                x_start = left_margin + int(start * ratio)
                x_end = left_margin + int(end * ratio)
                width = max(1, x_end - x_start)
                color_hex = colors.get(feature, "#ffffff")
                color_rgba = hex_to_rgba(color_hex)

                x1 = max(0, x_start)
                x2 = min(image_width, x_start + width)
                y1 = max(0, current_y)
                y2 = min(image_height, current_y + bar_width)
                _draw_rect_rgba(img_array, x1, y1, x2, y2, color_rgba)
        else:
            colored_feats = _build_colored_features(features, colors, min_width_exclude)
            max_feat_end = max((e for _, e, _ in features), default=0)
            bar_length = int(max_feat_end * ratio)
            rasterized = rasterize_features(colored_feats, bar_length, ratio,
                                            feature_mode, oversample, min_feature_width)
            if rasterized:
                for run in rasterized:
                    x1 = left_margin + int(run['scaled_start'])
                    x2 = left_margin + int(run['scaled_stop'] + 0.5)
                    color_rgba = hex_to_rgba(run['color'])
                    if run['fill_opacity'] < 1.0:
                        color_rgba = (color_rgba[0], color_rgba[1], color_rgba[2],
                                      int(color_rgba[3] * run['fill_opacity']))
                    x1 = max(0, x1)
                    x2 = min(image_width, x2)
                    y1 = max(0, current_y)
                    y2 = min(image_height, current_y + bar_width)
                    _draw_rect_rgba(img_array, x1, y1, x2, y2, color_rgba)

        if read_border:
            by1 = max(0, current_y)
            by2 = min(image_height, current_y + bar_width)
            max_feat_end = max((e for _, e, _ in features), default=0)
            read_right = min(image_width, left_margin + int(max_feat_end * ratio))
            read_left = left_margin
            if by2 > by1 and read_right > read_left:
                border_color = np.array(text_color + (255,), dtype=np.uint8)
                img_array[by1, read_left:read_right] = border_color
                img_array[min(by2, image_height - 1), read_left:read_right] = border_color
                img_array[by1:by2 + 1, read_left] = border_color
                img_array[by1:by2 + 1, min(read_right, image_width - 1)] = border_color

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
        font_size = config.get("font_size", 14)
        label_tiers = config.get("label_tiers", {})
        group_subgroup_order = config.get("group_subgroup_order", [])
        has_grouping = bool(label_tiers)

        if has_grouping:
            tier_offset = font_size + 10
            line_x = left_margin - int(5 * height_scale)

            # Tier 2 (inner): subgroup lines + labels
            for sample in sample_order:
                if sample in sample_y_start and sample in sample_y_end:
                    draw_ctx.line(
                        [(line_x, sample_y_start[sample]),
                         (line_x, sample_y_end[sample])],
                        fill=text_color, width=max(1, int(1.5 * height_scale))
                    )
                    _, subgroup = label_tiers.get(sample, (sample, None))
                    label_text = (subgroup or sample).replace("_", " ")
                    bbox = font.getbbox(label_text)
                    tw = bbox[2] - bbox[0]
                    th = bbox[3] - bbox[1]
                    txt_img = Image.new('RGBA', (tw + 4, th + 4), (0, 0, 0, 0))
                    txt_draw = ImageDraw.Draw(txt_img)
                    txt_draw.text((2, 2), label_text, fill=text_color + (255,), font=font)
                    txt_img = txt_img.rotate(90, expand=True)
                    paste_x = max(0, line_x - txt_img.width - int(5 * height_scale))
                    paste_y = sample_y_start[sample]
                    if paste_y + txt_img.height <= image_height:
                        img.paste(txt_img, (paste_x, paste_y), txt_img)

            # Tier 1 (outer): group lines + labels
            group_spans = compute_group_spans(
                sample_order, sample_y_start, sample_y_end, group_subgroup_order)
            tier1_x = line_x - tier_offset
            for group_name, gy_start, gy_end in group_spans:
                draw_ctx = ImageDraw.Draw(img)  # Refresh after paste
                draw_ctx.line(
                    [(tier1_x, gy_start), (tier1_x, gy_end)],
                    fill=text_color, width=max(1, int(1.5 * height_scale))
                )
                label_text = group_name.replace("_", " ")
                bbox = font.getbbox(label_text)
                tw = bbox[2] - bbox[0]
                th = bbox[3] - bbox[1]
                txt_img = Image.new('RGBA', (tw + 4, th + 4), (0, 0, 0, 0))
                txt_draw = ImageDraw.Draw(txt_img)
                txt_draw.text((2, 2), label_text, fill=text_color + (255,), font=font)
                txt_img = txt_img.rotate(90, expand=True)
                paste_x = max(0, tier1_x - txt_img.width - int(5 * height_scale))
                paste_y = gy_start
                if paste_y + txt_img.height <= image_height:
                    img.paste(txt_img, (paste_x, paste_y), txt_img)
        else:
            for sample in sample_order:
                if sample in sample_y_start and sample in sample_y_end:
                    # Vertical line at left spanning this sample
                    line_x = left_margin - int(5 * height_scale)
                    draw_ctx.line(
                        [(line_x, sample_y_start[sample]),
                         (line_x, sample_y_end[sample])],
                        fill=text_color, width=max(1, int(1.5 * height_scale))
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

    # Draw horizontal scale bar at top
    if config.get("draw_scale_bar", True):
        draw_ctx = ImageDraw.Draw(img)
        scale_bar_bp = config.get("scale_bar_bp") or 10000
        scale_bar_width_px = int(scale_bar_bp * ratio)
        if not config.get("scale_bar_bp") and scale_bar_width_px > max_width_px:
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

    _save_png(img, output_path, bg_color)

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
        scale_bar_bp = config.get("scale_bar_bp") or 10000
        scale_bar_width_px = int(scale_bar_bp * ratio)
        if not config.get("scale_bar_bp") and scale_bar_width_px > max_width_px:
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
        d.append(draw.Text(label, config.get("font_size", 14),
                           sb_x + scale_bar_width_px / 2, sb_y - 8,
                           fill=text_color, text_anchor="middle",
                           font_family="Basic Sans"))

    current_y = top_margin
    current_sample = None
    sample_y_start = {}
    sample_y_end = {}

    feature_mode = config.get("feature_mode", "raw")
    min_feature_width = config.get("min_feature_width", 0.5)
    min_width_exclude = config.get("min_width_exclude", [])
    oversample = config.get("oversample", 1)
    read_border = config.get("read_border", False)

    for sample, read_id, read_length, features in reads:
        if current_sample is not None and sample != current_sample:
            sample_y_end[current_sample] = current_y
            current_y += sample_spacing

        if sample not in sample_y_start:
            sample_y_start[sample] = current_y

        current_sample = sample

        if feature_mode == "raw":
            for start, end, feature in features:
                x_start = left_margin + int(start * ratio)
                width = max(1, int((end - start) * ratio))
                color = colors.get(feature, "#ffffff")
                d.append(draw.Rectangle(x_start, current_y, width, bar_width, fill=color))
        else:
            colored_feats = _build_colored_features(features, colors, min_width_exclude)
            max_feat_end = max((e for _, e, _ in features), default=0)
            bar_length = int(max_feat_end * ratio)
            rasterized = rasterize_features(colored_feats, bar_length, ratio,
                                            feature_mode, oversample, min_feature_width)
            if rasterized:
                for run in rasterized:
                    rx = left_margin + run['scaled_start']
                    rw = run['scaled_stop'] - run['scaled_start']
                    d.append(draw.Rectangle(rx, current_y, rw, bar_width,
                                            fill=run['color'],
                                            fill_opacity=run['fill_opacity']))

        if read_border:
            max_feat_end = max((e for _, e, _ in features), default=0)
            read_width_px = int(max_feat_end * ratio)
            d.append(draw.Rectangle(left_margin, current_y, read_width_px, bar_width,
                                    fill='none', stroke=text_color, stroke_width=0.5))

        current_y += bar_width + read_spacing

    if current_sample:
        sample_y_end[current_sample] = current_y

    # Draw labels and separators (unless --no-header)
    if not config.get("no_header", False):
        font_family = "Basic Sans"
        label_interval = 300
        read_height = bar_width + read_spacing
        font_size = config.get("font_size", 14)
        label_tiers = config.get("label_tiers", {})
        group_subgroup_order = config.get("group_subgroup_order", [])
        has_grouping = bool(label_tiers)

        if has_grouping:
            tier_offset = font_size + 10
            line_x = left_margin - 5
            tier_display_names = config.get("tier_display_names", {})

            # Tier 2 (inner): subgroup lines + labels
            for sample in sample_order:
                if sample in sample_y_start and sample in sample_y_end:
                    d.append(draw.Line(line_x, sample_y_start[sample],
                                       line_x, sample_y_end[sample],
                                       stroke=text_color, stroke_width=1.5))
                    _, subgroup = label_tiers.get(sample, (sample, None))
                    label_text = (subgroup or sample).replace("_", " ")
                    ly = sample_y_start[sample] + 50
                    lx = line_x - 10
                    d.append(draw.Text(
                        label_text, font_size, lx, ly,
                        fill=text_color, text_anchor="middle",
                        font_family=font_family,
                        transform=f"rotate(-90, {lx}, {ly})",
                    ))

            # Tier 2 display name (top of margin)
            if 1 in tier_display_names:
                d.append(draw.Text(
                    tier_display_names[1], font_size,
                    line_x, top_margin - 10,
                    fill=text_color, text_anchor="middle",
                    font_family=font_family,
                ))

            # Tier 1 (outer): group lines + labels
            group_spans = compute_group_spans(
                sample_order, sample_y_start, sample_y_end, group_subgroup_order)
            tier1_x = line_x - tier_offset
            for group_name, gy_start, gy_end in group_spans:
                d.append(draw.Line(tier1_x, gy_start, tier1_x, gy_end,
                                   stroke=text_color, stroke_width=1.5))
                ly = gy_start + 50
                d.append(draw.Text(
                    group_name.replace("_", " "), font_size, tier1_x - 10, ly,
                    fill=text_color, text_anchor="middle",
                    font_family=font_family,
                    transform=f"rotate(-90, {tier1_x - 10}, {ly})",
                ))

            # Tier 1 display name (top of margin)
            if 0 in tier_display_names:
                d.append(draw.Text(
                    tier_display_names[0], font_size,
                    tier1_x, top_margin - 10,
                    fill=text_color, text_anchor="middle",
                    font_family=font_family,
                ))
        else:
            for sample in sample_order:
                if sample in sample_y_start and sample in sample_y_end:
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
                            label_text, font_size, lx, ly,
                            fill=text_color, text_anchor="middle",
                            font_family=font_family,
                            transform=f"rotate(-90, {lx}, {ly})",
                        ))

    d.save_svg(output_path)
    logger.info("\nSaved: %s", output_path)
    logger.info(f"  Dimensions: {image_width} x {image_height} pixels")

    return image_width


def draw_legend_png(features_used, colors, config, horizontal=False, scale=1.0,
                    max_width=None, extra_items=None):
    """Render an auto-filtered color legend as a PIL Image.

    Uses compute_legend_layout() for grid-aligned columns with section headers.
    Optional extra_items are prepended (e.g. heatmap legend entries).

    Returns:
        PIL.Image: legend image ready for compositing
    """
    background = config["background"]
    bg_color = (0, 0, 0, 255) if background == "black" else (255, 255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    border_color = (255, 255, 255, 255) if background == "black" else (0, 0, 0, 255)
    color_sections = config.get("color_sections")

    base_font = config.get("font_size", 14)
    swatch_size = int(base_font * scale)
    font_size = int(base_font * scale)
    item_padding = int(3 * scale)

    font = _load_font(font_size)

    def measure_text(text, _fs):
        bbox = font.getbbox(text)
        return bbox[2] - bbox[0]

    filtered = _filter_legend_features(features_used, colors, color_sections)
    if extra_items:
        filtered = extra_items + filtered
    if not filtered:
        return None

    max_cols = 1 if horizontal else 6
    layout = compute_legend_layout(
        filtered, max_cols=max_cols, max_width=max_width,
        swatch_size=swatch_size, font_size=font_size,
        padding=int(10 * scale), col_gap=int(20 * scale),
        item_padding=item_padding, measure_text_fn=measure_text)

    legend_img = Image.new('RGBA', (layout.width, layout.height), bg_color)
    ctx = ImageDraw.Draw(legend_img)

    for item in layout.items:
        if item.is_header:
            ctx.text((int(item.x), int(item.y + 2)), item.label,
                     fill=text_color, font=font)
        else:
            # Color swatch with border
            sx = int(item.x) + layout.item_padding
            sy = int(item.y) + layout.item_padding
            rgba = hex_to_rgba(item.color)
            ctx.rectangle([sx, sy, sx + swatch_size, sy + swatch_size],
                          fill=rgba, outline=border_color, width=max(1, int(scale)))
            # Label text
            swatch_gap = int(swatch_size * 0.4)
            tx = sx + swatch_size + swatch_gap
            ctx.text((tx, sy), item.label, fill=text_color, font=font)

    return legend_img


LegendItem = namedtuple('LegendItem', ['x', 'y', 'label', 'color', 'is_header'])
LegendLayout = namedtuple('LegendLayout', ['items', 'width', 'height', 'col_width',
                                            'num_cols', 'swatch_size', 'font_size',
                                            'row_height', 'item_padding'])


def _build_heatmap_legend_items(config):
    """Convert heatmap color data to legend item tuples.

    Returns list of (display_name, color_hex, section_header) matching
    the format used by _filter_legend_features / compute_legend_layout.
    """
    items = []
    heatmap_colors = config.get("heatmap_colors", {})
    heatmap_display_names = config.get("heatmap_display_names", {})
    for col_name in config.get("metadata_columns", []):
        val_colors = heatmap_colors.get(col_name, {})
        header = heatmap_display_names.get(col_name, col_name)
        for val, color in val_colors.items():
            if val is None:
                continue
            display = str(val).replace("_", " ")
            items.append((display, color, header))
    return items


def _filter_legend_features(features_used, colors, color_sections=None):
    """Filter features for legend display.

    Deduplicates _specific variants, trims _multigroup suffixes,
    replaces underscores with spaces, and preserves section grouping.

    Returns:
        list of (display_name, color_hex, section_header_or_None) tuples
    """
    # Build set of base feature names actually in use
    used_bases = set()
    for feat in features_used:
        base = re.sub(r'_specific$', '', feat)
        used_bases.add(base)

    seen_display = set()
    result = []

    if color_sections:
        for header, section_feats in color_sections:
            for feat in section_feats:
                base = re.sub(r'_specific$', '', feat)
                if base not in used_bases:
                    continue
                display = re.sub(r'_multigroup\d+$', '', base).replace('_', ' ')
                if display in seen_display:
                    continue
                seen_display.add(display)
                color_hex = colors.get(feat, colors.get(base, None))
                if color_hex:
                    result.append((display, color_hex, header))
    else:
        for feat in sorted(features_used):
            base = re.sub(r'_specific$', '', feat)
            display = re.sub(r'_multigroup\d+$', '', base).replace('_', ' ')
            if display in seen_display:
                continue
            seen_display.add(display)
            color_hex = colors.get(feat, colors.get(base, None))
            if color_hex:
                result.append((display, color_hex, None))

    return result


def compute_legend_layout(filtered_items, *, max_width=None, max_cols=6,
                          swatch_size=12, font_size=12, padding=10,
                          col_gap=20, item_padding=3, measure_text_fn=None):
    """Compute a grid layout for legend items with optional section headers.

    Items fill top-to-bottom, then left-to-right. Section headers start a new
    column if the current column is non-empty.

    Returns:
        LegendLayout with positioned items and total dimensions.
    """
    if measure_text_fn is None:
        measure_text_fn = lambda text, fs: len(text) * fs * 0.48

    if not filtered_items:
        return LegendLayout([], 0, 0, 0, 0, swatch_size, font_size, 0, item_padding)

    swatch_gap = int(swatch_size * 0.4)
    row_height = swatch_size + 2 * item_padding + 4

    # Measure widest label to get uniform column width
    max_label_w = 0
    for display, _, _ in filtered_items:
        max_label_w = max(max_label_w, measure_text_fn(display, font_size))
    col_width = swatch_size + swatch_gap + int(max_label_w) + 2 * item_padding + 4

    # Determine number of columns
    if max_width:
        available = max_width - 2 * padding
        num_cols = min(max_cols, max(1, int((available + col_gap) / (col_width + col_gap))))
    else:
        num_cols = min(max_cols, max(1, len(filtered_items)))

    # Group items by section
    sections = []
    current_header = None
    current_items = []
    first = True
    for display, color_hex, header in filtered_items:
        if header != current_header and not first:
            sections.append((current_header, current_items))
            current_items = []
        current_header = header
        current_items.append((display, color_hex))
        first = False
    if current_items:
        sections.append((current_header, current_items))

    # Count total slots (header rows + item rows)
    total_slots = 0
    for header, items in sections:
        if header:
            total_slots += 1  # header row
        total_slots += len(items)

    num_rows = max(1, -(-total_slots // num_cols))  # ceil division

    # Place items column-first
    positioned = []
    col = 0
    row = 0
    for header, items in sections:
        # Start new column for new section if current column has items
        if header and row > 0:
            col += 1
            row = 0
            if col >= num_cols:
                col = 0
                # All columns full, just continue stacking

        if header:
            x = padding + col * (col_width + col_gap)
            y = padding + row * row_height
            positioned.append(LegendItem(x, y, header, None, True))
            row += 1

        for display, color_hex in items:
            if row >= num_rows and col + 1 < num_cols:
                col += 1
                row = 0
            x = padding + col * (col_width + col_gap)
            y = padding + row * row_height
            positioned.append(LegendItem(x, y, display, color_hex, False))
            row += 1

    # Compute actual dimensions from placed items
    max_x = max(item.x for item in positioned) + col_width if positioned else 0
    max_y = max(item.y for item in positioned) + row_height if positioned else 0
    total_width = max_x + padding
    total_height = max_y + padding

    actual_cols = (max(item.x for item in positioned) - padding) // (col_width + col_gap) + 1 if positioned else 0

    return LegendLayout(positioned, int(total_width), int(total_height),
                        col_width, int(actual_cols), swatch_size, font_size,
                        row_height, item_padding)


def _estimate_heatmap_legend_height(config, image_width, left_margin):
    """Pre-calculate height needed for the heatmap legend with row-wrapping."""
    metadata_columns = config.get("metadata_columns", [])
    heatmap_colors = config.get("heatmap_colors", {})
    heatmap_display_names = config.get("heatmap_display_names", {})
    font_size = config.get("font_size", 14)
    swatch_size = font_size
    col_gap = 20
    row_height = swatch_size + 10
    right_margin = 50
    max_x = image_width - right_margin

    x = left_margin
    num_rows = 1
    for col_name in metadata_columns:
        val_colors = heatmap_colors.get(col_name, {})
        display_name = heatmap_display_names.get(col_name, col_name)
        header_w = len(display_name) * (font_size * 0.6) + 10
        if x + header_w > max_x and x > left_margin:
            num_rows += 1
            x = left_margin
        x += header_w
        for val, color in val_colors.items():
            if val is None:
                continue
            item_w = swatch_size + len(str(val)) * (font_size * 0.6) + 15
            if x + item_w > max_x and x > left_margin:
                num_rows += 1
                x = left_margin
            x += item_w
        x += col_gap

    return num_rows * row_height


def estimate_featureset_legend_height(features_used, colors, config, image_width,
                                      extra_items=None):
    """Pre-calculate the height the SVG featureset legend will consume."""
    color_sections = config.get("color_sections")
    filtered = _filter_legend_features(features_used, colors, color_sections)
    if extra_items:
        filtered = extra_items + filtered
    if not filtered:
        return 0
    font_size = config.get("font_size", 14)
    layout = compute_legend_layout(
        filtered, max_width=image_width, swatch_size=font_size, font_size=font_size)
    return layout.height


def draw_legend_svg(d, features_used, colors, config, legend_y, image_width=None,
                    extra_items=None):
    """Draw color legend in SVG below reads.

    Uses compute_legend_layout() for grid-aligned columns with section headers.
    Optional extra_items are prepended (e.g. heatmap legend entries).

    Returns:
        int: total height consumed by the legend
    """
    background = config["background"]
    text_color = "#ffffff" if background == "black" else "#000000"
    border_color = "#ffffff" if background == "black" else "#000000"
    font_size = config.get("font_size", 14)
    font_family = "Basic Sans"
    color_sections = config.get("color_sections")

    filtered = _filter_legend_features(features_used, colors, color_sections)
    if extra_items:
        filtered = extra_items + filtered
    if not filtered:
        return 0

    layout = compute_legend_layout(
        filtered, max_width=image_width, swatch_size=font_size, font_size=font_size)

    for item in layout.items:
        ix = item.x
        iy = legend_y + item.y

        if item.is_header:
            d.append(draw.Text(
                item.label, font_size, ix, iy + font_size,
                fill=text_color, text_anchor="start", font_family=font_family,
                font_weight="bold",
            ))
        else:
            # Color swatch with border
            swatch_y = iy + layout.item_padding
            swatch_x = ix + layout.item_padding
            d.append(draw.Rectangle(
                swatch_x, swatch_y, layout.swatch_size, layout.swatch_size,
                fill=item.color, stroke=border_color, stroke_width=0.5,
            ))
            # Label text
            swatch_gap = int(layout.swatch_size * 0.4)
            text_x = swatch_x + layout.swatch_size + swatch_gap
            text_y = swatch_y + layout.swatch_size - 1
            d.append(draw.Text(
                item.label, font_size, text_x, text_y,
                fill=text_color, text_anchor="start", font_family=font_family,
            ))

    return layout.height


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
    logger.info(f"  Legend composited ({position}): {png_path}")


def load_cluster_reads(cluster_prefix, cluster_ids, results_dir,
                       database, featureset, smoothness, analysis):
    """Load all reads from specified clusters as plot_reads-format tuples.

    Returns:
        tuple: (reads_list, sample_order)
    """
    import pandas as pd

    assignments_file = f"{cluster_prefix}.sequence_assignments.tsv"
    if not os.path.exists(assignments_file):
        logger.error(f"Error: Cluster assignments file not found: {assignments_file}")
        return [], []

    logger.info(f"  Loading cluster assignments: {assignments_file}")
    df = pd.read_csv(assignments_file, sep='\t')
    if 'read' in df.columns:
        df = df.rename(columns={'read': 'sequence'})

    df = df[df['cluster'].isin(cluster_ids)]
    logger.info(f"  Found {len(df)} reads in clusters: {cluster_ids}")

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
        logger.info(f"  {label}: {found} reads with BED data")

    return result, sample_order


def draw_scale_bar_png(image_height, top_margin, ratio, background, output_path, scale_bar_bp_override=None,
                       height_scale=1.0, font_size=14):
    """Render a standalone vertical scale bar as PNG for animation overlay.

    Generated at source PNG dimensions with scaled ratio. The animation function
    will crop+resize this to match the viewport, keeping the bar proportional
    to the visible reads.
    """
    width = 60
    bg_color = (0, 0, 0) if background == "black" else (255, 255, 255)
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)

    scale_bar_bp = scale_bar_bp_override or 10000
    scale_bar_height = int(scale_bar_bp * ratio)
    max_height_px = image_height - top_margin - 40
    if not scale_bar_bp_override and scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    img = Image.new('RGB', (width, image_height), bg_color)
    draw_ctx = ImageDraw.Draw(img)

    # Load font — scale with height_scale so text survives downscale to viewport
    scaled_font_size = max(font_size, int(font_size * height_scale))
    font = _load_font(scaled_font_size)

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
    logger.info(f"  Scale bar saved for animation: {output_path}")
    logger.info(f"  Dimensions: {width} x {image_height} pixels")


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

    logger.info(f"\nGenerating animation: {mp4_path}")
    logger.info(f"  Direction: {direction}, Duration: {duration:.1f}s, FPS: {args.animate_fps}")
    logger.info(f"  Viewport: {vw}x{vh}, Zoom: {args.animate_zoom}")
    if legend:
        logger.info(f"  Legend: {legend}")
    if scale_bar:
        logger.info(f"  Scale bar: {scale_bar}")

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

    logger.info(f"  Animation: {mp4_path}")


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
        "--output", "-o", required=True, help="Output file path (extension auto-adjusted based on --format)"
    )
    parser.add_argument(
        "--scale-bar-output",
        default=None,
        help="Output path for separate scale bar SVG (PNG also rendered if rsvg-convert is available)",
    )
    parser.add_argument(
        "--read-list",
        default=None,
        help="File with read IDs to include (one per line, or TSV with IDs in first column). "
             "Optionally, 2nd and 3rd columns provide group and subgroup labels for "
             "two-level grouping (e.g., satellite type and cell line).",
    )
    parser.add_argument(
        "--filter-group", action="append", default=None,
        help="Filter --read-list to rows where the first label-tier column matches "
             "this value. Can be specified multiple times to include several groups.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Minimum read length in bp (shorter reads discarded)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="Maximum read length in bp (longer reads discarded)",
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
        default=5,
        help="Width of each read bar in pixels (default: 5)",
    )
    parser.add_argument(
        "--read-spacing",
        type=int,
        default=5,
        help="Horizontal spacing between reads (default: 5)",
    )
    parser.add_argument(
        "--sample-spacing",
        type=int,
        default=20,
        help="Horizontal spacing between sample groups (default: 20)",
    )
    parser.add_argument(
        "--subgroup-spacing",
        type=int,
        default=None,
        help="Horizontal spacing between subgroups within the same group "
             "(default: same as --sample-spacing)",
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
        choices=["white", "black", "both"],
        help="Background color (default: black). Use 'both' to generate black and white variants.",
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
        "--orient-satellite-top",
        action="store_true",
        help="Reorient reads so satellite-dense end is at top "
             "(telomere features first, satellite density as fallback)",
    )
    parser.add_argument(
        "--markers",
        help="TSV file with read_id, marker_start, marker_end columns. "
             "Draws arrowheads on the left side of each read at marker positions.",
    )
    parser.add_argument(
        "--marker-scale",
        type=float,
        default=1.0,
        help="Scale factor for marker arrowhead size (default: 1.0)",
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
        "--scale-bar-bp",
        type=int,
        default=None,
        help="Override scale bar size in bp (default: auto 10kb/5kb)",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        default=11,
        help="Font size for labels and scale bar text (default: 11)",
    )
    parser.add_argument(
        "--viewport-ratio",
        default="16:9",
        help="Aspect ratio W:H for output (default: 16:9). "
             "Width is calculated from content, height derived from ratio. "
             "Set to 'none' to disable and use raw content dimensions.",
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

    # ── Heatmap ──
    parser.add_argument(
        "--heatmap",
        action="store_true",
        help="Draw a categorical heatmap grid above read bars using extra columns "
             "(4th+) from --read-list TSV. Requires a header row.",
    )
    parser.add_argument(
        "--label-tier", action="append", default=None,
        help="Column from --read-list TSV to use as a tiered label, "
             "with optional display name: COLUMN:DISPLAY_NAME. "
             "Can be specified multiple times for multiple tiers (top to bottom order). "
             "Default: columns 2 and 3 (group, subgroup).",
    )
    parser.add_argument(
        "--heatmap-track", action="append", default=None,
        help="Column from --read-list TSV to use as a heatmap row, "
             "with optional display name: COLUMN:DISPLAY_NAME. "
             "Can be specified multiple times. "
             "Default with --heatmap: columns 4+ from header.",
    )
    parser.add_argument(
        "--max-read-length",
        type=int, default=None,
        help="Exclude reads longer than this many bp (default: no limit).",
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

    # ── Feature rendering ──
    parser.add_argument(
        "--feature-mode",
        default="smooth",
        choices=["raw", "transition", "smooth"],
        help="Feature rendering mode: smooth (default, windowed majority-vote "
             "downsampling), transition (direct scaling with min-width enforcement), "
             "raw (integer pixel stacking).",
    )
    parser.add_argument(
        "--min-feature-width",
        type=float,
        default=0.5,
        help="Min pixel width per feature in transition mode (default: 0.5).",
    )
    parser.add_argument(
        "--min-width-exclude",
        nargs="*",
        default=["novel", "*arm*", "ct*"],
        help="Glob patterns for features exempt from min-width inflation "
             "(default: novel *arm* ct*).",
    )
    parser.add_argument(
        "--oversample",
        type=int,
        default=1,
        help="Oversampling factor for smooth mode (default: 1).",
    )
    parser.add_argument(
        "--read-border",
        action="store_true",
        help="Draw thin black border around each read bar.",
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
    Parses '# comment' lines as section headers for legend grouping.

    Returns:
        tuple: (colors_dict, color_sections)
            - colors_dict: {feature_name: hex_color}
            - color_sections: [(section_header_or_None, [feature_names])]
    """
    colors = {"novel": "#ffffff"}  # Default for unknown features
    sections = []
    current_header = None
    current_features = []

    with open(colors_file, "r") as f:
        for line in f:
            stripped = line.strip()
            # Parse section headers from comment lines
            m = re.match(r'^#\s+(.+)', stripped)
            if m:
                # Save previous section if it has features
                if current_features:
                    sections.append((current_header, current_features))
                current_header = m.group(1).strip()
                current_features = []
                continue
            parts = stripped.split()
            if len(parts) < 2 or parts[0].lower() == "feature":
                continue
            feature, color = parts[0], parts[1]
            colors[feature] = color
            current_features.append(feature)
            # Also map without _specific suffix for smoothed BED files
            if feature.endswith("_specific"):
                colors[feature[:-9]] = color
            # Also map with _specific suffix
            if not feature.endswith("_specific") and not feature.endswith("_multigroup1"):
                colors[feature + "_specific"] = color

    # Save final section
    if current_features:
        sections.append((current_header, current_features))
    # If no sections found at all, wrap everything in a None-header section
    if not sections:
        sections = [(None, list(colors.keys()))]

    return colors, sections


def load_read_list(path):
    """Load read IDs and all column data from a TSV read-list file.

    Accepts:
      - 1 column: read IDs only (backward compatible)
      - 2+ columns: read_id + any number of data columns (TSV)

    Skips header lines (sequence/read/read_id) but captures column names.

    Returns:
        tuple: (read_ids, all_columns, read_data, read_rows)
        - read_ids: set of read ID strings
        - all_columns: ordered list of column names from the header (columns 2+)
        - read_data: dict mapping read_id -> {col_name: value} for all columns
          (last occurrence wins for duplicate read IDs)
        - read_rows: list of (read_id, {col_name: value}) in file order,
          preserving duplicates for correct group ordering
    """
    read_ids = set()
    all_columns = []
    read_data = {}
    read_rows = []

    open_func = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with open_func(path, mode) as f:
        for line in f:
            fields = line.strip().split("\t")
            if not fields or not fields[0]:
                continue
            rid = fields[0]
            if rid.lower() in ("sequence", "read", "read_id"):
                # Capture all column names from header (columns 2+)
                if len(fields) > 1:
                    all_columns = fields[1:]
                continue
            read_ids.add(rid)

            # Store all column values for this read
            if all_columns and len(fields) > 1:
                meta = {}
                for i, col_name in enumerate(all_columns):
                    idx = 1 + i
                    if idx < len(fields) and fields[idx]:
                        meta[col_name] = fields[idx]
                    else:
                        meta[col_name] = None
                read_data[rid] = meta
                read_rows.append((rid, meta))

    return read_ids, all_columns, read_data, read_rows


def apply_read_list_grouping(reads, read_groups, group_subgroup_order):
    """Replace sample field with composite group labels for two-level grouping.

    Args:
        reads: list of (sample, read_id, read_length, features) tuples
        read_groups: dict mapping read_id -> (group, subgroup)
        group_subgroup_order: ordered list of (group, subgroup) tuples

    Returns:
        tuple: (reads, sample_order, group_boundaries)
        - reads: updated read tuples with composite sample labels
        - sample_order: list of composite label strings in group order
        - group_boundaries: set of sample labels that are the LAST subgroup
          before a group change (used for thicker separators)
    """
    # Build composite labels
    def _label(group, subgroup):
        if subgroup:
            return f"{group} \u2014 {subgroup}"
        return group

    sample_order = [_label(g, s) for g, s in group_subgroup_order]

    # Relabel reads
    updated = []
    for sample, read_id, read_length, features in reads:
        if read_id in read_groups:
            group, subgroup = read_groups[read_id]
            new_sample = _label(group, subgroup)
            updated.append((new_sample, read_id, read_length, features))
        else:
            updated.append((sample, read_id, read_length, features))

    # Compute group boundaries: the LAST subgroup label before a group change
    group_boundaries = set()
    if len(group_subgroup_order) > 1:
        for i in range(len(group_subgroup_order) - 1):
            curr_group = group_subgroup_order[i][0]
            next_group = group_subgroup_order[i + 1][0]
            if curr_group != next_group:
                group_boundaries.add(sample_order[i])

    return updated, sample_order, group_boundaries


def compute_group_spans(sample_order, sample_starts, sample_ends, group_subgroup_order):
    """Compute positional spans for top-level groups from subgroup spans.

    Args:
        sample_order: list of composite label strings (e.g., "aSat — NHA p1")
        sample_starts: dict mapping composite label -> start position
        sample_ends: dict mapping composite label -> end position
        group_subgroup_order: list of (group, subgroup) tuples

    Returns:
        list of (group_label, span_start, span_end) in order.
        Only includes groups with at least one visible subgroup.
    """
    def _label(group, subgroup):
        if subgroup:
            return f"{group} \u2014 {subgroup}"
        return group

    result = []
    current_group = None
    group_start = None
    group_end = None

    for group, subgroup in group_subgroup_order:
        composite = _label(group, subgroup)
        if composite not in sample_starts:
            continue

        if group != current_group:
            if current_group is not None and group_start is not None:
                result.append((current_group, group_start, group_end))
            current_group = group
            group_start = sample_starts[composite]
            group_end = sample_ends.get(composite, group_start)
        else:
            group_end = sample_ends.get(composite, group_end)

    if current_group is not None and group_start is not None:
        result.append((current_group, group_start, group_end))

    return result


def _parse_bed_file(bed_path, sample):
    """Parse a BED file and return read tuples for a given sample.

    Returns:
        list of (sample, read_id, read_length, features) tuples
    """
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

    reads = []
    for read_id, features in read_features.items():
        read_length = max(end for _, end, _ in features)
        reads.append((sample, read_id, read_length, features))

    logger.info(f"  {sample}: {len(read_features)} reads loaded")
    return reads


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
                logger.warning(f"Warning: BED file not found for {sample}: {bed_path}")
                continue

        all_reads.extend(_parse_bed_file(bed_path, sample))

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
            logger.warning(f"Warning: BED file not found: {bed_path}")
            continue

        all_reads.extend(_parse_bed_file(bed_path, sample))

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


def compute_viewport_params(target_w, target_h, reads, config, horizontal=False):
    """Back-calculate bar_width, read_spacing, and ratio to fit a target viewport.

    For vertical mode: width governs read pitch, height governs ratio.
    For horizontal mode: width governs ratio, height governs read pitch.

    Returns:
        dict with updated bar_width, read_spacing, ratio values
    """
    sample_read_counts = defaultdict(int)
    for sample, _, _, _ in reads:
        sample_read_counts[sample] += 1
    num_samples = len(sample_read_counts)
    total_reads = len(reads)
    max_length = max(r[2] for r in reads)

    sample_spacing = config["sample_spacing"]
    top_margin = config["top_margin"]
    left_margin = config["left_margin"]
    bottom_margin = config["bottom_margin"]

    updates = {}

    if horizontal:
        # Width governs ratio (feature detail along x-axis)
        usable_w = target_w - left_margin - 50
        updates["ratio"] = usable_w / max_length if max_length > 0 else config["ratio"]

        # Height governs read pitch
        usable_h = target_h - top_margin - 50 - ((num_samples - 1) * sample_spacing)
        if total_reads > 0:
            pitch = max(2, usable_h / total_reads)
            updates["bar_width"] = max(1, int(pitch * 0.5))
            updates["read_spacing"] = max(1, int(pitch - updates["bar_width"]))
    else:
        # Width governs read pitch
        usable_w = target_w - left_margin - 50 - ((num_samples - 1) * sample_spacing)
        if total_reads > 0:
            pitch = max(2, usable_w / total_reads)
            updates["bar_width"] = max(1, int(pitch * 0.5))
            updates["read_spacing"] = max(1, int(pitch - updates["bar_width"]))

        # Height governs ratio
        usable_h = target_h - top_margin - bottom_margin
        updates["ratio"] = usable_h / max_length if max_length > 0 else config["ratio"]

    return updates


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

# Satellite features used for orientation fallback
SATELLITE_FEATURES = {
    'active', 'inactive', 'divergent', 'monomeric',
    'hsat1A', 'hsat1B', 'hsat2', 'hsat3',
    'bsat', 'gsat', 'censat',
    'noncentromeric',
}


def _orient_reads(reads, target_features, label):
    """Reorient reads so that target features are at the top (position 0).

    For each read, checks if any target features are closer to start or end.
    If closer to end, flips the read coordinates.

    Args:
        reads: List of (sample, read_id, read_length, features) tuples
        target_features: set of feature names to orient toward the top
        label: description for the print message (e.g. "telomere", "chromosome")

    Returns:
        Tuple of (list of reoriented reads, dict of {read_id: read_length} for flipped reads)
    """
    oriented_reads = []
    flipped_count = 0
    flipped_reads = {}

    for sample, read_id, read_length, features in reads:
        positions = []
        for start, end, feature in features:
            if feature in target_features:
                positions.extend([start, end])

        if positions:
            avg_pos = sum(positions) / len(positions)
            midpoint = read_length / 2

            if avg_pos > midpoint:
                flipped_features = [
                    (read_length - end, read_length - start, feature)
                    for start, end, feature in features
                ]
                flipped_features.sort(key=lambda x: x[0])
                oriented_reads.append((sample, read_id, read_length, flipped_features))
                flipped_count += 1
                flipped_reads[read_id] = read_length
            else:
                oriented_reads.append((sample, read_id, read_length, features))
        else:
            oriented_reads.append((sample, read_id, read_length, features))

    logger.info(f"  Reoriented {flipped_count} of {len(reads)} reads ({label} now at top)")
    return oriented_reads, flipped_reads


def orient_telomere_top(reads):
    """Reorient reads so telomere features are at the top (position 0)."""
    return _orient_reads(reads, TELOMERE_FEATURES, "telomere")


def orient_chromosome_top(reads):
    """Reorient reads so chromosome features are at the top (position 0)."""
    return _orient_reads(reads, CHROMOSOME_FEATURES, "chromosome")


def _orient_reads_with_fallback(reads, primary_features, fallback_features, label):
    """Reorient reads using primary features first, falling back to secondary.

    For each read:
    - If primary features found, orient by those (e.g. telomere at top).
    - Else if fallback features found, orient so fallback-dense end is at top.
    - Else leave unoriented.

    Args:
        reads: List of (sample, read_id, read_length, features) tuples
        primary_features: set of feature names to orient by first
        fallback_features: set of feature names to use as fallback
        label: description for the log message

    Returns:
        Tuple of (list of reoriented reads, dict of {read_id: read_length} for flipped reads)
    """
    oriented_reads = []
    flipped_count = 0
    primary_count = 0
    fallback_count = 0
    unoriented_count = 0
    flipped_reads = {}

    for sample, read_id, read_length, features in reads:
        # Try primary features first
        primary_positions = []
        for start, end, feature in features:
            if feature in primary_features:
                primary_positions.extend([start, end])

        if primary_positions:
            avg_pos = sum(primary_positions) / len(primary_positions)
            midpoint = read_length / 2
            if avg_pos > midpoint:
                flipped_features = [
                    (read_length - end, read_length - start, feature)
                    for start, end, feature in features
                ]
                flipped_features.sort(key=lambda x: x[0])
                oriented_reads.append((sample, read_id, read_length, flipped_features))
                flipped_count += 1
                flipped_reads[read_id] = read_length
            else:
                oriented_reads.append((sample, read_id, read_length, features))
            primary_count += 1
            continue

        # Try fallback features
        fallback_positions = []
        for start, end, feature in features:
            if feature in fallback_features:
                fallback_positions.extend([start, end])

        if fallback_positions:
            avg_pos = sum(fallback_positions) / len(fallback_positions)
            midpoint = read_length / 2
            if avg_pos > midpoint:
                flipped_features = [
                    (read_length - end, read_length - start, feature)
                    for start, end, feature in features
                ]
                flipped_features.sort(key=lambda x: x[0])
                oriented_reads.append((sample, read_id, read_length, flipped_features))
                flipped_count += 1
                flipped_reads[read_id] = read_length
            else:
                oriented_reads.append((sample, read_id, read_length, features))
            fallback_count += 1
            continue

        # No orientation features found
        oriented_reads.append((sample, read_id, read_length, features))
        unoriented_count += 1

    logger.info(f"  Reoriented {flipped_count} of {len(reads)} reads "
                f"({label}: {primary_count} by telomere, "
                f"{fallback_count} by satellite, {unoriented_count} unoriented)")
    return oriented_reads, flipped_reads


def orient_satellite_top(reads):
    """Reorient reads so satellite-dense end is at top.

    Uses telomere features as primary orientation signal, falling back to
    satellite density for reads without telomere.
    """
    return _orient_reads_with_fallback(
        reads, TELOMERE_FEATURES, SATELLITE_FEATURES, "satellite-dense"
    )


def assign_heatmap_colors(metadata_columns, read_metadata):
    """Assign colors to unique values in each metadata column.

    Returns:
        dict: {col_name: {value: hex_color}} mapping for each column.
    """
    palette = ["#40D392", "#60A5FA", "#F07167", "#FBBF24", "#EC4899", "#C4A9E8"]
    missing_color = "#545454"
    color_map = {}

    for col in metadata_columns:
        # Collect unique values preserving first-seen order
        seen = []
        seen_set = set()
        for rid, meta in read_metadata.items():
            val = meta.get(col)
            if val is not None and val not in seen_set:
                seen.append(val)
                seen_set.add(val)
        # Assign colors cyclically
        val_colors = {}
        for i, val in enumerate(seen):
            val_colors[val] = palette[i % len(palette)]
        val_colors[None] = missing_color
        color_map[col] = val_colors

    return color_map


def draw_heatmap_grid_svg(d, read_positions, config, text_color):
    """Draw heatmap grid of colored squares between labels and read bars (SVG)."""
    metadata_columns = config["metadata_columns"]
    read_metadata = config["read_metadata"]
    heatmap_colors = config["heatmap_colors"]
    heatmap_display_names = config.get("heatmap_display_names", {})
    bar_width = config["bar_width"]
    top_margin = config["top_margin"]
    heatmap_bottom_gap = config["heatmap_bottom_gap"]
    heatmap_row_gap = config["heatmap_row_gap"]
    background = config["background"]
    left_margin = config["left_margin"]
    font_size = config.get("font_size", 14)

    stroke_color = "#ffffff" if background == "black" else "#000000"
    n_rows = len(metadata_columns)

    for row_idx, col_name in enumerate(reversed(metadata_columns)):
        # Rows stack bottom-up from top_margin - heatmap_bottom_gap
        row_y = top_margin - heatmap_bottom_gap - (row_idx + 1) * bar_width - row_idx * heatmap_row_gap
        val_colors = heatmap_colors.get(col_name, {})

        for read_id, rx in read_positions:
            meta = read_metadata.get(read_id, {})
            val = meta.get(col_name)
            color = val_colors.get(val, val_colors.get(None, "#545454"))
            d.append(draw.Rectangle(
                rx, row_y, bar_width, bar_width,
                fill=color,
                stroke=stroke_color,
                stroke_width=0.5,
            ))

        # Row label at left margin (use display name if available)
        display_name = heatmap_display_names.get(col_name, col_name)
        label_y = row_y + bar_width / 2 + font_size * 0.35
        d.append(draw.Text(
            display_name, font_size,
            left_margin - 5, label_y,
            fill=text_color,
            text_anchor="end",
            font_family="Basic Sans",
        ))


def draw_heatmap_grid_png(img_array, read_positions, config, width_scale, height_scale):
    """Draw heatmap grid of colored squares between labels and read bars (PNG)."""
    metadata_columns = config["metadata_columns"]
    read_metadata = config["read_metadata"]
    heatmap_colors = config["heatmap_colors"]
    base_bar_width = config["bar_width"]
    base_top_margin = config["top_margin"]
    heatmap_bottom_gap = config["heatmap_bottom_gap"]
    heatmap_row_gap = config["heatmap_row_gap"]
    background = config["background"]

    bar_width = max(1, int(base_bar_width * width_scale))
    top_margin = int(base_top_margin * height_scale)
    scaled_bottom_gap = int(heatmap_bottom_gap * height_scale)
    scaled_row_gap = max(1, int(heatmap_row_gap * height_scale))
    box_h = max(1, int(base_bar_width * width_scale))

    border_color = np.array([255, 255, 255, 255] if background == "black"
                            else [0, 0, 0, 255], dtype=np.uint8)
    img_h, img_w = img_array.shape[:2]

    for row_idx, col_name in enumerate(reversed(metadata_columns)):
        row_y = top_margin - scaled_bottom_gap - (row_idx + 1) * box_h - row_idx * scaled_row_gap
        val_colors = heatmap_colors.get(col_name, {})

        for read_id, rx in read_positions:
            meta = read_metadata.get(read_id, {})
            val = meta.get(col_name)
            color_hex = val_colors.get(val, val_colors.get(None, "#545454"))
            color_rgba = hex_to_rgba(color_hex)

            x1 = max(0, rx)
            x2 = min(img_w, rx + bar_width)
            y1 = max(0, row_y)
            y2 = min(img_h, row_y + box_h)
            if x2 > x1 and y2 > y1:
                img_array[y1:y2, x1:x2] = color_rgba
                # Scaled border
                bt = max(1, int(width_scale))
                img_array[y1:y1 + bt, x1:x2] = border_color
                by = min(y2, img_h - 1)
                img_array[by - bt + 1:by + 1, x1:x2] = border_color
                img_array[y1:y2, x1:x1 + bt] = border_color
                bx = min(x2, img_w - 1)
                img_array[y1:y2, bx - bt + 1:bx + 1] = border_color


def draw_heatmap_legend_png(draw_ctx, config, top_margin, max_height_px, left_margin,
                            height_scale, width_scale):
    """Draw heatmap legend (colored swatches + value labels) below reads in PNG."""
    heatmap_colors = config["heatmap_colors"]
    metadata_columns = config["metadata_columns"]
    heatmap_display_names = config.get("heatmap_display_names", {})
    background = config["background"]
    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)

    legend_font_size = max(10, int(config.get("font_size", 14) * height_scale))
    font = _load_font(legend_font_size)
    swatch_size = legend_font_size
    col_gap = int(20 * width_scale)

    legend_x = left_margin
    legend_y = top_margin + max_height_px + int(10 * height_scale)

    for col_name in metadata_columns:
        val_colors = heatmap_colors.get(col_name, {})
        display_name = heatmap_display_names.get(col_name, col_name)
        header = display_name + ":"
        bold_font = _load_font(legend_font_size)
        draw_ctx.text((legend_x, legend_y), header, fill=text_color, font=bold_font)
        bbox = bold_font.getbbox(header)
        legend_x += bbox[2] - bbox[0] + int(10 * width_scale)
        for val, color in val_colors.items():
            if val is None:
                continue
            rgba = hex_to_rgba(color)
            draw_ctx.rectangle(
                [legend_x, legend_y, legend_x + swatch_size, legend_y + swatch_size],
                fill=rgba)
            label = str(val)
            draw_ctx.text((legend_x + swatch_size + int(3 * width_scale), legend_y),
                          label, fill=text_color, font=font)
            bbox = font.getbbox(label)
            legend_x += swatch_size + (bbox[2] - bbox[0]) + int(15 * width_scale)
        legend_x += col_gap

    return legend_font_size + int(10 * height_scale)


def draw_scale_bar(d, x, y, ratio, text_color, max_height_px, font_size=14, scale_bar_bp_override=None):
    """Draw a vertical scale bar showing read length scale."""
    scale_bar_bp = scale_bar_bp_override or 10000  # 10 kbp
    scale_bar_height = int(scale_bar_bp * ratio)

    # Don't draw if scale bar is too tall
    if not scale_bar_bp_override and scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    # Draw vertical bar
    d.append(draw.Line(x, y, x, y + scale_bar_height, stroke=text_color, stroke_width=1))

    # Draw tick marks at top and bottom
    d.append(draw.Line(x - 2, y, x + 6, y, stroke=text_color, stroke_width=1))
    d.append(
        draw.Line(
            x - 2,
            y + scale_bar_height,
            x + 6,
            y + scale_bar_height,
            stroke=text_color,
            stroke_width=1,
        )
    )

    # Draw label (rotated)
    label = f"{scale_bar_bp // 1000} Kbp"
    label_x = x - 5
    label_y = y + scale_bar_height / 2
    d.append(
        draw.Text(
            label,
            font_size,
            label_x,
            label_y,
            fill=text_color,
            text_anchor="middle",
            font_family="Basic Sans",
            transform=f"rotate(-90, {label_x}, {label_y})",
        )
    )


def draw_scale_bar_svg(output_path, image_height, top_margin, ratio, background, scale_bar_bp_override=None):
    """Draw a separate SVG file containing just the scale bar."""
    text_color = "#ffffff" if background == "black" else "#000000"
    max_height_px = image_height - top_margin - 40  # Approximate bottom margin

    # Calculate scale bar dimensions
    scale_bar_bp = scale_bar_bp_override or 10000  # 10 kbp
    scale_bar_height = int(scale_bar_bp * ratio)
    if not scale_bar_bp_override and scale_bar_height > max_height_px:
        scale_bar_bp = 5000
        scale_bar_height = int(scale_bar_bp * ratio)

    # Create a narrow SVG for the scale bar
    width = 60
    d = draw.Drawing(width, image_height, id_prefix="sb")
    d.append(draw.Rectangle(0, 0, width, image_height, fill=background))

    x = 45
    y = top_margin

    # Draw vertical bar
    d.append(draw.Line(x, y, x, y + scale_bar_height, stroke=text_color, stroke_width=1))

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
    logger.info(f"Saved scale bar: {output_path}")
    logger.info(f"  Dimensions: {width} x {image_height} pixels")

    # Also render PNG via rsvg-convert if available
    png_path = output_path.rsplit('.', 1)[0] + '.png'
    try:
        subprocess.run(
            ['rsvg-convert', '-o', png_path, output_path],
            check=True,
            capture_output=True,
        )
        logger.info(f"  Also rendered PNG via rsvg-convert: {png_path}")
    except FileNotFoundError:
        logger.info(f"  Note: rsvg-convert not found, skipping PNG render of scale bar")
    except subprocess.CalledProcessError:
        logger.warning(f"  Warning: rsvg-convert failed to render scale bar PNG")


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
    subgroup_spacing = config.get("subgroup_spacing", sample_spacing)
    group_boundaries = config.get("group_boundaries", set())
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
    markers = config.get("markers", {})
    if markers:
        _ms = config.get("marker_scale", 1.0)
        effective_left_margin += max(2, int(bar_width // 3 * _ms)) + 2  # Room for arrowheads
    n_group_gaps = sum(1 for s in sample_order[:-1] if s in group_boundaries) if len(sample_order) > 1 else 0
    n_subgroup_gaps = (num_samples - 1) - n_group_gaps
    image_width = (
        effective_left_margin
        + (total_reads * (bar_width + read_spacing))
        + (n_group_gaps * sample_spacing)
        + (n_subgroup_gaps * subgroup_spacing)
        + 50  # Right margin
    )
    # Pre-calculate legend heights to include in image_height
    features_used = {feat for _, _, _, feats in reads for _, _, feat in feats}
    has_heatmap = config.get("heatmap") and config.get("metadata_columns")
    has_featureset = bool(features_used)
    hm_items = _build_heatmap_legend_items(config) if has_heatmap else None
    legend_extra = 0
    if has_featureset or hm_items:
        legend_extra += estimate_featureset_legend_height(
            features_used, colors, config, image_width, extra_items=hm_items)

    image_height = top_margin + max_height_px + bottom_margin + legend_extra

    # Create drawing
    d = draw.Drawing(image_width, image_height, id_prefix="tr")
    text_color = "#ffffff" if background == "black" else "#000000"
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background))

    # Draw scale bar (left side) only if not creating separate file
    if draw_scale:
        draw_scale_bar(d, left_margin - 10, top_margin, ratio, text_color, max_height_px,
                        font_size=config.get("font_size", 14),
                        scale_bar_bp_override=config.get("scale_bar_bp"))

    # Use effective left margin for positioning reads
    left_margin = effective_left_margin

    # Track x position and sample boundaries
    current_x = left_margin
    current_sample = None
    sample_x_start = {}
    sample_x_end = {}

    feature_mode = config.get("feature_mode", "raw")
    min_feature_width = config.get("min_feature_width", 0.5)
    min_width_exclude = config.get("min_width_exclude", [])
    oversample = config.get("oversample", 1)
    read_border = config.get("read_border", False)
    markers = config.get("markers", {})
    marker_scale = config.get("marker_scale", 1.0)
    arrow_size = max(2, int(bar_width // 3 * marker_scale))
    read_positions = []

    for sample, read_id, read_length, features in reads:
        # Add sample spacing when sample changes
        if current_sample is not None and sample != current_sample:
            sample_x_end[current_sample] = current_x - read_spacing
            if current_sample in group_boundaries:
                current_x += sample_spacing
            else:
                current_x += subgroup_spacing

        if sample not in sample_x_start:
            sample_x_start[sample] = current_x

        current_sample = sample
        read_positions.append((read_id, current_x))

        if feature_mode == "raw":
            for start, end, feature in features:
                y_start = top_margin + int(start * ratio)
                height = max(1, int((end - start) * ratio))
                color = colors.get(feature, "#ffffff")
                d.append(
                    draw.Rectangle(current_x, y_start, bar_width, height, fill=color)
                )
        else:
            colored_feats = _build_colored_features(features, colors, min_width_exclude)
            max_feat_end = max((e for _, e, _ in features), default=0)
            bar_length = int(max_feat_end * ratio)
            rasterized = rasterize_features(colored_feats, bar_length, ratio,
                                            feature_mode, oversample, min_feature_width)
            if rasterized:
                for run in rasterized:
                    ry = top_margin + run['scaled_start']
                    rh = run['scaled_stop'] - run['scaled_start']
                    d.append(draw.Rectangle(current_x, ry, bar_width, rh,
                                            fill=run['color'],
                                            fill_opacity=run['fill_opacity']))

        if read_border:
            max_feat_end = max((e for _, e, _ in features), default=0)
            read_height_px = int(max_feat_end * ratio)
            d.append(draw.Rectangle(current_x, top_margin, bar_width, read_height_px,
                                    fill='none', stroke=text_color, stroke_width=0.5))

        if read_id in markers:
            for m_start, m_end in markers[read_id]:
                mid_y = top_margin + int((m_start + m_end) / 2 * ratio)
                d.append(draw.Lines(
                    current_x - arrow_size - 1, mid_y - arrow_size,
                    current_x - 1, mid_y,
                    current_x - arrow_size - 1, mid_y + arrow_size,
                    fill=text_color, close=True,
                ))

        current_x += bar_width + read_spacing

    # Record final sample boundary
    if current_sample:
        sample_x_end[current_sample] = current_x - read_spacing

    # Draw horizontal line, sample labels, and separators (unless --no-header)
    if not config.get("no_header", False):
        label_interval = 300  # Repeat label every N reads
        font_family = "Basic Sans"
        font_size = config.get("font_size", 14)
        label_tiers = config.get("label_tiers", {})
        group_subgroup_order = config.get("group_subgroup_order", [])
        has_grouping = bool(label_tiers)

        if has_grouping:
            tier_offset = font_size + 10
            hm_offset = config.get("heatmap_total", 0)
            tier_display_names = config.get("tier_display_names", {})

            # Tier 2 (bottom): subgroup lines + labels
            tier2_line_y = top_margin - 5 - hm_offset
            for sample in sample_order:
                if sample in sample_x_start and sample in sample_x_end:
                    d.append(draw.Line(
                        sample_x_start[sample], tier2_line_y,
                        sample_x_end[sample], tier2_line_y,
                        stroke=text_color, stroke_width=1.5,
                    ))
                    _, subgroup = label_tiers.get(sample, (sample, None))
                    label_text = (subgroup or sample).replace("_", " ")
                    d.append(draw.Text(
                        label_text, font_size,
                        sample_x_start[sample], top_margin - 10 - hm_offset,
                        fill=text_color, text_anchor="start",
                        font_family=font_family,
                    ))

            # Tier 2 display name (left margin)
            if 1 in tier_display_names:
                d.append(draw.Text(
                    tier_display_names[1], font_size,
                    left_margin - 5, top_margin - 10 - hm_offset,
                    fill=text_color, text_anchor="end",
                    font_family=font_family,
                ))

            # Tier 1 (top): group lines + labels spanning all subgroups
            group_spans = compute_group_spans(
                sample_order, sample_x_start, sample_x_end, group_subgroup_order)
            tier1_line_y = top_margin - 5 - tier_offset - hm_offset
            for group_name, gx_start, gx_end in group_spans:
                d.append(draw.Line(
                    gx_start, tier1_line_y, gx_end, tier1_line_y,
                    stroke=text_color, stroke_width=1.5,
                ))
                d.append(draw.Text(
                    group_name.replace("_", " "), font_size,
                    gx_start, tier1_line_y - 5,
                    fill=text_color, text_anchor="start",
                    font_family=font_family,
                ))

            # Tier 1 display name (left margin)
            if 0 in tier_display_names:
                d.append(draw.Text(
                    tier_display_names[0], font_size,
                    left_margin - 5, tier1_line_y - 5,
                    fill=text_color, text_anchor="end",
                    font_family=font_family,
                ))
        else:
            # Single-tier labels (no grouping)
            hm_offset = config.get("heatmap_total", 0)
            for sample in sample_order:
                if sample in sample_x_start and sample in sample_x_end:
                    d.append(draw.Line(
                        sample_x_start[sample], top_margin - 5 - hm_offset,
                        sample_x_end[sample], top_margin - 5 - hm_offset,
                        stroke=text_color, stroke_width=2,
                    ))
                    label_text = sample.replace("_", " ")
                    num_reads_in_sample = sample_read_counts[sample]
                    read_width = bar_width + read_spacing
                    for i in range(0, num_reads_in_sample, label_interval):
                        x_pos = sample_x_start[sample] + (i * read_width)
                        d.append(draw.Text(
                            label_text, font_size, x_pos, top_margin - 12 - hm_offset,
                            fill=text_color, text_anchor="start",
                            font_family=font_family,
                        ))

    # Draw heatmap grid if enabled
    if config.get("heatmap") and config.get("metadata_columns"):
        draw_heatmap_grid_svg(d, read_positions, config, text_color)

    # Draw unified color legend below reads (heatmap + featureset sections)
    if features_used or hm_items:
        feat_legend_y = top_margin + max_height_px + 10
        draw_legend_svg(d, features_used, colors, config, feat_legend_y, image_width,
                        extra_items=hm_items)

    # Save SVG
    d.save_svg(output_path)
    logger.info(f"\nSaved: {output_path}")
    logger.info(f"  Dimensions: {image_width} x {image_height} pixels")

    return image_height


def main():
    args = parse_args()
    t_start = time.time()

    # Hardcoded layout margins (not user-facing)
    args.top_margin = 30
    args.left_margin = 30
    args.bottom_margin = 15

    # Set up logging (console + file)
    log_path = setup_logging(args.output)

    logger.info("KaryoScope Read Visualization")
    logger.info("=" * 50)
    logger.info("Command: %s", " ".join(sys.argv))
    logger.debug("Output: %s", args.output)

    # If --animate is set, ensure PNG format
    fmt = args.format
    if args.animate and fmt == "svg":
        logger.info("  --animate requires PNG output, switching format to 'both'")
        fmt = "both"

    # Validate input arguments
    if args.bed:
        input_mode = "bed"
    elif args.clusters and args.cluster_prefix:
        input_mode = "cluster"
        if not args.results_dir:
            logger.error("Error: --clusters requires --results-dir")
            return
    elif args.samples and args.results_dir:
        input_mode = "samples"
    else:
        logger.error("Error: Must provide --bed, --samples+--results-dir, "
                      "or --clusters+--cluster-prefix+--results-dir")
        return

    # Load color mapping
    logger.info("\nLoading colors from: %s", args.colors)
    colors, color_sections = load_color_mapping(args.colors)
    logger.info("  Loaded %d color mappings", len(colors))

    # Load BED data
    if input_mode == "bed":
        logger.info("\nLoading feature data from BED files...")
        reads, sample_order = load_bed_files_direct(args.bed)
    elif input_mode == "cluster":
        cluster_ids = [int(c.strip()) for c in args.clusters.split(',')]
        logger.info("\nLoading reads for clusters: %s", cluster_ids)
        reads, sample_order = load_cluster_reads(
            args.cluster_prefix, cluster_ids, args.results_dir,
            args.database, args.featureset, args.smoothness, args.analysis)
    else:
        logger.info("\nLoading feature data from: %s", args.results_dir)
        reads = load_sample_bed_data(
            args.samples,
            args.results_dir,
            args.database,
            args.featureset,
            args.smoothness,
            args.analysis,
        )
        sample_order = args.samples

    logger.info("  Total reads: %d", len(reads))

    # Length filtering
    if args.min_length is not None or args.max_length is not None:
        before = len(reads)
        reads = [r for r in reads
                 if (args.min_length is None or r[2] >= args.min_length)
                 and (args.max_length is None or r[2] <= args.max_length)]
        logger.info("\nLength filter: %d -> %d reads (min=%s, max=%s)",
                     before, len(reads), args.min_length, args.max_length)

    # Filter to read list if provided (with optional grouping)
    all_columns = []
    read_data = {}
    read_groups = {}
    group_subgroup_order = []
    read_metadata = {}
    metadata_columns = []
    tier_specs = []   # [(column_name, display_name), ...]
    heatmap_specs = []  # [(column_name, display_name), ...]
    if args.read_list:
        read_ids, all_columns, read_data, read_rows = load_read_list(args.read_list)
        logger.info("\nFiltering to %d read IDs from: %s", len(read_ids), args.read_list)
        reads = [r for r in reads if r[1] in read_ids]
        logger.info("  Matched: %d reads", len(reads))
        if all_columns:
            logger.info("  Columns: %s", all_columns)

        # Parse --label-tier and --heatmap-track
        if args.label_tier:
            for spec in args.label_tier:
                col, _, name = spec.partition(":")
                tier_specs.append((col, name or col))
        elif all_columns and len(all_columns) >= 2:
            # Default: first two columns as tiers (group, subgroup)
            tier_specs = [(all_columns[0], all_columns[0]),
                          (all_columns[1], all_columns[1])]

        if args.heatmap_track:
            for spec in args.heatmap_track:
                col, _, name = spec.partition(":")
                heatmap_specs.append((col, name or col))
        elif args.heatmap and all_columns:
            # Default: remaining columns after tiers
            tier_cols = {t[0] for t in tier_specs}
            heatmap_specs = [(c, c) for c in all_columns if c not in tier_cols]

        # Validate column names exist in the header
        if all_columns:
            col_set = set(all_columns)
            for col, _ in tier_specs:
                if col not in col_set:
                    logger.warning("  --label-tier column '%s' not found in TSV header: %s", col, all_columns)
            for col, _ in heatmap_specs:
                if col not in col_set:
                    logger.warning("  --heatmap-track column '%s' not found in TSV header: %s", col, all_columns)

        # Build read_groups and group_subgroup_order from tier_specs
        # Iterate read_rows (raw file order) to preserve TSV ordering,
        # even when duplicate read IDs cause read_data dict entries to
        # be overwritten.
        if len(tier_specs) >= 2:
            group_col = tier_specs[0][0]
            subgroup_col = tier_specs[1][0]
            seen_pairs = set()
            for rid, meta in read_rows:
                group = meta.get(group_col)
                subgroup = meta.get(subgroup_col)
                if group and subgroup:
                    read_groups[rid] = (group, subgroup)
                    pair = (group, subgroup)
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        group_subgroup_order.append(pair)
                elif group:
                    read_groups[rid] = (group, None)
                    pair = (group, None)
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        group_subgroup_order.append(pair)
        elif len(tier_specs) == 1:
            group_col = tier_specs[0][0]
            seen_pairs = set()
            for rid, meta in read_rows:
                group = meta.get(group_col)
                if group:
                    read_groups[rid] = (group, None)
                    pair = (group, None)
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        group_subgroup_order.append(pair)

        # Build read_metadata and metadata_columns from heatmap_specs
        if heatmap_specs:
            metadata_columns = [col for col, _ in heatmap_specs]
            for rid in read_ids:
                src = read_data.get(rid, {})
                meta = {}
                for col, _ in heatmap_specs:
                    meta[col] = src.get(col)
                read_metadata[rid] = meta

        if read_groups:
            logger.info("  Grouping: %d group-subgroup pairs", len(group_subgroup_order))
        if metadata_columns:
            logger.info("  Heatmap columns: %s", metadata_columns)

        # --filter-group: keep only reads whose group matches specified values.
        # Must scan read_rows (not read_groups) because a read appearing in
        # multiple groups has its read_groups entry overwritten by the last one.
        if args.filter_group and tier_specs:
            allowed = set(args.filter_group)
            group_col = tier_specs[0][0]
            subgroup_col = tier_specs[1][0] if len(tier_specs) >= 2 else None
            # Collect IDs and rebuild read_groups from matching rows only
            keep_ids = set()
            filtered_groups = {}
            seen_pairs = set()
            filtered_order = []
            for rid, meta in read_rows:
                g = meta.get(group_col)
                if g not in allowed:
                    continue
                keep_ids.add(rid)
                s = meta.get(subgroup_col) if subgroup_col else None
                filtered_groups[rid] = (g, s)
                pair = (g, s)
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    filtered_order.append(pair)
            before = len(reads)
            reads = [r for r in reads if r[1] in keep_ids]
            read_groups = filtered_groups
            group_subgroup_order = filtered_order
            logger.info("  --filter-group %s: %d -> %d reads, %d groups",
                        ", ".join(sorted(allowed)), before, len(reads),
                        len(group_subgroup_order))

    # Sort group_subgroup_order: preserve first-seen order for both groups and subgroups
    if group_subgroup_order and len(tier_specs) >= 2:
        group_rank = {}
        subgroup_rank = {}
        for g, s in group_subgroup_order:
            if g not in group_rank:
                group_rank[g] = len(group_rank)
            if s not in subgroup_rank:
                subgroup_rank[s] = len(subgroup_rank)
        group_subgroup_order.sort(key=lambda pair: (group_rank.get(pair[0], 0), subgroup_rank.get(pair[1], 0)))

    # Validate heatmap option
    if args.heatmap or heatmap_specs:
        if heatmap_specs and not args.heatmap:
            args.heatmap = True  # --heatmap-track implies --heatmap
        if not metadata_columns:
            logger.warning("--heatmap requires heatmap columns in --read-list TSV. Disabling heatmap.")
            args.heatmap = False
        elif args.horizontal:
            logger.warning("--heatmap is not supported with --horizontal. Disabling heatmap.")
            args.heatmap = False

    if not reads:
        logger.error("No reads loaded. Check sample names and paths.")
        return

    # Orient reads so telomere is at top (if requested)
    flipped_reads = {}
    if args.orient_telomere_top:
        logger.info("\nOrienting reads (telomere at top)...")
        reads, flipped_reads = orient_telomere_top(reads)

    # Orient reads so chromosome is at top (if requested)
    if args.orient_chromosome_top:
        logger.info("\nOrienting reads (chromosome at top)...")
        reads, flipped_reads = orient_chromosome_top(reads)

    # Orient reads so satellite-dense end is at top (if requested)
    if args.orient_satellite_top:
        logger.info("\nOrienting reads (satellite-dense end at top)...")
        reads, flipped_reads = orient_satellite_top(reads)

    # Apply two-level grouping if read-list provided group/subgroup columns
    group_boundaries = set()
    label_tiers = {}
    tier_display_names = {}  # {tier_index: display_name}
    if read_groups and group_subgroup_order:
        reads, sample_order, group_boundaries = apply_read_list_grouping(
            reads, read_groups, group_subgroup_order)
        # Build label_tiers lookup: composite_label -> (group, subgroup)
        def _label(g, s):
            return f"{g} \u2014 {s}" if s else g
        label_tiers = {_label(g, s): (g, s) for g, s in group_subgroup_order}
        # Build tier display names from specs
        for i, (_, display_name) in enumerate(tier_specs):
            tier_display_names[i] = display_name
        # Add extra top margin for additional tiers beyond the first
        n_tiers = len(tier_specs)
        if n_tiers > 0:
            tier_offset = args.font_size + 10
            args.top_margin += (n_tiers - 1) * tier_offset
        logger.info("  Applied grouping: %d composite samples, %d tier(s)", len(sample_order), n_tiers)
        if tier_display_names:
            # Estimate left margin from longest tier display name
            max_label_len = max((len(dn) for _, dn in tier_specs), default=0)
            needed_margin = max(60, int(max_label_len * 8.5) + 10)
            args.left_margin = max(args.left_margin, needed_margin)

    # Build heatmap display name mapping
    heatmap_display_names = {}  # {column_name: display_name}
    for col, display_name in heatmap_specs:
        heatmap_display_names[col] = display_name

    # Heatmap layout: assign colors now; margin adjustment deferred until after
    # viewport ratio scaling so we use the final bar_width.
    heatmap_colors = {}
    heatmap_row_gap = 1
    heatmap_top_gap = 3
    heatmap_bottom_gap = 10
    heatmap_total = 0
    if args.heatmap and metadata_columns:
        args.bottom_margin = max(args.bottom_margin, 35)  # Room for heatmap + featureset legends
        heatmap_colors = assign_heatmap_colors(metadata_columns, read_metadata)
        for col, val_map in heatmap_colors.items():
            display = heatmap_display_names.get(col, col)
            cats = [f"{v}={c}" for v, c in val_map.items() if v is not None]
            logger.info("  %s: %s", display, ", ".join(cats))

    # Filter by max read length
    if args.max_read_length:
        before = len(reads)
        reads = [(s, rid, rl, f) for s, rid, rl, f in reads
                 if rl <= args.max_read_length]
        removed = before - len(reads)
        if removed:
            logger.info("\nFiltered %d read(s) exceeding %d bp", removed, args.max_read_length)

    # Load marker positions (for arrowheads) — multiple markers per read supported
    markers = defaultdict(list)
    if args.markers:
        n_markers = 0
        with open(args.markers) as mf:
            header = mf.readline().strip().split('\t')
            for line in mf:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    markers[parts[0]].append((int(parts[1]), int(parts[2])))
                    n_markers += 1
        # Flip marker positions for reads that were reoriented
        if flipped_reads:
            flipped_marker_count = 0
            for read_id in list(markers.keys()):
                if read_id in flipped_reads:
                    rl = flipped_reads[read_id]
                    markers[read_id] = [(rl - m_end, rl - m_start)
                                        for m_start, m_end in markers[read_id]]
                    flipped_marker_count += len(markers[read_id])
            if flipped_marker_count:
                logger.info("  Flipped %d marker(s) to match read orientation",
                            flipped_marker_count)
        logger.info("\nLoaded %d marker positions for %d reads from: %s",
                     n_markers, len(markers), args.markers)

    # Sort reads
    logger.info("\nSorting reads by sample order, then by length (descending)")
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
        "subgroup_spacing": args.subgroup_spacing if args.subgroup_spacing is not None else args.sample_spacing,
        "group_boundaries": group_boundaries,
        "sample_order": sample_order,
        "draw_scale_bar": args.scale_bar_output is None and not args.no_scale_bar,
        "scale_bar_bp": args.scale_bar_bp,
        "no_header": args.no_header,
        "png_scale": args.png_scale,
        "font_size": args.font_size,
        "label_tiers": label_tiers,
        "group_subgroup_order": group_subgroup_order,
        "tier_display_names": tier_display_names,
        "heatmap_display_names": heatmap_display_names,
        "feature_mode": args.feature_mode,
        "min_feature_width": args.min_feature_width,
        "min_width_exclude": args.min_width_exclude or [],
        "oversample": args.oversample,
        "read_border": args.read_border,
        "heatmap": args.heatmap,
        "read_metadata": read_metadata,
        "metadata_columns": metadata_columns,
        "heatmap_colors": heatmap_colors,
        "heatmap_row_gap": heatmap_row_gap,
        "heatmap_top_gap": heatmap_top_gap,
        "heatmap_bottom_gap": heatmap_bottom_gap,
        "heatmap_total": heatmap_total,
        "color_sections": color_sections,
        "markers": markers,
        "marker_scale": args.marker_scale,
    }

    # Apply viewport ratio overrides
    # Only the ratio (px/bp) changes with aspect ratio — bar_width, read_spacing,
    # heatmap squares, and labels stay constant regardless of ratio.
    if args.viewport_ratio and args.viewport_ratio.lower() != "none":
        vr_w, vr_h = [int(x) for x in args.viewport_ratio.split(':')]
        total_reads_count = len(sorted_reads)
        _sample_counts = defaultdict(int)
        for s, _, _, _ in sorted_reads:
            _sample_counts[s] += 1
        _num_samples = len(_sample_counts)

        bw = config["bar_width"]
        _so = config["sample_order"]
        _gb = config.get("group_boundaries", set())
        _n_group_gaps = sum(1 for s in _so[:-1] if s in _gb) if len(_so) > 1 else 0
        _n_subgroup_gaps = (_num_samples - 1) - _n_group_gaps
        natural_w = (
            config["left_margin"]
            + total_reads_count * (bw + config["read_spacing"])
            + _n_group_gaps * config["sample_spacing"]
            + _n_subgroup_gaps * config.get("subgroup_spacing", config["sample_spacing"])
            + 50
        )

        # Reserve heatmap space before computing ratio (needs final bar_width)
        if args.heatmap and metadata_columns:
            n_rows = len(metadata_columns)
            heatmap_total = (n_rows * bw
                             + (n_rows - 1) * heatmap_row_gap
                             + heatmap_top_gap + heatmap_bottom_gap)
            config["top_margin"] += heatmap_total
            config["heatmap_total"] = heatmap_total
            logger.info("\nHeatmap: %d row(s), %dpx added to top margin", n_rows, heatmap_total)

        # Derive height from aspect ratio; only adjust ratio (px/bp)
        natural_h = int(natural_w * vr_h / vr_w)
        max_length = max(r[2] for r in sorted_reads)
        usable_h = natural_h - config["top_margin"] - config["bottom_margin"]
        config["ratio"] = usable_h / max_length if max_length > 0 else config["ratio"]

        logger.info("\nViewport ratio %d:%d -> %dx%d: bar_width=%d, read_spacing=%d, ratio=%.6f",
                    vr_w, vr_h, natural_w, natural_h,
                    config['bar_width'], config['read_spacing'], config['ratio'])
    else:
        # No viewport ratio — still need heatmap margin if applicable
        if args.heatmap and metadata_columns:
            final_bw = config["bar_width"]
            n_rows = len(metadata_columns)
            heatmap_total = (n_rows * final_bw
                             + (n_rows - 1) * heatmap_row_gap
                             + heatmap_top_gap + heatmap_bottom_gap)
            config["top_margin"] += heatmap_total
            config["heatmap_total"] = heatmap_total
            logger.info("\nHeatmap: %d row(s), %dpx added to top margin", n_rows, heatmap_total)

    # Log full config at DEBUG level
    logger.debug("Config: %s", {k: v for k, v in config.items()
                 if k not in ('read_metadata',)})
    for sample in sample_order:
        count = sum(1 for s, _, _, _ in sorted_reads if s == sample)
        logger.debug("  %s: %d reads", sample, count)

    # Determine background theme(s) to render
    bg_themes = ['black', 'white'] if args.background == 'both' else [args.background]
    base_output = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
    animation_done = False

    for bg_color in bg_themes:
        config["background"] = bg_color
        theme_suffix = f"_{bg_color}" if len(bg_themes) > 1 else ""

        # Determine output paths for this theme
        svg_output = f"{base_output}{theme_suffix}.svg"
        png_output = f"{base_output}{theme_suffix}.png"

        image_size = None  # height for vertical, width for horizontal

        if args.batch_size:
            num_batches = (len(sorted_reads) + args.batch_size - 1) // args.batch_size
            logger.info("\nGenerating %d batched outputs (%d reads each) [%s]",
                        num_batches, args.batch_size, bg_color)

            for i in range(0, len(sorted_reads), args.batch_size):
                batch = sorted_reads[i:i + args.batch_size]
                batch_num = (i // args.batch_size) + 1
                batch_base = f"{base_output}{theme_suffix}.batch{batch_num:03d}"

                if args.horizontal:
                    if fmt in ("svg", "both"):
                        batch_svg = f"{batch_base}.svg"
                        logger.info("\nGenerating SVG (horizontal): %s", batch_svg)
                        image_size = draw_reads_horizontal_svg(batch, colors, batch_svg, config)
                    if fmt in ("png", "both"):
                        batch_png = f"{batch_base}.png"
                        logger.info("\nGenerating PNG (horizontal): %s", batch_png)
                        image_size = draw_reads_horizontal_png(batch, colors, batch_png, config)
                else:
                    if fmt in ("svg", "both"):
                        batch_svg = f"{batch_base}.svg"
                        logger.info("\nGenerating SVG: %s", batch_svg)
                        t_draw = time.time()
                        image_size = draw_reads_svg(batch, colors, batch_svg, config)
                        logger.debug("SVG generation took %.2fs", time.time() - t_draw)
                    if fmt in ("png", "both"):
                        batch_png = f"{batch_base}.png"
                        logger.info("\nGenerating PNG: %s", batch_png)
                        t_draw = time.time()
                        image_size = draw_reads_png(batch, colors, batch_png, config)
                        logger.debug("PNG generation took %.2fs", time.time() - t_draw)
        else:
            if args.horizontal:
                if fmt in ("svg", "both"):
                    logger.info("\nGenerating SVG (horizontal): %s", svg_output)
                    image_size = draw_reads_horizontal_svg(sorted_reads, colors, svg_output, config)
                if fmt in ("png", "both"):
                    logger.info("\nGenerating PNG (horizontal): %s", png_output)
                    image_size = draw_reads_horizontal_png(sorted_reads, colors, png_output, config)
            else:
                if fmt in ("svg", "both"):
                    logger.info("\nGenerating SVG: %s", svg_output)
                    t_draw = time.time()
                    image_size = draw_reads_svg(sorted_reads, colors, svg_output, config)
                    logger.debug("SVG generation took %.2fs", time.time() - t_draw)
                if fmt in ("png", "both"):
                    logger.info("\nGenerating PNG: %s", png_output)
                    t_draw = time.time()
                    image_size = draw_reads_png(sorted_reads, colors, png_output, config)
                    logger.debug("PNG generation took %.2fs", time.time() - t_draw)

        # Handle legend: save separately for animation, or composite onto PNG
        legend_png_path = None
        if fmt in ("png", "both"):
            features_used = {feat for _, _, _, feats in sorted_reads for _, _, feat in feats}
            height_scale, width_scale = calculate_png_params(sorted_reads, config)
            png_img = Image.open(png_output)
            png_width = png_img.size[0]
            png_img.close()
            hm_items = _build_heatmap_legend_items(config) if config.get("heatmap") else None
            legend_img = draw_legend_png(features_used, colors, config,
                                         horizontal=args.horizontal,
                                         scale=height_scale,
                                         max_width=png_width,
                                         extra_items=hm_items)
            if legend_img is not None:
                if args.animate:
                    # Save as separate file — animation function will overlay it
                    legend_png_path = f"{base_output}{theme_suffix}.legend.png"
                    if legend_img.mode == 'RGBA':
                        bg_c = (0, 0, 0) if bg_color == "black" else (255, 255, 255)
                        bg_img = Image.new('RGB', legend_img.size, bg_c)
                        bg_img.paste(legend_img, mask=legend_img.split()[3])
                        bg_img.save(legend_png_path)
                    else:
                        legend_img.save(legend_png_path)
                    logger.info("  Legend saved for animation: %s", legend_png_path)
                else:
                    # Composite directly onto the reads image
                    position = "right" if args.horizontal else "below"
                    composite_legend(png_output, legend_img, position=position)

        # Scale bar: adaptive horizontal panning draws its own; vertical needs a PNG overlay
        scale_bar_png_path = None
        if args.animate and not args.horizontal and not args.no_scale_bar and image_size:
            scale_bar_png_path = f"{base_output}{theme_suffix}.scalebar.png"
            height_scale, _ = calculate_png_params(sorted_reads, config)
            scaled_ratio = args.ratio * height_scale
            scaled_top_margin = int(args.top_margin * height_scale)
            draw_scale_bar_png(image_size, scaled_top_margin, scaled_ratio,
                               bg_color, scale_bar_png_path,
                               scale_bar_bp_override=config.get("scale_bar_bp"),
                               height_scale=height_scale,
                               font_size=config.get("font_size", 11))

        # Generate animation if requested (only once, for the first theme)
        if args.animate and not animation_done:
            if fmt not in ("png", "both"):
                logger.info("\nGenerating PNG for animation: %s", png_output)
                if args.horizontal:
                    draw_reads_horizontal_png(sorted_reads, colors, png_output, config)
                else:
                    draw_reads_png(sorted_reads, colors, png_output, config)
            _run_animation(png_output, len(sorted_reads), args, config,
                           legend_path=legend_png_path, scale_bar_path=scale_bar_png_path)
            animation_done = True

        # Generate separate scale bar if requested
        if args.scale_bar_output and image_size and not args.horizontal:
            sb_output = args.scale_bar_output
            if len(bg_themes) > 1:
                sb_base = sb_output.rsplit('.', 1)[0] if '.' in sb_output else sb_output
                sb_ext = sb_output.rsplit('.', 1)[1] if '.' in sb_output else 'svg'
                sb_output = f"{sb_base}_{bg_color}.{sb_ext}"
            draw_scale_bar_svg(
                sb_output,
                image_size,
                args.top_margin,
                args.ratio,
                bg_color,
                scale_bar_bp_override=config.get("scale_bar_bp"),
            )

    elapsed = time.time() - t_start
    logger.info("\nDone! (%.1fs)", elapsed)
    logger.info("Log file: %s", log_path)


if __name__ == "__main__":
    main()
