#!/usr/bin/env python3
"""
create_panning_animation.py

Create panning animations from wide or tall images.
Supports SVG and PNG inputs, fixed and adaptive zoom modes.
Designed for PowerPoint 16:9 widescreen presentations.

Zoom modes:
- fixed: classic seamless looping pan with optional crop ratio
- adaptive: dynamically zooms to keep tallest visible read filling viewport

Duration modes:
- --duration N: fixed duration in seconds
- --num-reads N --reads-per-second R: auto-calculate from read count (default)

Usage:
    # Fixed zoom, auto-duration from reads
    python create_panning_animation.py \\
        --input reads.svg --output pan.mp4 \\
        --num-reads 35412 --reads-per-second 30

    # Adaptive zoom with legend
    python create_panning_animation.py \\
        --input reads.svg --output pan.mp4 \\
        --zoom adaptive --num-reads 35412 \\
        --legend legend.svg
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

# Allow very large images
Image.MAX_IMAGE_PIXELS = None

# Playwright browser instance (lazy loaded)
_playwright = None
_browser = None


def _get_playwright_browser():
    """Get or create a Playwright browser instance for SVG rendering."""
    global _playwright, _browser
    if _browser is None:
        from playwright.sync_api import sync_playwright
        _playwright = sync_playwright().start()
        _browser = _playwright.chromium.launch()
    return _browser


def _close_playwright():
    """Close Playwright browser instance."""
    global _playwright, _browser
    if _browser:
        _browser.close()
        _browser = None
    if _playwright:
        _playwright.stop()
        _playwright = None


def render_svg_region_playwright(svg_path, output_path, x_offset, y_offset, width, height, scale=1.0):
    """Render a region of an SVG using Playwright/Chromium.

    Args:
        svg_path: Path to SVG file
        output_path: Output PNG path
        x_offset: Left offset in SVG coordinates
        y_offset: Top offset in SVG coordinates
        width: Width to render in SVG coordinates
        height: Height to render in SVG coordinates
        scale: Scale factor (zoom)

    Returns True on success, False on failure.
    """
    global _playwright, _browser
    page = None
    try:
        browser = _get_playwright_browser()
        page = browser.new_page()

        # Get SVG dimensions
        svg_w, svg_h = get_svg_dimensions(svg_path)
        if svg_w is None:
            svg_w, svg_h = 1000, 1000

        # Scale the viewport to get the desired output size
        viewport_w = int(width * scale)
        viewport_h = int(height * scale)

        # Use file:// URL to load SVG directly
        abs_svg_path = os.path.abspath(svg_path)

        # Set viewport large enough to fit the scaled SVG
        full_w = int(svg_w * scale)
        full_h = int(svg_h * scale)
        page.set_viewport_size({"width": max(viewport_w, full_w), "height": max(viewport_h, full_h)})

        # Load SVG directly via file:// URL
        page.goto(f"file://{abs_svg_path}", wait_until="load", timeout=120000)

        # Wait for SVG to render
        page.wait_for_timeout(3000)

        # Apply zoom via CSS transform on the root SVG element
        page.evaluate(f'''
            const svg = document.querySelector('svg');
            if (svg) {{
                svg.style.transform = 'scale({scale})';
                svg.style.transformOrigin = '0 0';
            }}
        ''')
        page.wait_for_timeout(500)

        # Take a clip screenshot of the region we want
        page.screenshot(
            path=output_path,
            type="png",
            clip={
                "x": int(x_offset * scale),
                "y": int(y_offset * scale),
                "width": viewport_w,
                "height": viewport_h
            }
        )
        page.close()
        page = None

        return True
    except Exception as e:
        print(f"    Playwright render error: {e}")
        # Try to restart the browser if it crashed
        if page:
            try:
                page.close()
            except:
                pass
        if _browser:
            try:
                _browser.close()
            except:
                pass
        _browser = None
        return False


def render_svg_to_png_playwright(svg_path, output_path, target_width=None, target_height=None):
    """Render an entire SVG to PNG using Playwright/Chromium.

    This is used as a fallback for SVGs that exceed rsvg-convert's element limit.
    For very wide SVGs, renders in horizontal tiles and stitches them together.

    Args:
        svg_path: Path to SVG file
        output_path: Output PNG path
        target_width: Target width (optional, defaults to max safe width)
        target_height: Target height (optional, defaults to native)

    Returns True on success, False on failure.
    """
    global _playwright, _browser

    # Get SVG dimensions
    svg_w, svg_h = get_svg_dimensions(svg_path)
    if svg_w is None:
        svg_w, svg_h = 1000, 1000

    # Chromium's max texture size
    max_dim = 16384

    # Calculate scale to fit height within limits (prioritize vertical resolution)
    if target_height and target_height <= max_dim:
        scale = target_height / svg_h
    elif target_width and target_width <= max_dim:
        scale = target_width / svg_w
    else:
        # Scale so height fits within max_dim
        scale = min(max_dim / svg_h, 1.0)

    output_w = int(svg_w * scale)
    output_h = min(int(svg_h * scale), max_dim)

    # If width fits in one tile, render directly
    if output_w <= max_dim:
        return _render_svg_single_tile(svg_path, output_path, output_w, output_h, scale)

    # Otherwise, render in horizontal tiles and stitch
    print(f"    Playwright: rendering {svg_w:.0f}x{svg_h:.0f} in tiles (scale={scale:.3f})")

    tile_width = max_dim
    num_tiles = (output_w + tile_width - 1) // tile_width
    tiles = []

    try:
        browser = _get_playwright_browser()

        for i in range(num_tiles):
            tile_start_px = i * tile_width  # Start position in output pixels
            tile_end_px = min((i + 1) * tile_width, output_w)
            tile_w = tile_end_px - tile_start_px

            # Convert to SVG coordinates
            svg_x_start = tile_start_px / scale
            svg_x_end = tile_end_px / scale

            # Create a page for this tile
            page = browser.new_page()
            page.set_viewport_size({"width": tile_w, "height": output_h})

            abs_svg_path = os.path.abspath(svg_path)
            page.goto(f"file://{abs_svg_path}", wait_until="load", timeout=600000)  # 10 min for large SVGs
            page.wait_for_timeout(3000)

            # Scale and position SVG
            page.evaluate(f'''
                const svg = document.querySelector('svg');
                if (svg) {{
                    svg.style.width = '{output_w}px';
                    svg.style.height = '{output_h}px';
                    svg.style.marginLeft = '-{tile_start_px}px';
                }}
            ''')
            page.wait_for_timeout(500)

            # Screenshot this tile
            tile_path = f"{output_path}.tile{i}.png"
            page.screenshot(path=tile_path, type="png", full_page=False)
            page.close()

            tiles.append((tile_path, tile_start_px, tile_w))
            print(f"      Tile {i+1}/{num_tiles}: x={tile_start_px}-{tile_end_px}")

        # Stitch tiles together
        result = Image.new("RGB", (output_w, output_h), (0, 0, 0))
        for tile_path, x_offset, tile_w in tiles:
            tile_img = Image.open(tile_path)
            result.paste(tile_img, (x_offset, 0))
            tile_img.close()
            os.unlink(tile_path)

        result.save(output_path)
        result.close()

        print(f"    Stitched {num_tiles} tiles → {output_w}x{output_h}")
        return True

    except Exception as e:
        print(f"    Playwright render error: {e}")
        # Clean up any partial tiles
        for tile_path, _, _ in tiles:
            if os.path.exists(tile_path):
                os.unlink(tile_path)
        if _browser:
            try:
                _browser.close()
            except:
                pass
        _browser = None
        return False


def _render_svg_single_tile(svg_path, output_path, output_w, output_h, scale):
    """Render a single-tile SVG to PNG."""
    global _playwright, _browser
    page = None
    try:
        browser = _get_playwright_browser()
        page = browser.new_page()

        abs_svg_path = os.path.abspath(svg_path)
        page.set_viewport_size({"width": output_w, "height": output_h})
        page.goto(f"file://{abs_svg_path}", wait_until="load", timeout=180000)
        page.wait_for_timeout(5000)

        # Apply scale via CSS
        page.evaluate(f'''
            const svg = document.querySelector('svg');
            if (svg) {{
                svg.style.width = '{output_w}px';
                svg.style.height = '{output_h}px';
            }}
        ''')
        page.wait_for_timeout(1000)

        page.screenshot(path=output_path, type="png", full_page=False)
        page.close()

        print(f"    Playwright: rendered single tile {output_w}x{output_h}")
        return True
    except Exception as e:
        print(f"    Playwright render error: {e}")
        if page:
            try:
                page.close()
            except:
                pass
        if _browser:
            try:
                _browser.close()
            except:
                pass
        _browser = None
        return False


# ─── SVG helpers ──────────────────────────────────────────────────────────────

def get_svg_dimensions(svg_path):
    """Parse SVG width and height from file header."""
    with open(svg_path) as f:
        header = f.read(4000)
    w_match = re.search(r'width="([\d.]+)"', header)
    h_match = re.search(r'height="([\d.]+)"', header)
    if w_match and h_match:
        return float(w_match.group(1)), float(h_match.group(1))
    return None, None


def svg_to_png(svg_path, height=None, width=None):
    """Convert SVG to temporary PNG using rsvg-convert or cairosvg fallback.

    Args:
        svg_path: Path to SVG file
        height: Target height in pixels (optional)
        width: Target width in pixels (optional)

    Returns (temp_png_path, scale_factor).
    """
    svg_w, svg_h = get_svg_dimensions(svg_path)

    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    tmp.close()

    print(f"  Converting SVG → PNG: {os.path.basename(svg_path)}")
    if svg_w and svg_h:
        print(f"    SVG dimensions: {svg_w:.0f} x {svg_h:.0f}")
    if height:
        print(f"    Target height: {height}")
    if width:
        print(f"    Target width: {width}")

    # Try rsvg-convert first
    cmd = ["rsvg-convert", "-o", tmp.name]
    if height:
        cmd.extend(["-h", str(height)])
    if width:
        cmd.extend(["-w", str(width)])
    cmd.append(svg_path)

    success = False
    try:
        env = os.environ.copy()
        env["RSVG_MAX_LOADED_ELEMENTS"] = "100000000"
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        success = True
    except subprocess.CalledProcessError as e:
        # Fallback to cairosvg for large SVGs
        print(f"    rsvg-convert failed, trying cairosvg...")
        try:
            import cairosvg
            if width and svg_w:
                scale_factor = width / svg_w
            elif height and svg_h:
                scale_factor = height / svg_h
            else:
                scale_factor = 1.0
            cairosvg.svg2png(
                url=svg_path,
                write_to=tmp.name,
                scale=scale_factor,
            )
            success = True
        except ImportError:
            pass
        except Exception as cairo_err:
            pass

    # Fallback to Playwright for SVGs with many elements
    if not success:
        print(f"    Trying Playwright renderer...")
        success = render_svg_to_png_playwright(
            svg_path, tmp.name,
            target_width=width,
            target_height=height
        )
        if not success:
            raise RuntimeError(
                f"All renderers (rsvg-convert, cairosvg, playwright) failed for {svg_path}"
            )

    # Calculate scale factor
    if width and svg_w:
        scale = width / svg_w
    elif height and svg_h:
        scale = height / svg_h
    else:
        scale = 1.0

    img = Image.open(tmp.name)
    print(f"    PNG dimensions: {img.size[0]} x {img.size[1]}")
    img.close()

    return tmp.name, scale


def ensure_png(path, svg_height=None, svg_width=None):
    """Return (png_path, is_temp, scale_factor).

    If upscaling fails (element limits), falls back to native resolution.
    """
    if path.lower().endswith(".svg"):
        try:
            png_path, scale = svg_to_png(path, height=svg_height, width=svg_width)
            return png_path, True, scale
        except RuntimeError as e:
            # If we were trying to upscale and it failed, try native resolution
            if svg_height or svg_width:
                print(f"    Upscaling failed, falling back to native resolution...")
                png_path, scale = svg_to_png(path, height=None, width=None)
                return png_path, True, scale
            raise
    return path, False, 1.0


# ─── Font / scale bar helpers ────────────────────────────────────────────────

def _load_font(size=14):
    """Load Basic Sans font with fallback."""
    font_dir = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
    font_path = font_dir / "BasicSans-Regular.otf"
    if font_path.exists():
        return ImageFont.truetype(str(font_path), size)
    try:
        return ImageFont.truetype("Arial", size)
    except (OSError, IOError):
        return ImageFont.load_default()


def _nice_scale_value(bp):
    """Round bp to a nice human-readable scale bar value."""
    nice = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]
    for v in nice:
        if v >= bp * 0.4:
            return v
    return nice[-1]


def draw_dynamic_scale_bar(frame, zoom, ratio, padding_x=20, padding_y=20,
                           target_height=80, text_color=(255, 255, 255)):
    """Draw a scale bar on a video frame that reflects the current zoom level.

    The bar has a fixed approximate visual height; the label snaps to a nice
    round number of bp.
    """
    draw_ctx = ImageDraw.Draw(frame)

    # Physical length represented by target_height viewport pixels
    source_px = target_height / zoom
    bp = source_px / ratio
    nice_bp = _nice_scale_value(bp)

    # Actual visual height for the nice value
    actual_source_px = nice_bp * ratio
    actual_height = int(actual_source_px * zoom)
    actual_height = max(10, min(actual_height, frame.height - 2 * padding_y))

    bar_x = padding_x
    bar_y = padding_y

    # Vertical bar
    draw_ctx.rectangle(
        [bar_x, bar_y, bar_x + 3, bar_y + actual_height], fill=text_color
    )
    # Ticks
    draw_ctx.line([(bar_x - 3, bar_y), (bar_x + 6, bar_y)], fill=text_color, width=1)
    draw_ctx.line(
        [(bar_x - 3, bar_y + actual_height), (bar_x + 6, bar_y + actual_height)],
        fill=text_color,
        width=1,
    )

    label = f"{nice_bp // 1000} Kbp" if nice_bp >= 1000 else f"{nice_bp} bp"
    font = _load_font(14)
    draw_ctx.text((bar_x + 10, bar_y + 2), label, fill=text_color, font=font)


# ─── Adaptive zoom helpers ───────────────────────────────────────────────────

def detect_content_heights(img_array, top_margin, left_margin=0, bg_color="black"):
    """Vectorised column-wise content height detection below *top_margin*.

    Args:
        img_array: numpy array of shape (H, W, 3)
        top_margin: rows above this are ignored (header area)
        left_margin: columns before this are set to 0 (y-axis labels area)
        bg_color: "black" or "white" background

    Returns:
        1D array of heights per column
    """
    h, w, _ = img_array.shape
    if top_margin >= h:
        return np.zeros(w, dtype=np.float64)

    region = img_array[top_margin:]
    is_content = (
        np.any(region != 0, axis=2)
        if bg_color == "black"
        else np.any(region != 255, axis=2)
    )

    row_indices = np.arange(is_content.shape[0]).reshape(-1, 1)
    content_rows = np.where(is_content, row_indices, -1)
    max_rows = content_rows.max(axis=0)
    heights = np.where(max_rows >= 0, max_rows + 1, 0).astype(np.float64)

    # Zero out the left margin (y-axis labels, not actual reads)
    if left_margin > 0:
        heights[:left_margin] = 0

    return heights


def build_height_profile(raw_heights):
    """Smooth and enforce monotonically non-increasing envelope.

    The profile should smoothly decrease (or stay flat) as we pan right,
    reflecting reads sorted by decreasing length.
    """
    profile = raw_heights.copy()

    # Find first non-zero index (content start)
    nonzero_idx = np.nonzero(profile > 0)[0]
    if len(nonzero_idx) == 0:
        return profile  # No content found

    first_content = nonzero_idx[0]

    # Forward-fill zeros within content region (fills spacing gaps between reads)
    last = 0.0
    for i in range(first_content, len(profile)):
        if profile[i] > 0:
            last = profile[i]
        elif last > 0:
            profile[i] = last

    # Enforce monotonically non-increasing starting from first content
    # (Since reads are sorted longest→shortest, heights should only decrease)
    for i in range(first_content + 1, len(profile)):
        if profile[i] > profile[i - 1]:
            profile[i] = profile[i - 1]

    # Fill leading zeros with the first content height (for smooth start)
    if first_content > 0:
        profile[:first_content] = profile[first_content]

    return profile


def _find_content_bounds(profile):
    """Return (first_x, last_x) with non-zero profile values."""
    nonzero = np.nonzero(profile > 0)[0]
    if len(nonzero) == 0:
        return 0, len(profile)
    return int(nonzero[0]), int(nonzero[-1]) + 1


# ─── Adaptive zoom panning ───────────────────────────────────────────────────

def _smooth_profile(profile, window=100):
    """Apply exponential moving average for smoother zoom transitions."""
    smoothed = np.zeros_like(profile)
    alpha = 2 / (window + 1)
    smoothed[0] = profile[0]
    for i in range(1, len(profile)):
        smoothed[i] = alpha * profile[i] + (1 - alpha) * smoothed[i - 1]
    return smoothed


# ─── Strip-based pre-rendering ────────────────────────────────────────────────

def compute_strip_zoom_levels(profile, num_strips, available_height, max_zoom, content_padding=10):
    """Divide SVG into strips, compute max zoom needed for each.

    Args:
        profile: Content height profile (1D array, height per x-position)
        num_strips: Number of strips to divide into
        available_height: Viewport height available for content
        max_zoom: Maximum zoom factor
        content_padding: Padding around content

    Returns:
        List of (x_start, x_end, zoom_level) tuples
    """
    total_width = len(profile)
    strip_width = max(1, total_width // num_strips)
    strips = []

    for i in range(num_strips):
        x_start = i * strip_width
        x_end = min((i + 1) * strip_width, total_width)
        if x_start >= total_width:
            break

        # Max content height in this strip (use max to ensure we render enough)
        strip_profile = profile[x_start:x_end]
        if len(strip_profile) == 0:
            continue
        max_content = max(strip_profile.max(), 20)  # minimum content height

        # Zoom needed to fill available_height
        zoom = min(available_height / (max_content + content_padding), max_zoom)
        strips.append((x_start, x_end, zoom))

    return strips


def prerender_strips(svg_path, strips, svg_h, top_margin, available_height, output_dir, use_playwright=False):
    """Render each strip at fixed output height using rsvg-convert or Playwright.

    Each strip is rendered at `available_height` tall, with width proportional
    to the strip's zoom level. This means strips at higher zoom are wider.

    Args:
        svg_path: Path to SVG file
        strips: List of (x_start, x_end, zoom) tuples
        svg_h: SVG height in pixels
        top_margin: Top margin to skip (header area)
        available_height: Target output height for all strips
        output_dir: Directory to store strip PNGs
        use_playwright: If True, use Playwright instead of rsvg-convert

    Returns:
        List of (path, strip_x_start, strip_x_end, zoom) tuples
    """
    os.makedirs(output_dir, exist_ok=True)
    strip_info = []
    env = os.environ.copy()
    env["RSVG_MAX_LOADED_ELEMENTS"] = "100000000"

    content_height = svg_h - top_margin
    rsvg_failed = False

    for i, (x_start, x_end, zoom) in enumerate(strips):
        strip_w = x_end - x_start
        out_path = os.path.join(output_dir, f"strip_{i:03d}.png")

        # Output dimensions:
        # - Height = available_height (fixed)
        # - Width = strip_w * zoom (varies by zoom)
        page_h = available_height
        page_w = int(strip_w * zoom)

        success = False

        # Try rsvg-convert first (unless it already failed or Playwright forced)
        if not use_playwright and not rsvg_failed:
            # The SVG is scaled uniformly by zoom factor, then offset
            # to bring the strip region into view
            left_offset = -x_start * zoom
            top_offset = -top_margin * zoom

            cmd = [
                "rsvg-convert",
                "--unlimited",
                f"--zoom={zoom}",
                f"--left={left_offset}",
                f"--top={top_offset}",
                f"--page-width={page_w}",
                f"--page-height={page_h}",
                "-o", out_path,
                svg_path,
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
                success = True
            except subprocess.CalledProcessError as e:
                if "cannot load more than" in e.stderr:
                    # Element limit exceeded, switch to Playwright for all remaining strips
                    if not rsvg_failed:
                        print(f"  SVG has too many elements, switching to Playwright renderer...")
                    rsvg_failed = True
                else:
                    print(f"  Warning: rsvg-convert failed for strip {i}: {e.stderr}")

        # Fall back to Playwright
        if not success:
            success = render_svg_region_playwright(
                svg_path, out_path,
                x_offset=x_start,
                y_offset=top_margin,
                width=strip_w,
                height=content_height,
                scale=zoom
            )

        if success and os.path.exists(out_path):
            strip_info.append((out_path, x_start, x_end, zoom))
        else:
            print(f"  Warning: Failed to render strip {i}")
            strip_info.append((None, x_start, x_end, zoom))

    return strip_info


def load_strip_images(strip_info):
    """Load pre-rendered strip PNGs into memory.

    Returns:
        List of (numpy_array, x_start, x_end, zoom) tuples
    """
    images = []
    for path, x_start, x_end, zoom in strip_info:
        if path and os.path.exists(path):
            img = Image.open(path).convert("RGB")
            images.append((np.array(img), x_start, x_end, zoom))
        else:
            images.append((None, x_start, x_end, zoom))
    return images


def compose_frame_from_strips(
    strip_data, x_pos, panning_width, available_height, bg_color
):
    """Compose a frame by extracting from pre-rendered strips.

    The strips are pre-rendered at the final output height, so we just need
    to extract the horizontal region and potentially resize width only.

    Args:
        strip_data: List of (numpy_array, x_start, x_end, zoom) tuples
        x_pos: Current pan x position (in source SVG coordinates)
        panning_width: Width of panning viewport in output pixels
        available_height: Height available for content in output pixels (should match strip height)
        bg_color: Background color tuple (R, G, B)

    Returns:
        PIL Image of the composed panning content
    """
    x_end_pos = x_pos + panning_width

    # Create output frame
    frame = Image.new("RGB", (panning_width, available_height), bg_color)
    out_x = 0

    for strip_arr, strip_start, strip_end, strip_zoom in strip_data:
        if strip_arr is None:
            continue

        # Check if this strip overlaps with visible region
        if strip_end <= x_pos or strip_start >= x_end_pos:
            continue

        # Calculate overlap in source coordinates
        overlap_start = max(x_pos, strip_start)
        overlap_end = min(x_end_pos, strip_end)
        overlap_width = overlap_end - overlap_start

        if overlap_width <= 0:
            continue

        # Position within the strip image (in strip's zoomed coordinates)
        strip_local_x = int((overlap_start - strip_start) * strip_zoom)
        strip_local_w = int(overlap_width * strip_zoom)

        strip_h, strip_w_px, _ = strip_arr.shape

        # Clamp to strip bounds
        strip_local_x = max(0, min(strip_local_x, strip_w_px - 1))
        strip_local_w = min(strip_local_w, strip_w_px - strip_local_x)

        if strip_local_w <= 0:
            continue

        # Extract the region (full height, cropped horizontally)
        crop_arr = strip_arr[:, strip_local_x:strip_local_x + strip_local_w]
        crop_img = Image.fromarray(crop_arr)

        # Output width for this portion (1:1 in source coordinates)
        out_w = int(overlap_width)

        # Resize width to match 1:1 source mapping (strip was rendered at zoom, so compress back)
        # Height is already correct (available_height)
        if crop_img.width != out_w or crop_img.height != available_height:
            crop_resized = crop_img.resize((out_w, available_height), Image.Resampling.LANCZOS)
        else:
            crop_resized = crop_img

        # Paste into frame
        frame.paste(crop_resized, (out_x, 0))
        out_x += out_w

    return frame


def _render_svg_frame(svg_path, x_pos, y_pos, zoom, page_width, page_height):
    """Render a region of an SVG at given zoom level.

    Uses rsvg-convert with --left/--top offsets to render a specific region.

    Args:
        svg_path: Path to SVG file
        x_pos: X position in SVG coordinates (left edge of viewport)
        y_pos: Y position in SVG coordinates (top edge of viewport)
        zoom: Zoom factor (uniform scaling)
        page_width: Output width in pixels
        page_height: Output height in pixels

    Returns:
        PIL Image of the rendered frame
    """
    # Calculate offsets (negative to bring content into view)
    left_offset = -x_pos * zoom
    top_offset = -y_pos * zoom

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = tmp.name

    cmd = [
        "rsvg-convert",
        f"--zoom={zoom}",
        f"--left={left_offset}",
        f"--top={top_offset}",
        f"--page-width={page_width}",
        f"--page-height={page_height}",
        "-o", tmp_path,
        svg_path,
    ]

    env = os.environ.copy()
    env["RSVG_MAX_LOADED_ELEMENTS"] = "10000000"

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        img = Image.open(tmp_path).convert("RGB")
        return img
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def create_adaptive_horizontal_panning_svg(
    svg_path,
    output_path,
    duration,
    fps,
    viewport_width,
    viewport_height,
    legend_path,
    background,
    top_margin,
    left_margin,
    ratio,
    max_zoom,
    scale_bar_padding,
    scale_bar_width=100,
    zoom_smoothing=200,
    num_strips=50,
):
    """Horizontal panning with adaptive zoom using strip-based pre-rendering.

    Pre-renders the SVG in strips at appropriate zoom levels, then composes
    frames from the pre-rendered strips. Much faster than per-frame SVG rendering.
    """
    print("Creating adaptive-zoom horizontal panning animation (strip-based)")
    print(f"  Input:    {svg_path}")
    print(f"  Output:   {output_path}")
    print(f"  Duration: {duration:.1f}s  FPS: {fps}")
    print(f"  Viewport: {viewport_width} x {viewport_height}")

    # Get SVG dimensions
    svg_w, svg_h = get_svg_dimensions(svg_path)
    if svg_w is None or svg_h is None:
        raise ValueError("Could not parse SVG dimensions")
    print(f"  SVG:      {svg_w:.0f} x {svg_h:.0f}")

    # Render a low-res version for content height detection
    # Limit to 32000px max dimension to stay under librsvg limit
    max_render_dim = 32000
    preview_scale = 1.0
    if svg_w > max_render_dim or svg_h > max_render_dim:
        preview_scale = min(max_render_dim / svg_w, max_render_dim / svg_h)
        preview_w = int(svg_w * preview_scale)
        preview_h = int(svg_h * preview_scale)
        print(f"  Rendering preview at {preview_scale:.3f}x ({preview_w} x {preview_h}) for content detection...")
    else:
        preview_w, preview_h = int(svg_w), int(svg_h)
        print("  Rendering preview for content detection...")

    # Try to render preview for content detection
    profile = None
    preview_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            preview_path = tmp.name
        cmd = ["rsvg-convert", "--unlimited", "-o", preview_path]
        if preview_scale < 1.0:
            cmd.extend(["-w", str(preview_w)])
        cmd.append(svg_path)
        env = os.environ.copy()
        env["RSVG_MAX_LOADED_ELEMENTS"] = "100000000"
        subprocess.run(cmd, check=True, capture_output=True, env=env)
        preview_img = Image.open(preview_path)
        preview_array = np.array(preview_img)
        img_h, img_w = preview_array.shape[:2]

        # Content height profile
        # Scale margins for preview if needed
        scaled_top_margin = int(top_margin * preview_scale)
        scaled_left_margin = int(left_margin * preview_scale)
        print(f"  Detecting content heights (top_margin={top_margin}, left_margin={left_margin})…")
        raw_heights = detect_content_heights(preview_array, scaled_top_margin, scaled_left_margin, background)
        profile = build_height_profile(raw_heights)

        # Scale profile back to original SVG coordinates
        if preview_scale < 1.0:
            # Scale height values (content heights) back to SVG coordinates
            profile = profile / preview_scale
            # Resample profile to match original SVG width
            original_width = int(svg_w - left_margin)
            x_new = np.linspace(0, len(profile) - 1, original_width)
            profile = np.interp(x_new, np.arange(len(profile)), profile)

        # Smooth the profile for gradual zoom transitions
        if zoom_smoothing > 0:
            profile = _smooth_profile(profile, window=zoom_smoothing)
            print(f"  Applied zoom smoothing (window={zoom_smoothing})")

        max_c = profile.max()
        nonzero = profile[profile > 0]
        min_c = nonzero.min() if len(nonzero) else 1
        print(f"  Content height range: {max_c:.0f} – {min_c:.0f} px (source)")

    except subprocess.CalledProcessError as e:
        # Preview rendering failed (likely too many elements)
        # Fall back to uniform zoom across all strips
        print(f"  Warning: Preview rendering failed (SVG may have too many elements)")
        print(f"  Using uniform zoom based on SVG height...")

        # Estimate content height from SVG height minus margins
        content_height = svg_h - top_margin
        # Create uniform profile with SVG height
        original_width = int(svg_w - left_margin)
        profile = np.full(original_width, content_height, dtype=float)
        # Set img_w/img_h from SVG dimensions since we didn't get preview
        img_w = int(svg_w)
        img_h = int(svg_h)

        max_c = content_height
        min_c = content_height
        print(f"  Using uniform content height: {max_c:.0f} px")

    finally:
        if preview_path and os.path.exists(preview_path):
            os.unlink(preview_path)

    content_start, content_end = _find_content_bounds(profile)
    print(f"  Content x-range: {content_start} – {content_end}")

    # Layout
    panning_width = viewport_width - scale_bar_width
    if panning_width < 400:
        scale_bar_width = max(60, viewport_width // 10)
        panning_width = viewport_width - scale_bar_width
    print(f"  Layout: scale_bar={scale_bar_width}px + panning={panning_width}px")

    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    bg_tuple = (0, 0, 0) if background == "black" else (255, 255, 255)
    content_padding = 10
    min_content_h = 20
    reads_top_y = 20
    available_height = viewport_height - reads_top_y

    # Compute strip zoom levels
    strips = compute_strip_zoom_levels(
        profile, num_strips, available_height, max_zoom, content_padding
    )
    print(f"  Strips:   {len(strips)} (zoom range: {strips[0][2]:.2f}x – {strips[-1][2]:.2f}x)")

    # Pre-render strips
    strip_dir = tempfile.mkdtemp(prefix="panning_strips_")
    try:
        print(f"  Pre-rendering {len(strips)} strips...")
        strip_info = prerender_strips(svg_path, strips, int(svg_h), top_margin, available_height, strip_dir)
        print(f"  Loading strips into memory...")
        strip_data = load_strip_images(strip_info)

        # Legend
        legend_img = None
        legend_array = None
        legend_h = 0
        if legend_path and os.path.exists(legend_path):
            if legend_path.lower().endswith(".svg"):
                legend_png, _ = svg_to_png(legend_path, width=viewport_width)
                legend_img = Image.open(legend_png).convert("RGB")
                os.unlink(legend_png)
            else:
                legend_img = Image.open(legend_path).convert("RGB")
            lw, lh = legend_img.size
            if lw != viewport_width:
                s = viewport_width / lw
                legend_img = legend_img.resize(
                    (viewport_width, int(lh * s)), Image.Resampling.LANCZOS
                )
            legend_h = legend_img.size[1]
            if legend_h % 2 != 0:
                legend_h += 1
                padded = Image.new("RGB", (viewport_width, legend_h), background)
                padded.paste(legend_img, (0, 0))
                legend_img = padded
            legend_array = np.array(legend_img)
            print(f"  Legend:   {legend_img.size[0]} x {legend_img.size[1]}")

        # Ensure even dimensions for H.264
        total_h = viewport_height + legend_h
        if viewport_width % 2 != 0:
            viewport_width += 1
            panning_width = viewport_width - scale_bar_width
        if total_h % 2 != 0:
            total_h += 1
            legend_h = total_h - viewport_height

        total_frames = int(duration * fps)
        pan_start = max(0, content_start - 10)
        pan_end = content_end
        pan_distance = max(1, pan_end - pan_start)

        print(f"  Frames:   {total_frames}")
        print(f"  Max zoom: {max_zoom}")

        # Start ffmpeg encoder
        ffmpeg_cmd = [
            "ffmpeg", "-y",
            "-f", "rawvideo",
            "-pix_fmt", "rgb24",
            "-s", f"{viewport_width}x{total_h}",
            "-r", str(fps),
            "-i", "-",
            "-c:v", "libx264",
            "-preset", "medium",
            "-crf", "18",
            "-pix_fmt", "yuv420p",
            output_path,
        ]
        process = subprocess.Popen(ffmpeg_cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        print(f"\nRendering {total_frames} frames …")
        last_pct = -1

        for fidx in range(total_frames):
            t = fidx / max(1, total_frames - 1)

            # Source x position (left edge of panning viewport)
            x_pos = int(pan_start + t * pan_distance)
            x_pos = min(x_pos, img_w - panning_width)
            x_pos = max(0, x_pos)

            # Content height at current x position
            ch = max(min_content_h, profile[min(x_pos, len(profile) - 1)])

            # Vertical zoom to fit content in available viewport height
            v_zoom = min(available_height / (ch + content_padding), max_zoom)

            # Compose frame from pre-rendered strips
            panning_frame = compose_frame_from_strips(
                strip_data, x_pos, panning_width, available_height, bg_tuple
            )

            # Create full frame with scale bar area
            frame = Image.new("RGB", (viewport_width, viewport_height), bg_tuple)

            # Paste panning content at (scale_bar_width, reads_top_y)
            frame.paste(panning_frame, (scale_bar_width, reads_top_y))

            # Draw scale bar
            if ratio > 0:
                draw_ctx = ImageDraw.Draw(frame)
                bar_height = min(150, available_height - 20)

                source_px = bar_height / v_zoom
                bp = source_px / ratio
                nice_bp = _nice_scale_value(bp)

                actual_source_px = nice_bp * ratio
                actual_height = int(actual_source_px * v_zoom)
                actual_height = max(10, min(actual_height, available_height - 20))

                bar_x = scale_bar_width - 12
                bar_y = reads_top_y

                draw_ctx.rectangle(
                    [bar_x, bar_y, bar_x + 4, bar_y + actual_height], fill=text_color
                )
                draw_ctx.line([(bar_x - 4, bar_y), (bar_x + 8, bar_y)], fill=text_color, width=2)
                draw_ctx.line(
                    [(bar_x - 4, bar_y + actual_height), (bar_x + 8, bar_y + actual_height)],
                    fill=text_color, width=2,
                )

                label = f"{nice_bp // 1000} Kbp" if nice_bp >= 1000 else f"{nice_bp} bp"
                font = _load_font(36)

                txt_img = Image.new("RGBA", (200, 50), (0, 0, 0, 0))
                txt_draw = ImageDraw.Draw(txt_img)
                bbox = txt_draw.textbbox((0, 0), label, font=font)
                txt_w = bbox[2] - bbox[0]
                txt_h = bbox[3] - bbox[1]
                txt_draw.text(((200 - txt_w) // 2, (50 - txt_h) // 2), label, fill=text_color + (255,), font=font)
                txt_rotated = txt_img.rotate(90, expand=True)

                bar_center_y = bar_y + actual_height // 2
                label_x = bar_x - 55
                label_y = bar_center_y - txt_rotated.height // 2
                label_y = max(0, label_y)
                frame.paste(txt_rotated, (label_x, label_y), txt_rotated)

            # Assemble output frame
            frame_arr = np.array(frame)
            if legend_array is not None:
                out = np.empty((total_h, viewport_width, 3), dtype=np.uint8)
                out[:viewport_height] = frame_arr
                out[viewport_height:viewport_height + legend_array.shape[0]] = legend_array
                if viewport_height + legend_array.shape[0] < total_h:
                    out[viewport_height + legend_array.shape[0]:] = bg_tuple
            else:
                out = frame_arr

            process.stdin.write(out.tobytes())

            pct = int(100 * fidx / total_frames)
            if pct != last_pct and pct % 5 == 0:
                last_pct = pct
                print(f"  {pct}%  (frame {fidx}/{total_frames})", flush=True)

        process.stdin.close()
        process.wait()

        out_size = os.path.getsize(output_path) / (1024 * 1024)
        print(f"\n  Saved: {output_path}")
        print(f"  Size:  {out_size:.1f} MB")

    finally:
        # Clean up strip directory
        if os.path.exists(strip_dir):
            shutil.rmtree(strip_dir)
        # Clean up Playwright if used
        _close_playwright()


def create_adaptive_horizontal_panning(
    input_path,
    output_path,
    duration,
    fps,
    viewport_width,
    viewport_height,
    legend_path,
    background,
    top_margin,
    left_margin,
    ratio,
    max_zoom,
    scale_bar_padding,
    scale_bar_width=100,
    zoom_smoothing=200,
):
    """Horizontal panning with adaptive vertical-only zoom.

    Key features:
    - Vertical-only zoom: keeps constant number of reads visible (horizontal extent fixed)
    - Smooth zoom transitions via exponential moving average
    - Scale bar in dedicated left margin (not overlapping content)
    - Legend composited below panning area
    """

    print("Creating adaptive-zoom horizontal panning animation")
    print(f"  Input:    {input_path}")
    print(f"  Output:   {output_path}")
    print(f"  Duration: {duration:.1f}s  FPS: {fps}")
    print(f"  Viewport: {viewport_width} x {viewport_height}")

    # Load image into numpy array
    img = Image.open(input_path)
    img_array = np.array(img)
    img_h, img_w, _ = img_array.shape
    print(f"  Image:    {img_w} x {img_h}")

    # Content height profile
    print(f"  Detecting content heights (top_margin={top_margin}, left_margin={left_margin})…")
    raw_heights = detect_content_heights(img_array, top_margin, left_margin, background)
    profile = build_height_profile(raw_heights)

    # Smooth the profile for gradual zoom transitions
    if zoom_smoothing > 0:
        profile = _smooth_profile(profile, window=zoom_smoothing)
        print(f"  Applied zoom smoothing (window={zoom_smoothing})")

    max_c = profile.max()
    nonzero = profile[profile > 0]
    min_c = nonzero.min() if len(nonzero) else 1
    print(f"  Content height range: {max_c:.0f} – {min_c:.0f} px (source)")

    content_start, content_end = _find_content_bounds(profile)
    print(f"  Content x-range: {content_start} – {content_end}")

    # Allocate layout:
    # [scale_bar_width] [panning_area] = viewport_width
    panning_width = viewport_width - scale_bar_width
    if panning_width < 400:
        print(f"  Warning: panning area too narrow, reducing scale bar width")
        scale_bar_width = max(60, viewport_width // 10)
        panning_width = viewport_width - scale_bar_width
    print(f"  Layout: scale_bar={scale_bar_width}px + panning={panning_width}px")

    # Legend
    legend_img = None
    legend_array = None
    legend_h = 0
    if legend_path and os.path.exists(legend_path):
        legend_img = Image.open(legend_path).convert("RGB")
        lw, lh = legend_img.size
        if lw != viewport_width:
            s = viewport_width / lw
            legend_img = legend_img.resize(
                (viewport_width, int(lh * s)), Image.Resampling.LANCZOS
            )
        legend_h = legend_img.size[1]
        if legend_h % 2 != 0:
            legend_h += 1
            padded = Image.new("RGB", (viewport_width, legend_h), background)
            padded.paste(legend_img, (0, 0))
            legend_img = padded
        legend_array = np.array(legend_img)
        print(f"  Legend:   {legend_img.size[0]} x {legend_img.size[1]}")

    # Ensure even dimensions for H.264
    total_h = viewport_height + legend_h
    if viewport_width % 2 != 0:
        viewport_width += 1
        panning_width = viewport_width - scale_bar_width
    if total_h % 2 != 0:
        total_h += 1
        legend_h = total_h - viewport_height

    total_frames = int(duration * fps)
    pan_start = max(0, content_start - 10)
    pan_end = content_end
    pan_distance = max(1, pan_end - pan_start)

    print(f"  Frames:   {total_frames}")
    print(f"  Max zoom: {max_zoom}")

    text_color = (255, 255, 255) if background == "black" else (0, 0, 0)
    bg_tuple = (0, 0, 0) if background == "black" else (255, 255, 255)
    content_padding = 10
    min_content_h = 20

    # Determine fixed horizontal crop width in source pixels
    # This ensures constant number of reads visible throughout
    fixed_crop_w = panning_width  # 1:1 mapping for horizontal

    # Start ffmpeg encoder
    ffmpeg_cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", f"{viewport_width}x{total_h}",
        "-r", str(fps),
        "-i", "-",
        "-c:v", "libx264",
        "-preset", "medium",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        output_path,
    ]
    process = subprocess.Popen(ffmpeg_cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    print(f"\nRendering {total_frames} frames …")
    last_pct = -1

    # Fixed y-position for read tops in viewport (small padding from top)
    reads_top_y = 20

    for fidx in range(total_frames):
        t = fidx / max(1, total_frames - 1)

        # Source x position (left edge of panning viewport)
        x_pos = int(pan_start + t * pan_distance)
        x_pos = min(x_pos, img_w - fixed_crop_w)
        x_pos = max(0, x_pos)

        # Content height at current x position (below top_margin)
        ch = max(min_content_h, profile[min(x_pos, len(profile) - 1)])

        # Vertical zoom: scale content to fit in available viewport height
        available_height = viewport_height - reads_top_y
        v_zoom = min(available_height / (ch + content_padding), max_zoom)

        # Crop region in source - start at top_margin (where reads begin)
        # This anchors the TOP of reads at a fixed position
        crop_y0 = top_margin
        crop_y1 = int(top_margin + available_height / v_zoom)
        crop_y1 = min(crop_y1, img_h)

        crop_x0 = x_pos
        crop_x1 = min(crop_x0 + fixed_crop_w, img_w)

        crop = img_array[crop_y0:crop_y1, crop_x0:crop_x1]

        # Resize: horizontal stays 1:1, vertical is scaled
        if crop.shape[0] > 0 and crop.shape[1] > 0:
            new_w = crop.shape[1]  # Keep horizontal pixels
            new_h = available_height
            panning_frame = Image.fromarray(crop).resize(
                (new_w, new_h), Image.Resampling.BILINEAR  # BILINEAR is faster than LANCZOS
            )
            # If crop was narrower than panning_width, pad with background
            if new_w < panning_width:
                padded = Image.new("RGB", (panning_width, available_height), bg_tuple)
                padded.paste(panning_frame, (0, 0))
                panning_frame = padded
        else:
            panning_frame = Image.new("RGB", (panning_width, available_height), bg_tuple)

        # Create frame with scale bar directly adjacent to panning content
        frame = Image.new("RGB", (viewport_width, viewport_height), bg_tuple)

        # Paste panning content at (scale_bar_width, reads_top_y)
        frame.paste(panning_frame, (scale_bar_width, reads_top_y))

        # Draw scale bar flush against the panning content, aligned with reads top
        if ratio > 0:
            draw_ctx = ImageDraw.Draw(frame)
            bar_height = min(150, available_height - 20)

            # Physical length represented by bar_height viewport pixels
            source_px = bar_height / v_zoom
            bp = source_px / ratio
            nice_bp = _nice_scale_value(bp)

            # Actual visual height for the nice value
            actual_source_px = nice_bp * ratio
            actual_height = int(actual_source_px * v_zoom)
            actual_height = max(10, min(actual_height, available_height - 20))

            # Position: flush against right edge of scale bar area, TOP aligned with reads
            bar_x = scale_bar_width - 12
            bar_y = reads_top_y  # Align with top of reads

            # Vertical bar (thicker)
            draw_ctx.rectangle(
                [bar_x, bar_y, bar_x + 4, bar_y + actual_height], fill=text_color
            )
            # Ticks
            draw_ctx.line([(bar_x - 4, bar_y), (bar_x + 8, bar_y)], fill=text_color, width=2)
            draw_ctx.line(
                [(bar_x - 4, bar_y + actual_height), (bar_x + 8, bar_y + actual_height)],
                fill=text_color, width=2,
            )

            # Vertical label (rotated 90 degrees, large font)
            label = f"{nice_bp // 1000} Kbp" if nice_bp >= 1000 else f"{nice_bp} bp"
            font = _load_font(36)

            # Create temporary image for rotated text
            txt_img = Image.new("RGBA", (200, 50), (0, 0, 0, 0))
            txt_draw = ImageDraw.Draw(txt_img)
            # Get text bounding box for centering
            bbox = txt_draw.textbbox((0, 0), label, font=font)
            txt_w = bbox[2] - bbox[0]
            txt_h = bbox[3] - bbox[1]
            # Draw text centered in the temporary image
            txt_draw.text(((200 - txt_w) // 2, (50 - txt_h) // 2), label, fill=text_color + (255,), font=font)
            txt_rotated = txt_img.rotate(90, expand=True)

            # Position: center of label aligns with center of scale bar
            # But clamp so top of label doesn't go above frame top
            bar_center_y = bar_y + actual_height // 2
            label_x = bar_x - 55
            label_y = bar_center_y - txt_rotated.height // 2
            label_y = max(0, label_y)  # Don't let label go above top of frame
            frame.paste(txt_rotated, (label_x, label_y), txt_rotated)

        # Assemble output frame (panning + legend)
        frame_arr = np.array(frame)
        if legend_array is not None:
            out = np.empty((total_h, viewport_width, 3), dtype=np.uint8)
            out[:viewport_height] = frame_arr
            out[viewport_height:viewport_height + legend_array.shape[0]] = legend_array
            if viewport_height + legend_array.shape[0] < total_h:
                out[viewport_height + legend_array.shape[0]:] = bg_tuple
        else:
            out = frame_arr

        process.stdin.write(out.tobytes())

        pct = int(100 * fidx / total_frames)
        if pct != last_pct and pct % 5 == 0:
            last_pct = pct
            print(f"  {pct}%  (frame {fidx}/{total_frames})", flush=True)

    process.stdin.close()
    process.wait()

    if process.returncode != 0:
        stderr = process.stderr.read().decode()
        print(f"  ffmpeg error:\n{stderr[-500:]}")
        raise RuntimeError("ffmpeg encoding failed")

    fsize = os.path.getsize(output_path)
    print(f"\n  Saved: {output_path}")
    print(f"  Size:  {fsize / 1024 / 1024:.1f} MB")
    return output_path


# ─── Fixed-zoom panning (original) ──────────────────────────────────────────

def create_horizontal_panning(
    input_path, output_path, duration, fps, crop_ratio,
    viewport_width, viewport_height,
    scale_bar_path, legend_path, background, scale_bar_padding,
):
    """Seamless horizontal panning animation (fixed zoom)."""

    print("Creating horizontal panning animation (fixed zoom)")
    print(f"  Input:    {input_path}")
    print(f"  Output:   {output_path}")
    print(f"  Duration: {duration:.1f}s  FPS: {fps}")

    img = Image.open(input_path)
    orig_width, orig_height = img.size
    print(f"  Original size: {orig_width} x {orig_height}")

    # Crop to top portion
    crop_height = int(orig_height * crop_ratio)
    img_cropped = img.crop((0, 0, orig_width, crop_height))
    print(f"  Cropped to top {crop_ratio*100:.0f}%: {orig_width} x {crop_height}")

    # Scale to viewport height
    scale_factor = viewport_height / crop_height
    scaled_width = int(orig_width * scale_factor)
    img_scaled = img_cropped.resize(
        (scaled_width, viewport_height), Image.Resampling.LANCZOS
    )
    print(f"  Scaled to viewport height: {scaled_width} x {viewport_height}")

    # Seamless tile
    tiled_width = scaled_width * 2
    img_tiled = Image.new("RGB", (tiled_width, viewport_height), background)
    img_tiled.paste(img_scaled, (0, 0))
    img_tiled.paste(img_scaled, (scaled_width, 0))
    print(f"  Tiled for seamless loop: {tiled_width} x {viewport_height}")

    temp_files = []
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tiled_path = tmp.name
        temp_files.append(tiled_path)
        img_tiled.save(tiled_path, "PNG")

    pan_distance = scaled_width
    total_frames = int(duration * fps)
    px_per_frame = pan_distance / total_frames

    print(f"  Pan distance: {pan_distance}px over {total_frames} frames")
    print(f"  Speed: {px_per_frame:.2f} px/frame")

    # ffmpeg filter chain
    filter_parts = []
    inputs = ["-loop", "1", "-i", tiled_path]
    input_idx = 0

    crop_filter = (
        f"[{input_idx}:v]crop={viewport_width}:{viewport_height}:"
        f"'mod(t*{pan_distance}/{duration},iw-{viewport_width})':0[panning]"
    )
    filter_parts.append(crop_filter)
    last_stream = "panning"
    input_idx += 1

    if scale_bar_path and os.path.exists(scale_bar_path):
        inputs.extend(["-i", scale_bar_path])
        overlay_filter = (
            f"[{last_stream}][{input_idx}:v]overlay="
            f"{scale_bar_padding}:{scale_bar_padding}[with_scale]"
        )
        filter_parts.append(overlay_filter)
        last_stream = "with_scale"
        input_idx += 1
        print(f"  Adding scale bar: {scale_bar_path}")

    if legend_path and os.path.exists(legend_path):
        with tempfile.NamedTemporaryFile(suffix=".mp4", delete=False) as tmp:
            panning_video_path = tmp.name
            temp_files.append(panning_video_path)

        filter_complex = ";".join(filter_parts)
        run_ffmpeg([
            "ffmpeg", "-y", *inputs,
            "-filter_complex", filter_complex,
            "-map", f"[{last_stream}]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            panning_video_path,
        ])

        legend_img = Image.open(legend_path)
        lw, lh = legend_img.size
        if lw != viewport_width:
            s = viewport_width / lw
            legend_img = legend_img.resize(
                (viewport_width, int(lh * s)), Image.Resampling.LANCZOS
            )

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            legend_scaled_path = tmp.name
            temp_files.append(legend_scaled_path)
            legend_img.save(legend_scaled_path, "PNG")

        run_ffmpeg([
            "ffmpeg", "-y",
            "-i", panning_video_path,
            "-loop", "1", "-i", legend_scaled_path,
            "-filter_complex", "[0:v][1:v]vstack=inputs=2[out]",
            "-map", "[out]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            output_path,
        ])
    else:
        filter_complex = ";".join(filter_parts)
        run_ffmpeg([
            "ffmpeg", "-y", *inputs,
            "-filter_complex", filter_complex,
            "-map", f"[{last_stream}]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            output_path,
        ])

    for f in temp_files:
        if os.path.exists(f):
            os.unlink(f)

    fsize = os.path.getsize(output_path)
    print(f"  Saved: {output_path}")
    print(f"  File size: {fsize / 1024 / 1024:.1f} MB")
    return output_path


def create_vertical_panning(
    input_path, output_path, duration, fps,
    viewport_width, viewport_height,
    scale_bar_path, legend_path, background, scale_bar_padding,
):
    """Seamless vertical panning animation (fixed zoom)."""

    print("Creating vertical panning animation (fixed zoom)")
    print(f"  Input:    {input_path}")
    print(f"  Output:   {output_path}")
    print(f"  Duration: {duration:.1f}s  FPS: {fps}")

    img = Image.open(input_path)
    orig_width, orig_height = img.size
    print(f"  Original size: {orig_width} x {orig_height}")

    scale_factor = viewport_width / orig_width
    scaled_height = int(orig_height * scale_factor)
    img_scaled = img.resize(
        (viewport_width, scaled_height), Image.Resampling.LANCZOS
    )
    print(f"  Scaled to viewport width: {viewport_width} x {scaled_height}")

    min_tiled_height = viewport_height + scaled_height
    num_tiles = max(2, (min_tiled_height // scaled_height) + 1)
    tiled_height = scaled_height * num_tiles
    img_tiled = Image.new("RGB", (viewport_width, tiled_height), background)
    for i in range(num_tiles):
        img_tiled.paste(img_scaled, (0, i * scaled_height))
    print(f"  Tiled: {viewport_width} x {tiled_height} ({num_tiles} tiles)")

    temp_files = []
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tiled_path = tmp.name
        temp_files.append(tiled_path)
        img_tiled.save(tiled_path, "PNG")

    pan_distance = scaled_height
    total_frames = int(duration * fps)
    px_per_frame = pan_distance / total_frames

    print(f"  Pan distance: {pan_distance}px over {total_frames} frames")
    print(f"  Speed: {px_per_frame:.2f} px/frame")

    filter_parts = []
    inputs = ["-loop", "1", "-i", tiled_path]
    input_idx = 0

    crop_filter = (
        f"[{input_idx}:v]crop={viewport_width}:{viewport_height}:"
        f"0:'mod(t*{pan_distance}/{duration},{pan_distance})'[panning]"
    )
    filter_parts.append(crop_filter)
    last_stream = "panning"
    input_idx += 1

    if scale_bar_path and os.path.exists(scale_bar_path):
        inputs.extend(["-i", scale_bar_path])
        overlay_filter = (
            f"[{last_stream}][{input_idx}:v]overlay="
            f"{scale_bar_padding}:{scale_bar_padding}[with_scale]"
        )
        filter_parts.append(overlay_filter)
        last_stream = "with_scale"
        input_idx += 1
        print(f"  Adding scale bar: {scale_bar_path}")

    if legend_path and os.path.exists(legend_path):
        with tempfile.NamedTemporaryFile(suffix=".mp4", delete=False) as tmp:
            panning_video_path = tmp.name
            temp_files.append(panning_video_path)

        filter_complex = ";".join(filter_parts)
        run_ffmpeg([
            "ffmpeg", "-y", *inputs,
            "-filter_complex", filter_complex,
            "-map", f"[{last_stream}]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            panning_video_path,
        ])

        legend_img = Image.open(legend_path)
        lw, lh = legend_img.size
        if lh != viewport_height:
            s = viewport_height / lh
            legend_img = legend_img.resize(
                (int(lw * s), viewport_height), Image.Resampling.LANCZOS
            )

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            legend_scaled_path = tmp.name
            temp_files.append(legend_scaled_path)
            legend_img.save(legend_scaled_path, "PNG")

        run_ffmpeg([
            "ffmpeg", "-y",
            "-i", panning_video_path,
            "-loop", "1", "-i", legend_scaled_path,
            "-filter_complex", "[0:v][1:v]hstack=inputs=2[out]",
            "-map", "[out]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            output_path,
        ])
    else:
        filter_complex = ";".join(filter_parts)
        run_ffmpeg([
            "ffmpeg", "-y", *inputs,
            "-filter_complex", filter_complex,
            "-map", f"[{last_stream}]",
            "-t", str(int(duration)),
            "-r", str(fps),
            "-c:v", "libx264", "-preset", "medium", "-crf", "18",
            "-pix_fmt", "yuv420p",
            output_path,
        ])

    for f in temp_files:
        if os.path.exists(f):
            os.unlink(f)

    fsize = os.path.getsize(output_path)
    print(f"  Saved: {output_path}")
    print(f"  File size: {fsize / 1024 / 1024:.1f} MB")
    return output_path


# ─── Utilities ───────────────────────────────────────────────────────────────

def run_ffmpeg(cmd):
    """Run ffmpeg with error handling."""
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"  ffmpeg error: {e.stderr[-500:]}")
        raise


# ─── CLI ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Create panning animations from wide/tall images",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Input / output
    p.add_argument("--input", "-i", required=True,
                   help="Input image (PNG or SVG)")
    p.add_argument("--output", "-o", required=True,
                   help="Output MP4 video")
    p.add_argument("--direction", choices=["horizontal", "vertical"],
                   default="horizontal", help="Panning direction (default: horizontal)")

    # Duration
    dur = p.add_argument_group("duration (choose one)")
    dur.add_argument("--duration", "-d", type=float, default=None,
                     help="Fixed duration in seconds")
    dur.add_argument("--num-reads", type=int, default=None,
                     help="Number of reads in the image (for auto-duration)")
    dur.add_argument("--reads-per-second", type=float, default=18.75,
                     help="Target reads/s scroll rate (default: 18.75)")

    # Zoom
    zm = p.add_argument_group("zoom mode")
    zm.add_argument("--zoom", choices=["fixed", "adaptive"], default="fixed",
                    help="Zoom mode (default: fixed)")
    zm.add_argument("--crop-ratio", type=float, default=0.5,
                    help="Height crop ratio, fixed zoom only (default: 0.5)")
    zm.add_argument("--max-zoom", type=float, default=8.0,
                    help="Max zoom factor, adaptive only (default: 8)")
    zm.add_argument("--top-margin", type=int, default=80,
                    help="Source image top margin px (default: 80, from plot_reads)")
    zm.add_argument("--left-margin", type=int, default=60,
                    help="Source image left margin px (default: 60, from plot_reads)")
    zm.add_argument("--ratio", type=float, default=1 / 300,
                    help="bp-to-pixel ratio in source (default: 1/300)")
    zm.add_argument("--scale-bar-width", type=int, default=100,
                    help="Width of scale bar area in output px (default: 100)")
    zm.add_argument("--zoom-smoothing", type=int, default=200,
                    help="Zoom smoothing window size (0=none, default: 200)")

    # Viewport
    p.add_argument("--viewport-width", type=int, default=None,
                   help="Viewport width px (default: 1920 horiz, 640 vert)")
    p.add_argument("--viewport-height", type=int, default=None,
                   help="Viewport height px (default: 1080 horiz, 864 vert)")
    p.add_argument("--fps", type=int, default=30,
                   help="Frames per second (default: 30)")

    # SVG conversion
    p.add_argument("--svg-height", type=int, default=None,
                   help="Target height for SVG→PNG conversion (default: native)")

    # Overlays
    p.add_argument("--scale-bar", type=str, default=None,
                   help="Static scale bar PNG overlay (fixed zoom only)")
    p.add_argument("--legend", type=str, default=None,
                   help="Legend image (PNG or SVG) below/beside panning area")
    p.add_argument("--background", default="black",
                   help="Background colour (default: black)")
    p.add_argument("--scale-bar-padding", type=int, default=10,
                   help="Padding around scale bar/overlays (default: 10)")

    return p.parse_args()


def main():
    args = parse_args()

    # ── Resolve duration ──────────────────────────────────────────────────
    if args.num_reads:
        duration = args.num_reads / args.reads_per_second
        print(f"Auto-duration: {args.num_reads} reads / {args.reads_per_second} reads·s⁻¹"
              f" = {duration:.1f}s ({duration / 60:.1f} min)")
    elif args.duration:
        duration = args.duration
    else:
        duration = 60
        print("No --duration or --num-reads given; using default 60s")

    # ── Convert SVG inputs to PNG ─────────────────────────────────────────
    temp_files = []

    # Note: Strip-based SVG rendering is disabled - Playwright-based full conversion is used instead
    use_svg_strip_rendering = False

    # Determine viewport width early for legend scaling
    if args.zoom == "adaptive" and args.direction == "horizontal":
        viewport_width = args.viewport_width or 1920
    elif args.direction == "horizontal":
        viewport_width = args.viewport_width or 1920
    else:
        viewport_width = args.viewport_width or 640

    # For strip-based rendering, we handle input differently
    if use_svg_strip_rendering:
        input_path = args.input
        inp_scale = 1.0
        legend_path = args.legend  # Keep as SVG, will be handled in create_adaptive_horizontal_panning_svg
    else:
        # For adaptive zoom with SVG, render at max safe resolution (under 32767 limit)
        svg_height = args.svg_height
        if svg_height is None and args.zoom == "adaptive" and args.input.lower().endswith(".svg"):
            native_w, native_h = get_svg_dimensions(args.input)
            if native_w and native_h:
                # Calculate max scale that stays under 32767 on both dimensions
                max_scale = min(32767 / native_w, 32767 / native_h, 4.0)  # Cap at 4x
                if max_scale > 1.0:
                    svg_height = int(native_h * max_scale)
                    print(f"Scaling SVG to {max_scale:.2f}x for better zoom quality")

        input_path, inp_is_temp, inp_scale = ensure_png(
            args.input, svg_height=svg_height
        )
        if inp_is_temp:
            temp_files.append(input_path)

        legend_path = args.legend
        if legend_path:
            # Render legend SVG at viewport width for crisp display
            legend_path, leg_is_temp, _ = ensure_png(legend_path, svg_width=viewport_width)
            if leg_is_temp:
                temp_files.append(legend_path)

    scale_bar_path = args.scale_bar
    if scale_bar_path:
        scale_bar_path, sb_is_temp, _ = ensure_png(scale_bar_path)
        if sb_is_temp:
            temp_files.append(scale_bar_path)

    # Adjust ratio for SVG→PNG scaling
    effective_ratio = args.ratio * inp_scale

    # ── Dispatch ──────────────────────────────────────────────────────────
    try:
        if use_svg_strip_rendering:
            # Use strip-based SVG rendering for large SVGs
            viewport_width = args.viewport_width or 1920
            viewport_height = args.viewport_height or 1080

            create_adaptive_horizontal_panning_svg(
                svg_path=input_path,
                output_path=args.output,
                duration=duration,
                fps=args.fps,
                viewport_width=viewport_width,
                viewport_height=viewport_height,
                legend_path=legend_path,
                background=args.background,
                top_margin=args.top_margin,
                left_margin=args.left_margin,
                ratio=args.ratio,
                max_zoom=args.max_zoom,
                scale_bar_padding=args.scale_bar_padding,
                scale_bar_width=args.scale_bar_width,
                zoom_smoothing=args.zoom_smoothing,
            )
        elif args.zoom == "adaptive" and args.direction == "horizontal":
            viewport_width = args.viewport_width or 1920
            viewport_height = args.viewport_height or 1080

            create_adaptive_horizontal_panning(
                input_path=input_path,
                output_path=args.output,
                duration=duration,
                fps=args.fps,
                viewport_width=viewport_width,
                viewport_height=viewport_height,
                legend_path=legend_path,
                background=args.background,
                top_margin=int(args.top_margin * inp_scale),
                left_margin=int(args.left_margin * inp_scale),
                ratio=effective_ratio,
                max_zoom=args.max_zoom,
                scale_bar_padding=args.scale_bar_padding,
                scale_bar_width=args.scale_bar_width,
                zoom_smoothing=args.zoom_smoothing,
            )
        elif args.direction == "horizontal":
            viewport_width = args.viewport_width or 1920
            viewport_height = args.viewport_height or 432

            create_horizontal_panning(
                input_path=input_path,
                output_path=args.output,
                duration=duration,
                fps=args.fps,
                crop_ratio=args.crop_ratio,
                viewport_width=viewport_width,
                viewport_height=viewport_height,
                scale_bar_path=scale_bar_path,
                legend_path=legend_path,
                background=args.background,
                scale_bar_padding=args.scale_bar_padding,
            )
        else:  # vertical
            viewport_width = args.viewport_width or 640
            viewport_height = args.viewport_height or 864

            create_vertical_panning(
                input_path=input_path,
                output_path=args.output,
                duration=duration,
                fps=args.fps,
                viewport_width=viewport_width,
                viewport_height=viewport_height,
                scale_bar_path=scale_bar_path,
                legend_path=legend_path,
                background=args.background,
                scale_bar_padding=args.scale_bar_padding,
            )
    finally:
        # Clean up temp PNG conversions
        for f in temp_files:
            if os.path.exists(f):
                os.unlink(f)

    print("\nDone!")


if __name__ == "__main__":
    main()
