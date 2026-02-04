#!/usr/bin/env python3
"""
create_panning_animation.py

Create a seamless looping panning animation from a wide or tall image.
Supports both horizontal (left-to-right) and vertical (top-to-bottom) panning.
Designed for PowerPoint 16:9 widescreen presentations.

Optional features:
- Static scale bar overlay (stays fixed while content pans)
- Static legend image (added below for horizontal, right for vertical)

Usage:
    # Horizontal panning (default)
    python create_panning_animation.py \\
        --input wide_image.png \\
        --output animation.mp4 \\
        --direction horizontal

    # Vertical panning with scale bar and legend
    python create_panning_animation.py \\
        --input tall_image.png \\
        --output animation.mp4 \\
        --direction vertical \\
        --scale-bar scale_bar.png \\
        --legend legend.png
"""

import argparse
import os
import subprocess
import tempfile
from PIL import Image


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create seamless looping panning animation from wide or tall images"
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input PNG image file"
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output MP4 video file"
    )
    parser.add_argument(
        "--direction", choices=["horizontal", "vertical"], default="horizontal",
        help="Panning direction (default: horizontal)"
    )
    parser.add_argument(
        "--duration", "-d", type=int, default=60,
        help="Duration of one pan cycle in seconds (default: 60)"
    )
    parser.add_argument(
        "--fps", type=int, default=30,
        help="Frames per second (default: 30)"
    )
    parser.add_argument(
        "--crop-ratio", type=float, default=0.5,
        help="Ratio of image height to keep from top, horizontal mode only (default: 0.5)"
    )
    parser.add_argument(
        "--viewport-width", type=int, default=None,
        help="Viewport width in pixels (default: 1920 for horizontal, 640 for vertical)"
    )
    parser.add_argument(
        "--viewport-height", type=int, default=None,
        help="Viewport height in pixels (default: 432 for horizontal, 864 for vertical)"
    )
    parser.add_argument(
        "--scale-bar", type=str, default=None,
        help="Static scale bar PNG to overlay at top-left (optional)"
    )
    parser.add_argument(
        "--legend", type=str, default=None,
        help="Legend PNG to add below (horizontal) or right (vertical) of panning area (optional)"
    )
    parser.add_argument(
        "--background", type=str, default="white",
        help="Background color for padding (default: white)"
    )
    parser.add_argument(
        "--scale-bar-padding", type=int, default=10,
        help="Padding around scale bar in pixels (default: 10)"
    )
    return parser.parse_args()


def create_horizontal_panning(
    input_path,
    output_path,
    duration,
    fps,
    crop_ratio,
    viewport_width,
    viewport_height,
    scale_bar_path,
    legend_path,
    background,
    scale_bar_padding,
):
    """Create a seamless horizontal panning animation from a wide image."""

    print(f"Creating horizontal panning animation")
    print(f"  Input: {input_path}")
    print(f"  Output: {output_path}")
    print(f"  Duration: {duration}s, FPS: {fps}")

    # Load image
    img = Image.open(input_path)
    orig_width, orig_height = img.size
    print(f"  Original size: {orig_width} x {orig_height}")

    # Step 1: Crop to top portion
    crop_height = int(orig_height * crop_ratio)
    img_cropped = img.crop((0, 0, orig_width, crop_height))
    print(f"  Cropped to top {crop_ratio*100:.0f}%: {orig_width} x {crop_height}")

    # Step 2: Scale to viewport height while maintaining aspect ratio
    scale_factor = viewport_height / crop_height
    scaled_width = int(orig_width * scale_factor)
    img_scaled = img_cropped.resize((scaled_width, viewport_height), Image.Resampling.LANCZOS)
    print(f"  Scaled to viewport height: {scaled_width} x {viewport_height}")

    # Step 3: Create seamless tile (duplicate horizontally)
    tiled_width = scaled_width * 2
    img_tiled = Image.new('RGB', (tiled_width, viewport_height), background)
    img_tiled.paste(img_scaled, (0, 0))
    img_tiled.paste(img_scaled, (scaled_width, 0))
    print(f"  Tiled for seamless loop: {tiled_width} x {viewport_height}")

    # Save tiled image to temp file
    temp_files = []
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        tiled_path = tmp.name
        temp_files.append(tiled_path)
        img_tiled.save(tiled_path, 'PNG')

    # Calculate pan parameters
    pan_distance = scaled_width
    total_frames = duration * fps
    pixels_per_frame = pan_distance / total_frames

    print(f"  Pan distance: {pan_distance}px over {total_frames} frames")
    print(f"  Speed: {pixels_per_frame:.2f} px/frame")

    # Build ffmpeg filter chain
    filter_parts = []
    inputs = ['-loop', '1', '-i', tiled_path]
    input_idx = 0

    # Base panning crop filter
    crop_filter = (
        f"[{input_idx}:v]crop={viewport_width}:{viewport_height}:"
        f"'mod(t*{pan_distance}/{duration},iw-{viewport_width})':0[panning]"
    )
    filter_parts.append(crop_filter)
    last_stream = "panning"
    input_idx += 1

    # Add scale bar overlay if provided
    if scale_bar_path and os.path.exists(scale_bar_path):
        inputs.extend(['-i', scale_bar_path])
        overlay_filter = (
            f"[{last_stream}][{input_idx}:v]overlay={scale_bar_padding}:{scale_bar_padding}[with_scale]"
        )
        filter_parts.append(overlay_filter)
        last_stream = "with_scale"
        input_idx += 1
        print(f"  Adding scale bar: {scale_bar_path}")

    # Create intermediate panning video or final output
    if legend_path and os.path.exists(legend_path):
        # First create panning video, then vstack with legend
        with tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) as tmp:
            panning_video_path = tmp.name
            temp_files.append(panning_video_path)

        filter_complex = ";".join(filter_parts)
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            *inputs,
            '-filter_complex', filter_complex,
            '-map', f'[{last_stream}]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            panning_video_path
        ]

        print(f"\nGenerating panning video...")
        run_ffmpeg(ffmpeg_cmd)

        # Now vstack with legend
        print(f"  Adding legend: {legend_path}")

        # Load legend and scale to viewport width
        legend_img = Image.open(legend_path)
        legend_width, legend_height = legend_img.size
        if legend_width != viewport_width:
            scale = viewport_width / legend_width
            new_height = int(legend_height * scale)
            legend_img = legend_img.resize((viewport_width, new_height), Image.Resampling.LANCZOS)

        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            legend_scaled_path = tmp.name
            temp_files.append(legend_scaled_path)
            legend_img.save(legend_scaled_path, 'PNG')

        # vstack panning video with legend
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            '-i', panning_video_path,
            '-loop', '1', '-i', legend_scaled_path,
            '-filter_complex', '[0:v][1:v]vstack=inputs=2[out]',
            '-map', '[out]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            output_path
        ]

        print(f"Combining with legend...")
        run_ffmpeg(ffmpeg_cmd)
    else:
        # No legend, output directly
        filter_complex = ";".join(filter_parts)
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            *inputs,
            '-filter_complex', filter_complex,
            '-map', f'[{last_stream}]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            output_path
        ]

        print(f"\nGenerating video...")
        run_ffmpeg(ffmpeg_cmd)

    # Clean up temp files
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.unlink(temp_file)

    # Get output file size
    file_size = os.path.getsize(output_path)
    print(f"  Saved: {output_path}")
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")

    return output_path


def create_vertical_panning(
    input_path,
    output_path,
    duration,
    fps,
    viewport_width,
    viewport_height,
    scale_bar_path,
    legend_path,
    background,
    scale_bar_padding,
):
    """Create a seamless vertical panning animation from a tall image."""

    print(f"Creating vertical panning animation")
    print(f"  Input: {input_path}")
    print(f"  Output: {output_path}")
    print(f"  Duration: {duration}s, FPS: {fps}")

    # Load image
    img = Image.open(input_path)
    orig_width, orig_height = img.size
    print(f"  Original size: {orig_width} x {orig_height}")

    # Step 1: Scale to viewport width while maintaining aspect ratio
    scale_factor = viewport_width / orig_width
    scaled_height = int(orig_height * scale_factor)
    img_scaled = img.resize((viewport_width, scaled_height), Image.Resampling.LANCZOS)
    print(f"  Scaled to viewport width: {viewport_width} x {scaled_height}")

    # Step 2: Create seamless tile (duplicate vertically)
    min_tiled_height = viewport_height + scaled_height
    num_tiles = max(2, (min_tiled_height // scaled_height) + 1)
    tiled_height = scaled_height * num_tiles
    img_tiled = Image.new('RGB', (viewport_width, tiled_height), background)
    for i in range(num_tiles):
        img_tiled.paste(img_scaled, (0, i * scaled_height))
    print(f"  Tiled for seamless loop: {viewport_width} x {tiled_height} ({num_tiles} tiles)")

    # Save tiled image to temp file
    temp_files = []
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        tiled_path = tmp.name
        temp_files.append(tiled_path)
        img_tiled.save(tiled_path, 'PNG')

    # Calculate pan parameters
    pan_distance = scaled_height
    total_frames = duration * fps
    pixels_per_frame = pan_distance / total_frames

    print(f"  Pan distance: {pan_distance}px over {total_frames} frames")
    print(f"  Speed: {pixels_per_frame:.2f} px/frame")

    # Build ffmpeg filter chain
    filter_parts = []
    inputs = ['-loop', '1', '-i', tiled_path]
    input_idx = 0

    # Base panning crop filter (vertical)
    crop_filter = (
        f"[{input_idx}:v]crop={viewport_width}:{viewport_height}:"
        f"0:'mod(t*{pan_distance}/{duration},{pan_distance})'[panning]"
    )
    filter_parts.append(crop_filter)
    last_stream = "panning"
    input_idx += 1

    # Add scale bar overlay if provided
    if scale_bar_path and os.path.exists(scale_bar_path):
        inputs.extend(['-i', scale_bar_path])
        overlay_filter = (
            f"[{last_stream}][{input_idx}:v]overlay={scale_bar_padding}:{scale_bar_padding}[with_scale]"
        )
        filter_parts.append(overlay_filter)
        last_stream = "with_scale"
        input_idx += 1
        print(f"  Adding scale bar: {scale_bar_path}")

    # Create intermediate panning video or final output
    if legend_path and os.path.exists(legend_path):
        # First create panning video, then hstack with legend
        with tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) as tmp:
            panning_video_path = tmp.name
            temp_files.append(panning_video_path)

        filter_complex = ";".join(filter_parts)
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            *inputs,
            '-filter_complex', filter_complex,
            '-map', f'[{last_stream}]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            panning_video_path
        ]

        print(f"\nGenerating panning video...")
        run_ffmpeg(ffmpeg_cmd)

        # Now hstack with legend
        print(f"  Adding legend: {legend_path}")

        # Load legend and scale to viewport height
        legend_img = Image.open(legend_path)
        legend_width, legend_height = legend_img.size
        if legend_height != viewport_height:
            scale = viewport_height / legend_height
            new_width = int(legend_width * scale)
            legend_img = legend_img.resize((new_width, viewport_height), Image.Resampling.LANCZOS)

        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            legend_scaled_path = tmp.name
            temp_files.append(legend_scaled_path)
            legend_img.save(legend_scaled_path, 'PNG')

        # hstack panning video with legend
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            '-i', panning_video_path,
            '-loop', '1', '-i', legend_scaled_path,
            '-filter_complex', '[0:v][1:v]hstack=inputs=2[out]',
            '-map', '[out]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            output_path
        ]

        print(f"Combining with legend...")
        run_ffmpeg(ffmpeg_cmd)
    else:
        # No legend, output directly
        filter_complex = ";".join(filter_parts)
        ffmpeg_cmd = [
            'ffmpeg', '-y',
            *inputs,
            '-filter_complex', filter_complex,
            '-map', f'[{last_stream}]',
            '-t', str(duration),
            '-r', str(fps),
            '-c:v', 'libx264',
            '-preset', 'medium',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            output_path
        ]

        print(f"\nGenerating video...")
        run_ffmpeg(ffmpeg_cmd)

    # Clean up temp files
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.unlink(temp_file)

    # Get output file size
    file_size = os.path.getsize(output_path)
    print(f"  Saved: {output_path}")
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")

    return output_path


def run_ffmpeg(cmd):
    """Run ffmpeg command with error handling."""
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"  ffmpeg error: {e.stderr}")
        raise


def main():
    args = parse_args()

    # Set direction-specific defaults
    if args.direction == "horizontal":
        viewport_width = args.viewport_width or 1920
        viewport_height = args.viewport_height or 432

        create_horizontal_panning(
            input_path=args.input,
            output_path=args.output,
            duration=args.duration,
            fps=args.fps,
            crop_ratio=args.crop_ratio,
            viewport_width=viewport_width,
            viewport_height=viewport_height,
            scale_bar_path=args.scale_bar,
            legend_path=args.legend,
            background=args.background,
            scale_bar_padding=args.scale_bar_padding,
        )
    else:  # vertical
        viewport_width = args.viewport_width or 640
        viewport_height = args.viewport_height or 864

        create_vertical_panning(
            input_path=args.input,
            output_path=args.output,
            duration=args.duration,
            fps=args.fps,
            viewport_width=viewport_width,
            viewport_height=viewport_height,
            scale_bar_path=args.scale_bar,
            legend_path=args.legend,
            background=args.background,
            scale_bar_padding=args.scale_bar_padding,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
