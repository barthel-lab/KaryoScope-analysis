# Audit: create_panning_animation.py (1195 lines)

## 1. Purpose
Turns a tall/wide PNG (a rendered read plot) into a panning MP4 for 16:9
PowerPoint presentations. Supports a seamless looping **fixed-zoom** pan
(horizontal or vertical) and an **adaptive-zoom** horizontal pan that
dynamically zooms so the tallest currently-visible read fills the viewport
(docstring lines 2-35). Optionally composites a legend and draws a dynamic
vertical scale bar that re-labels itself as zoom changes.

## 2. CLI surface
`argparse` in `parse_args` (lines 991-1062), `RawTextHelpFormatter`. Single
script; `main` (1065-1191) dispatches by `(--zoom, --direction)`.

- I/O: `--input/-i` (PNG only, 998), `--output/-o` (MP4, 1000),
  `--direction` {horizontal,vertical} (1002).
- Duration group (1006-1012): `--duration/-d` OR `--num-reads` +
  `--reads-per-second` (default 20.625). Resolution priority: num-reads →
  duration → default 60 s (lines 1069-1077).
- Zoom group (1015-1042): `--zoom` {fixed,adaptive} (1016), `--crop-ratio`
  (0.5, fixed only, 1018), `--max-zoom` (auto if unset, 1020),
  `--scale-bar-max-bp` (20000, 1022), `--scale-bar-fraction` (0.1, 1024),
  `--top-bias` (0.0, 1026), `--top-margin` (80, 1029), `--left-margin` (60,
  1031), `--ratio` (1/300, 1033), `--scale-bar-width` (100, 1035),
  `--zoom-smoothing` (200, 1037), `--uniform-zoom` (default True, 1039),
  `--vertical-only-zoom` (legacy, 1041).
- Viewport: `--viewport-width`/`--viewport-height` (auto by mode, 1045-1048),
  `--fps` (30, 1049).
- Overlays: `--scale-bar` (static PNG, fixed only, 1053), `--legend`
  (PNG/SVG, 1055), `--background` (black, 1057), `--scale-bar-padding` (10,
  1059).

## 3. Inputs & outputs
Inputs:
- Main input MUST be PNG; non-PNG raises `ValueError` (lines 1080-1084).
- Legend may be PNG or SVG; SVG is converted via `ensure_png`/`svg_to_png`
  (lines 1098-1101, 132-148, 72-129).
- Static `--scale-bar` PNG/SVG (fixed-zoom path only, 1103-1107).

Outputs:
- One MP4 (`output_path`), H.264 / yuv420p, CRF 18, preset medium (encoder
  args lines 478-490 adaptive; 766-817 / 794-805 fixed). No audio.

## 4. Pipeline / control flow (key functions + line numbers)
- `main` (1065): resolve duration (1069-1077) → validate PNG (1080-1084) →
  pre-compute viewport width for legend scaling (1088-1093) → convert legend /
  scale-bar SVGs (1096-1107) → dispatch (1110-1184) → `finally` cleanup of temp
  PNGs (1185-1189).
- Dispatch branches: adaptive+horizontal →
  `create_adaptive_horizontal_panning` (computes auto `max_zoom`, 1118-1131);
  fixed+horizontal → `create_horizontal_panning`; vertical →
  `create_vertical_panning`.
- `create_adaptive_horizontal_panning` (355-669): the core renderer.
  - Load PNG → numpy (398-401).
  - `detect_content_heights` (224-256): vectorised per-column content height
    below `top_margin`, zeroing `left_margin`.
  - `build_height_profile` (259-292): forward-fill gaps, enforce monotonic
    non-increasing envelope (reads sorted longest→shortest).
  - `_smooth_profile` EMA (305-312) when `--zoom-smoothing>0`.
  - `_find_content_bounds` (295-300) sets pan range.
  - Legend load/resize/even-pad (434-449); even-dim enforcement for H.264
    (451-458).
  - Spawns **persistent ffmpeg** via `Popen` reading rawvideo from stdin
    (478-491); **background thread drains stderr** to avoid pipe deadlock
    (492-499).
  - Per-frame loop (511-655): compute content height at pan x, vertical zoom
    (519-531, with `top_bias` boost), crop (534-558, header band via
    `include_header`), resize with cv2 or PIL (568-579), assemble frame
    (581-584), draw dynamic scale bar with rotated-label cache (586-637),
    stack legend (639-648), `process.stdin.write(out.tobytes())` (650).
  - Close stdin, `wait()`, join stderr, check returncode (657-664).
- `create_horizontal_panning` (674-826): crop top `crop_ratio`, scale to
  viewport height, **tile x2 for seamless loop** (704-709), drive ffmpeg
  `crop=...:'mod(t*dist/dur,iw-w)':0` filter (729-732); optional scale bar via
  `draw_dynamic_scale_bar` overlay, optional legend via a **second ffmpeg pass**
  + `vstack` (761-805).
- `create_vertical_panning` (829-975): scale to viewport width, tile vertically
  (853-859), `crop=...:0:'mod(...)'` (878-881), optional legend via `hstack`
  (910-954).
- `run_ffmpeg` (980-986): one-shot `subprocess.run(check=True)` with stderr
  tail on error.
- SVG helpers `get_svg_dimensions` (61-69), `svg_to_png` (72-129),
  `ensure_png` (132-148).
- Font/scale-bar helpers `_load_font` (153-162), `_nice_scale_value` (165-171),
  `draw_dynamic_scale_bar` (174-219).
- Strip helper `compute_strip_zoom_levels` (317-350) — **defined but never
  called** (dead code).

## 5. Key design decisions (cite lines + WHY if stated)
- **Streaming rawvideo to a persistent ffmpeg** in the adaptive path (478-491)
  instead of writing PNG frames to disk: avoids materialising thousands of
  frames; lets per-frame zoom vary. WHY (implied): adaptive zoom can't be done
  with a single ffmpeg crop expression.
- **Stderr drained in a daemon thread** (492-499, comment line 492): explicitly
  to "avoid pipe buffer deadlock" — ffmpeg's stderr would otherwise fill and
  block when we never read it during a long stdin write.
- **cv2 used when available** (50-54, 568-573) "cv2 is ~5-10x faster" (comment
  568) than PIL resize; PIL bilinear fallback (574-577).
- **`Image.MAX_IMAGE_PIXELS = None`** (line 57) "Allow very large images" —
  disables Pillow decompression-bomb guard for the multi-GB read PNGs.
- **Even dimensions enforced** for H.264/yuv420p (legend pad 443-447; viewport
  451-458; legend even-w/h 784, 933). WHY: libx264 yuv420p requires even W/H.
- **Rotated scale-bar label cache** (509, 623-632): rotation is expensive;
  cache keyed by label string — big win since label changes rarely.
- **Auto max_zoom from scale-bar geometry** (1118-1131): derives zoom so
  `scale_bar_max_bp` fills `scale_bar_fraction` of the viewport at max zoom.
- **Monotonic non-increasing height profile** (282-287): exploits that reads are
  pre-sorted longest→shortest, giving a smooth zoom-out as the pan advances.
- **Seamless loop via tiling** in fixed mode (704-709, 853-859): duplicate the
  scaled image so the `mod()` crop wraps without a visible seam.

## 6. Assumptions (frame/encoding assumptions)
- Input PNG layout matches `plot_reads`/telogator: `top_margin` rows are
  header, `left_margin` cols are axis labels (defaults 80/60, 1029-1032), reads
  start below/right of those. `detect_content_heights` zeroes those regions
  (224-256). Wrong margins → wrong zoom / clipped content.
- Background is pure black `(0,0,0)` or pure white `(255,255,255)`; content
  detection is `any(pixel != 0)` / `!= 255` (241-245). Anti-aliased near-bg
  pixels count as content; a non-pure background breaks detection.
- Reads are sorted longest→shortest (monotonic profile, 282-287); an unsorted
  image yields a clamped/incorrect envelope.
- `--ratio` (1/300) must match the source PNG's actual bp→px scale, including
  any `png_scale` applied upstream; `plot_reads._run_animation` multiplies
  `ratio * height_scale` before calling (see §7). Manual coupling — no
  validation that the value is correct.
- ffmpeg, and `rsvg-convert` (legend SVG only), are on PATH. cv2 optional.
- H.264 even-dimension requirement; odd legend/viewport sizes are padded.
- numpy array channel count assumed 3 (RGB): `detect_content_heights` unpacks
  `h, w, _` (236) and frames are RGB24; an RGBA or grayscale PNG would break
  (the adaptive path does **not** `.convert("RGB")` the main input, unlike the
  legend at line 435).

## 7. Dependencies
External libs: `numpy` (46), `Pillow` (`Image, ImageDraw, ImageFont`, 47),
**optional `cv2`/opencv** (50-54, `HAS_CV2`), `argparse`, `math` (imported,
line 38 — appears **unused**), `os`, `re`, `subprocess`, `sys`, `tempfile`,
`pathlib`, plus inline `threading, io` (493 — `io` unused).

External tools — **ffmpeg** (critical):
- Adaptive: persistent `Popen(["ffmpeg","-y","-f","rawvideo","-pix_fmt",
  "rgb24","-s", WxH,"-r",fps,"-i","-","-c:v","libx264","-preset","medium",
  "-crf","18","-pix_fmt","yuv420p", out])` reading raw RGB frames from stdin
  (478-491); stderr drained in a thread (492-499).
- Fixed horizontal: `crop` filter with `mod()` pan expr (729-732), optional
  scale-bar `overlay` (752-755), optional legend as a **2nd ffmpeg pass** +
  `vstack` (766-805). `-t int(duration)`, `-loop 1` on the tiled PNG input.
- Fixed vertical: same shape, `crop` x=0 + y-`mod()` (878-881), legend via
  `hstack` (947). Encoder args duplicated in 4 places (libx264/medium/crf18/
  yuv420p).
- `rsvg-convert` for legend/scale-bar SVG→PNG with
  `RSVG_MAX_LOADED_ELEMENTS=100000000` env (99-109).

karyoplot usage: **none currently.** This script predates the migration and
keeps its own `_load_font` and SVG/scale-bar helpers.

Inter-script / shared deps: **`KaryoScope_plot_reads.py` imports this module
directly** via a `sys.path.insert(0, script_dir)` hack
(`_run_animation`, plot_reads lines 1589-1594), calling
`create_horizontal_panning`, `create_vertical_panning`,
`create_adaptive_horizontal_panning`. plot_reads passes a `scaled_ratio =
args.ratio * height_scale` and scaled top/left margins (plot_reads 1620-1644)
— so this module is effectively a library, not just a script. Its PNG input is
produced by `plot_reads --format png` or `KaryoScope_telogator_reads_viz.py`.

## 8. Proposed home in new layout
- Subcommand: `karyoscope-analysis animate` (or `panning-animation`). Thin
  `commands/animate.py` (click) → `core/animation/` package.
- Decomposition of `core/animation/`:
  - `profile.py`: `detect_content_heights`, `build_height_profile`,
    `_find_content_bounds`, `_smooth_profile` (pure numpy — easily testable).
  - `adaptive.py`: `create_adaptive_horizontal_panning` (frame loop + ffmpeg
    streaming).
  - `fixed.py`: `create_horizontal_panning`, `create_vertical_panning`.
  - `encoder.py`: an ffmpeg wrapper (consolidate the 4 duplicated encoder arg
    lists + the stderr-draining Popen pattern + `run_ffmpeg`).
  - Replace the `plot_reads → sys.path.insert` import hack with a normal
    package import once both live in `karyoscope_analysis`.
- karyoplot push-down candidates + module:
  - `_load_font` (153-162) → **delete**, use `karyoplot.core.fonts.pil_font`
    (it is a duplicate of that function's logic; telogator already migrated).
  - `svg_to_png` / `ensure_png` / `get_svg_dimensions` (61-148) →
    `karyoplot.svg.export` already has `svg_to_png` +
    `is_rsvg_convert_available`; extend it with the `width/height/scale-factor`
    + native-fallback behaviour and delete the local copies.
  - `_nice_scale_value` + `draw_dynamic_scale_bar` (165-219) and the inline
    per-frame scale-bar block (586-637) → a shared raster scale-bar helper in
    `karyoplot.svg`/a new `karyoplot.raster`, unified with the telogator
    script's two scale-bar functions (four implementations total today).
  - ffmpeg encoder wrapper is a strong **new shared module** candidate (no home
    in karyoplot yet — it has zero video code); could live in
    `karyoplot` (e.g. `karyoplot.video`) or stay analysis-local if considered
    app-specific. Flag for user decision.

## 9. Smells / risks / dead code / duplication (line-cited)
- **Dead code**: `compute_strip_zoom_levels` (317-350) is never called. The
  "Strip-based pre-rendering" section header (315) is vestigial.
- **Unused imports**: `math` (38), `io` (493). `ImageFont`/`ImageDraw` are used.
- **PIL pixel limit deliberately disabled** (line 57) — necessary here, but note
  the matching telogator PNG writer does **not** set it (cross-script
  inconsistency).
- **subprocess deadlock**: handled in the adaptive path (stderr thread,
  492-499) but the **fixed** paths use `run_ffmpeg`/`subprocess.run(
  capture_output=True)` (983) which buffers all output in memory — fine for the
  short fixed-mode commands but inconsistent.
- **No ffmpeg/rsvg presence check** before invocation (unlike
  `karyoplot.svg.export.is_rsvg_convert_available`); missing ffmpeg surfaces as
  a raw `FileNotFoundError` from `Popen`/`run`, and on the adaptive path the
  stderr thread + stdin write can raise `BrokenPipe` confusingly.
- **Massive duplication** of encoder arg lists (libx264/preset/crf/pix_fmt) in
  ~4 spots (485-489, 772-775, 800-803, 920-923, 950-953) and of the
  legend-resize-to-even logic (780-792 vs 927-941).
- **`create_adaptive_horizontal_panning` is ~315 lines** (355-669) with a deep
  per-frame loop — hard to read/test; mixes layout math, cropping, drawing, and
  encoding.
- **No input `.convert("RGB")`** on the main image (398-399); an RGBA/grayscale
  PNG would break `detect_content_heights`/frame assembly (legend *is* converted
  at 435).
- **`--vertical-only-zoom` vs `--uniform-zoom`**: `--uniform-zoom` defaults True
  and `store_true` can't turn it off; only `--vertical-only-zoom` flips it
  (1116). The `--uniform-zoom` flag is effectively a no-op/confusing.
- **`--scale-bar-padding`, `--max-zoom`, `--top-bias`, `--scale-bar-*`** apply
  only to specific zoom modes; no validation warns when a flag is irrelevant to
  the chosen mode (e.g. `--crop-ratio` ignored in adaptive).
- **Temp-file leakage risk**: fixed-mode temp files cleaned in a plain loop at
  function end (819-821, 968-970); if ffmpeg raises mid-way they leak (no
  `try/finally` there, unlike `main`'s outer `finally` for legend/scale-bar at
  1185-1189).
- **`-t int(duration)`** truncates fractional seconds (771, 800, 920, 949),
  while the adaptive path uses `int(duration*fps)` frames — slight duration
  mismatch between modes.

## 10. Testability notes
- Pure functions are highly testable: `detect_content_heights`,
  `build_height_profile`, `_find_content_bounds`, `_smooth_profile`,
  `_nice_scale_value`, `compute_strip_zoom_levels` (numpy in/out) — unit test
  with small synthetic arrays.
- `get_svg_dimensions` testable with a tiny SVG fixture.
- The encoder paths need ffmpeg; in CI either skip (mark `requires_ffmpeg`) or
  test the *command construction* by extracting an encoder-arg builder
  (currently inlined) and asserting the arg list. A short 1-2 s smoke encode is
  feasible as an integration test if ffmpeg is present.
- The persistent-Popen + stderr-thread logic is the trickiest part; extracting
  it into an `encoder` object with an injectable runner would make it unit
  testable without real ffmpeg.
- Font/rsvg dependencies must degrade or skip on CI (Barthel fonts absent;
  `pil_font` already handles this once migrated).
- Output validation: assert MP4 exists, non-zero size, and (with ffprobe)
  duration/dimensions — good golden-ish integration checks.

## 11. Open questions for the user
- Should the ffmpeg encoder wrapper live in `karyoplot` (new `video` module) or
  stay analysis-local? karyoplot currently has no video code.
- Is the **fixed-zoom** mode still used, or has **adaptive** superseded it? If
  legacy, we could drop `create_horizontal/vertical_panning`,
  `--crop-ratio`, `--vertical-only-zoom`, and the static `--scale-bar` overlay,
  removing most of the duplication.
- Can `compute_strip_zoom_levels` (dead, 317-350) be deleted?
- Should the `plot_reads → create_panning_animation` `sys.path.insert` import
  (plot_reads 1589-1594) be replaced by a proper package import as part of this
  reorg (they'd share `karyoscope_analysis.core.animation`)?
- Are the manually-coupled defaults (`--ratio` 1/300, `--top-margin` 80,
  `--left-margin` 60) safe to centralise so telogator/plot_reads/animation can't
  drift?
- Confirm only black/white pure backgrounds are supported (content detection
  assumption), or do we need a tolerance-based detector?
