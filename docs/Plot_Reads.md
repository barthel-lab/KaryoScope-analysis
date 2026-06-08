# Plot Reads & Panning Animations

This page documents the two scripts that turn per-read feature BED files into static read-stack figures and panning MP4 animations:

- `scripts/KaryoScope_plot_reads.py` — the top-level CLI. Loads BED files, lays out reads, draws the figure (SVG/PNG), composes the legend, and optionally drives the animation.
- `scripts/create_panning_animation.py` — the panning/zoom engine. Given a wide PNG, produces an MP4 by cropping a moving viewport across the image.

When `--animate` is set on the plot script, the two are wired together end-to-end: the plot script renders the PNG and immediately calls into `create_panning_animation.create_adaptive_horizontal_panning(...)` with all calibration parameters (`ratio`, `top_margin`, `left_margin`, `png_scale`) passed through. **Running the animation as a separate manual step is supported but error-prone** — see [Calibration pitfalls](#calibration-pitfalls).

---

## Quick start

```bash
python scripts/KaryoScope_plot_reads.py \
  --bed SAMPLE1:path/to/sample1.bed.gz SAMPLE2:path/to/sample2.bed.gz \
  --colors resources/KS_human_CHM13_v2/KS_human_CHM13_v2.region_subtelomere_flat.colors.txt \
  --read-list reads_of_interest.readlist.tsv \
  --output figure.png \
  --format png \
  --legend \
  --orient-telomere-top \
  --background black
```

Add `--animate --animate-viewport 1920x1080` to also produce `figure.mp4` next to the PNG.

---

## Inputs

### `--bed` (required)

One or more BED files in the form `SAMPLE:PATH` (sample name explicit) or just `PATH` (sample name derived from filename). Files may be gzipped. Each BED must have the standard 4 columns: `chrom, start, end, feature` — `feature` strings are looked up in the colors file.

Alternative input modes:

- `--samples NAME1 NAME2 … --results-dir DIR` — auto-discovers BED files under a results directory.
- `--cluster-prefix PREFIX --clusters 22,69,112 --results-dir DIR` — picks BEDs by cluster ID from a clustering run.

### `--colors` (required)

Path to a feature-to-color mapping file. See [Colors file format](#colors-file-format) below.

### `--read-list` (optional)

Plain TSV (gzip supported) selecting which reads from the BEDs to draw and (optionally) how to group them.

- **1 column** — read IDs only. All reads are drawn in one block per sample.
- **2 columns** — `read_id <tab> sample`. Same as above; the sample column is informational.
- **3+ columns** — `read_id <tab> sample <tab> group [<tab> subgroup …]`. Reads are grouped two-level for visual separation (e.g. satellite type / cell line).
- Header rows that begin with `sequence`, `read`, or `read_id` are skipped automatically.
- Extra trailing columns become heatmap tracks when `--heatmap` is set.

---

## CLI reference

Arguments are grouped here by purpose. All argparse definitions live in `scripts/KaryoScope_plot_reads.py:1655-2010`.

### Layout & sizing

| Flag | Default | Effect |
|---|---|---|
| `--bar-width N` | — | Horizontal width of each read bar in pixels (source units). |
| `--read-spacing N` | — | Pixel gap between adjacent reads. |
| `--sample-spacing N` | — | Pixel gap between sample blocks. |
| `--subgroup-spacing N` | — | Pixel gap between subgroup blocks (when `--read-list` provides a 3rd/4th column). |
| `--ratio FLOAT` | `1/300` | bp-to-pixel ratio for read **height**. Drives both the read drawing and the scale bar. Changing this without regenerating the source PNG breaks scale-bar calibration in the animation. |
| `--viewport-ratio W:H` | `16:9` | Aspect ratio for the output figure. `none` disables and uses raw content dimensions. |
| `--png-scale FLOAT` | `8.0` | PNG resolution multiplier. Used both when rasterizing the SVG and as the `height_scale` factor that the animation passes back into the panning math. |
| `--font-size N` | `11` | Label / scale-bar text size. |
| `--horizontal` | off | Draw reads horizontally (1 read per row) instead of the default vertical layout (1 read per column). |
| `--no-header` | off | Skip sample labels and separator lines above the reads. |
| `--no-scale-bar` | off | Skip the left-margin scale bar in the source PNG. |
| `--scale-bar-bp N` | auto | Override the scale-bar physical size in bp (default picks 10 kb / 5 kb based on read lengths). |
| `--background COLOR` | `black` | `white`, `black`, or `both`. `both` renders one PNG per background. |

### Orientation (one of)

| Flag | Behavior |
|---|---|
| `--orient-telomere-top` | Flip reads so canonical/noncanonical telomere features end up at the top. |
| `--orient-chromosome-top` | Flip reads so chromosome-specific (arm) features end up at the top. |
| `--orient-satellite-top` | Flip reads so the satellite-dense end is at the top (telomere first, satellite density as fallback). |

### Output

| Flag | Default | Effect |
|---|---|---|
| `--output PATH` / `-o PATH` | required | Output path. Extension is auto-adjusted to `--format`. |
| `--format {svg,png,both}` | `svg` | `--animate` implies `png`. |
| `--scale-bar-output PATH` | — | Write a separate scale-bar SVG/PNG. |

### Legend

| Flag | Effect |
|---|---|
| `--legend` | Enables an auto-filtered, section-grouped color legend. Vertical reads → legend below in a multi-column grid; horizontal reads → legend to the right in a single column. See [Legend behavior](#legend-behavior) below. |

The legend is **not** controlled by additional flags — column count, subheading text, column ordering, and label display strings are all derived from the colors file plus the rendered reads. See [Legend behavior](#legend-behavior).

### Animation

| Flag | Default | Effect |
|---|---|---|
| `--animate` | off | Renders an MP4 alongside the PNG. |
| `--animate-direction {horizontal,vertical}` | auto | Pan direction. Auto: horizontal pan for vertical reads, vertical pan for horizontal reads. |
| `--animate-duration SEC` | auto | Override total duration. |
| `--animate-reads-per-second FLOAT` | `20.625` | Used to auto-compute duration when `--animate-duration` is not set. Halve to slow the pan to 50 %, quarter for 25 %. |
| `--animate-viewport WxH` | auto | Viewport in pixels. PPT widescreen is `1920x1080`; a "16:4" wide-strip variant is `1920x480`. |
| `--animate-fps N` | `30` | Output FPS. |
| `--animate-zoom {fixed,adaptive}` | `fixed` | `adaptive` adjusts vertical zoom per-frame to fit local read heights into the viewport. `fixed` uses a single global zoom based on the tallest read in the dataset, so the camera doesn't visibly "breathe" as it pans. |
| `--animate-crop-ratio FLOAT` | `0.5` | Height crop ratio for fixed zoom. |
| `--animate-legend PATH` | — | Composite an external legend PNG/SVG into the animation. Ignored when `--legend` is set (the built-in legend is composited automatically). |
| `--animate-include-header` / `--no-animate-include-header` | on | Include the header band (sample labels, separator lines) above the reads in each frame. Off restores the older behavior that crops to the reads top. |

### Feature rendering

| Flag | Default | Effect |
|---|---|---|
| `--feature-mode {raw,transition,smooth}` | `smooth` | `smooth` = windowed majority-vote downsampling (recommended); `transition` = direct scaling with min-width enforcement; `raw` = integer pixel stacking. |
| `--min-feature-width FLOAT` | `0.5` | Minimum pixel width per feature in `transition` mode. |
| `--min-width-exclude PAT [PAT …]` | `novel *arm* ct*` | Glob patterns for features exempt from min-width inflation. |
| `--oversample N` | `1` | Oversampling factor for `smooth` mode. |
| `--read-border` | off | Draw a thin black border around each read bar. |

### Filtering

| Flag | Effect |
|---|---|
| `--filter-group VAL` | Restrict to reads whose group column matches (repeatable). |
| `--min-length BP` | Drop reads shorter than this. |
| `--max-length BP` | Drop reads longer than this. |

---

## Colors file format

A plain text file with two whitespace-separated columns plus optional comment-line section markers. The first non-comment row may optionally be a `feature color` header — it is detected case-insensitively and skipped.

```text
feature	color
# Subtelomere
canonical_telomere	#008000
noncanonical_telomere	#9370DB
ITS	#FFFF00
TAR1	#FF00FF
telomere_like	#FFA500
nonsubtelomeric	#808080
# Region (satellite)
centromeric	#B0C4DE
aSat	#E60000
alpha_hor	#FF2020
…
# Region (structural)
ct	#B8AFA2
arm	#808080
…
# Other
novel	#ffffff
```

Parsing happens in `karyoplot.core.colors.load_palette_file()` at `KaryoScope-plotlib/src/karyoplot/core/colors.py:124-166`. The function recognizes three line shapes:

1. `# <text>` — section marker. The text after `# ` becomes the section header for all subsequent features until the next marker.
2. `feature color` — feature-to-color mapping.
3. Blank lines and the `feature color` header row are ignored.

### Suffix conventions

Two suffix conventions are baked into the parser and the rest of the pipeline:

- **`_specific`** — a feature like `canonical_telomere_specific` is treated as an alias for `canonical_telomere`. The parser writes both keys into the palette pointing at the same color (line 152-153).
- **`_multigroup1`, `_multigroup2`, …** — used when the same feature appears in multiple section groups. Only the trailing digits change. The legend display strips this suffix (see below).

When `suffix_both_ways=True` is set (the default in `KaryoScope_plot_reads.py`), the loader also creates `feature_specific` → same color for every bare `feature` that doesn't already carry one of these suffixes.

---

## Legend behavior

This is the most commonly misunderstood part of the pipeline. The legend is **not** a separate config — every aspect of it (which entries appear, what they're labeled, what subheadings they sit under, how columns wrap) is derived from the colors file plus the rendered reads.

### Subheading derivation

Subheading text **is the text of each `# ...` comment line in the colors file, verbatim**.

- Parser: `karyoplot/core/colors.py:127-135`. The regex is `^#\s+(.+)` — anything after `# ` (one or more spaces) becomes the header string for the section that follows.
- The section ends at the next `#` comment line; everything between two comments belongs to the first comment's header.
- If a colors file has no `#` comments at all, the loader wraps everything in a single section with `header=None` (line 162-164), and the legend renders with no subheadings.
- A feature line with no preceding comment lands in a leading `header=None` section.

In the v2 region-subtelomere-flat file shown above this produces the subheadings:

| Subheading | Features in section |
|---|---|
| Subtelomere | `canonical_telomere`, `noncanonical_telomere`, `ITS`, `TAR1`, `telomere_like`, `nonsubtelomeric` |
| Region (satellite) | `centromeric`, `aSat`, `alpha_hor`, `active_hor`, `dhor`, `hor`, `mixedAlpha`, `mon`, `bSat`, `cenSat`, `gSat`, `HSat`, `HSat1`, `HSat1A`, `HSat1B`, `HSat2`, `HSat3`, `rDNA` |
| Region (structural) | `ct`, `arm`, `p_arm`, `q_arm` |
| Other | `novel` |

**Implication:** to rename a subheading, change the comment in the colors file. To add one, insert a `# New header` line above the features that should sit under it. To remove all subheadings, delete every `#` comment.

### Feature visibility (auto-filter)

The legend only shows features that **actually appear in the rendered reads** — not every feature in the colors file. Implemented in `_filter_legend_features()` at `KaryoScope_plot_reads.py:1165-1207`:

```python
used_bases = set()
for feat in features_used:
    base = re.sub(r'_specific$', '', feat)
    used_bases.add(base)

…
for header, section_feats in color_sections:
    for feat in section_feats:
        base = re.sub(r'_specific$', '', feat)
        if base not in used_bases:
            continue
        display = re.sub(r'_multigroup\d+$', '', base).replace('_', ' ')
        …
```

The order in which sections and entries appear in the legend is the **order they appear in the colors file**. Within a section, features that aren't present in the BED files are skipped silently.

### Display label munging

Three transformations are applied in order to produce the human-readable label (`KaryoScope_plot_reads.py:1186-1189`):

1. Strip the trailing `_specific` if present.
2. Strip the trailing `_multigroup1`, `_multigroup2`, … if present.
3. Replace remaining underscores with spaces.

Examples:

| Colors-file feature | Legend label |
|---|---|
| `canonical_telomere` | `canonical telomere` |
| `array_multigroup1` | `array` |
| `telomere_like_multigroup1` | `telomere like` |
| `noncanonical_telomere_specific` | `noncanonical telomere` |
| `HSat1A` | `HSat1A` (no underscore, no suffix) |

After the label is computed, a `seen_display` set deduplicates collisions — so if both `canonical_telomere` and `canonical_telomere_specific` appear in the BEDs, only one entry is emitted.

### Column layout

`compute_legend_layout()` (`KaryoScope_plot_reads.py:1210-1328`) lays out the filtered items as a grid. The important rules:

1. **Items fill top-to-bottom, then left-to-right** within a column.
2. **Each section header starts a new column** (lines 1248-1256 build per-section blocks). Sections do not share columns — this is why a 4-section colors file always produces at least 4 columns in the legend.
3. The maximum column count defaults to **6** (`max_cols=6` parameter at line 1210). If the available width permits, the layout uses up to that many columns; otherwise it uses fewer.
4. Column width is uniform and sized to the longest label across all sections.
5. The number of rows is initially `ceil(total_slots / num_cols)`. If a single section would need more columns than allowed, `num_rows` grows until the layout fits; if it still doesn't fit, `num_cols` expands beyond `max_cols` rather than overlay (lines 1281-1284).

In horizontal-read mode, the legend is forced to a **single column** placed to the right of the reads. In vertical-read mode it's the multi-column grid placed below.

### Colors-file duplicates

The parser does **not** detect duplicate hex codes across features. If two features share the same color in the file (e.g. `telomere_like_multigroup1 #FFA500` and `DJ #FFA500`), both entries appear in the legend with the same swatch — visually ambiguous but not flagged. Audit colors files manually before committing.

---

## Source PNG anatomy

When `--format png` is used, the rendered PNG has the following structure (sizes given relative to the `--png-scale` multiplier and `--ratio`):

```
┌─────────────────────────────────────────────────────────┐
│                  Top margin (header band)               │  top_margin
│   ┌────┐ Sample labels, separator lines                 │
├───┤ SB ├─────────────────────────────────────────────────┤
│   │    │                                                 │
│ L │    │           Reads (vertical bars)                 │  content
│ M │    │                                                 │  height
│   └────┘                                                 │
├─────────────────────────────────────────────────────────┤
│                       Legend                             │
└─────────────────────────────────────────────────────────┘
```

- **`top_margin`** — header band with sample labels and horizontal separator lines. Controlled by `--no-header` (skip entirely) and `--png-scale` (height scaling). The animation engine needs to know this height to position the dynamic scale bar correctly.
- **`left_margin` (LM)** — left padding holding the embedded scale bar ("SB"). Drawn by `draw_scale_bar()` inside the rasterizer (`KaryoScope_plot_reads.py:2670` and surrounding). Controlled by `--no-scale-bar` (skip) and `--scale-bar-bp` (override physical size).
- **Reads** — vertical (default) or horizontal (with `--horizontal`) bars colored by per-position feature lookup.
- **Legend** — composited below (vertical mode) or to the right (horizontal mode) when `--legend` is set.

The embedded scale bar lives entirely inside `left_margin`. The animation engine ignores `[0, left_margin]` x-coordinates when computing the panning range — see [Calibration pitfalls](#calibration-pitfalls).

---

## Animation pipeline

When `--animate` is set, `KaryoScope_plot_reads.py:1633-1644` calls into `create_adaptive_horizontal_panning()`:

```python
create_adaptive_horizontal_panning(
    png_path, mp4_path, duration, args.animate_fps, vw, vh,
    legend_path=legend, background=bg,
    top_margin=scaled_top_margin,
    left_margin=scaled_left_margin,
    ratio=scaled_ratio,
    max_zoom=1.0, scale_bar_padding=10,
    include_header=args.animate_include_header,
    zoom_mode=args.animate_zoom)
```

Each parameter is **derived from the same arguments and config that generated the source PNG** — that's how the math stays consistent.

### Viewport composition

The animation viewport is split into:

```
[ scale_bar_width (100 px) | panning_area ]
                           = viewport_width
```

The scale bar is drawn live by the animation code on the left 100 px; the panning area shows the cropped/resized source PNG.

Total frame height is `viewport_height + legend_height` (the legend, if any, is composited below the panning area as a static strip on every frame).

### Pan range

```python
content_start, content_end = _find_content_bounds(profile)
pan_start = max(left_margin, content_start - 10)
pan_end   = content_end
```

The `max(left_margin, …)` clamp ensures the camera never lingers on the source PNG's embedded scale bar — without it, frame 0 shows two scale bars stacked next to each other (the dynamic one in the viewport's left strip plus the embedded one bleeding in from the source).

### Zoom modes (`--animate-zoom`)

`adaptive` (legacy): per-frame zoom

```python
ch = max(min_content_h, profile[min(x_preview, len(profile) - 1)])
v_zoom = min(available_height / (ch + content_padding), max_zoom)
h_zoom = v_zoom        # uniform zoom
```

`ch` (content height at the current x position) is sampled per frame, so the zoom adapts to local read heights. At a 1080-px viewport this is usually clamped to `max_zoom=1.0` because the reads fit; at a 480-px viewport the per-frame variation becomes visible as "zoom in/out" while panning.

`fixed`: global zoom

```python
global_max_ch = max(min_content_h, float(profile.max()))
ch = global_max_ch
```

A single zoom level is computed from the tallest read in the dataset and used for every frame. The camera pans without "breathing." This is the default in this codebase.

### Scale bar (live, in the viewport)

The dynamic scale bar is drawn each frame at `bar_y = reads_top_y + header_vp_h`, where `header_vp_h` is the viewport projection of `top_margin`:

```python
header_vp_h = int(top_margin * available_height / max(1, top_margin + crop_h)) \
              if include_header else 0
```

Its physical length is computed from `ratio`:

```python
source_px      = bar_height / v_zoom
bp             = source_px / ratio
nice_bp        = _nice_scale_value(bp)
actual_height  = int(nice_bp * ratio * v_zoom)
```

For this math to produce a scale bar that visually matches the rendered reads, both `ratio` and `top_margin` **must be the values used when the source PNG was rendered**. See below.

---

## Calibration pitfalls

The animation engine is calibrated by three parameters:

| Parameter | Used for | Default |
|---|---|---|
| `--ratio` | bp ↔ pixel conversion in the source PNG and the live scale bar | `1/300` |
| `--top-margin` (through `scaled_top_margin`) | vertical position of the live scale bar within the viewport | derived from `args.top_margin × height_scale` |
| `--png-scale` (== `height_scale`) | scales `ratio` and `top_margin` for the rasterized PNG's pixel space | `8.0` |

When you run `KaryoScope_plot_reads.py --animate ...` end-to-end, all three values pass through automatically and the math is self-consistent.

**If you split the steps** — render a PNG with one set of arguments, then call `create_panning_animation.py` (or re-invoke `KaryoScope_plot_reads.py` against an existing PNG) with a different set — calibration breaks in visible ways:

| Symptom | Likely cause |
|---|---|
| Dynamic scale bar drawn **above** the sample labels instead of beside the reads | `top_margin × png_scale` doesn't match the actual header height in the source PNG. |
| "10 Kbp" tick spans far fewer (or far more) pixels than the embedded scale bar in the source PNG | `ratio × png_scale` is wrong. |
| Two `10 Kbp` markers visible side by side at frame 0 | `pan_start < left_margin`. Fixed by the `max(left_margin, …)` clamp; pre-fix renders show this. |
| Legend feature labels use different colors than the colors file says | The colors file path passed to the animation differs from the one used by the source PNG. The legend is regenerated from `--colors` regardless of what was baked into the PNG. |

**Recommended workflow:** drive everything from `KaryoScope_plot_reads.py --animate`. If you must split steps, log the exact `--ratio`, `--top-margin`, `--png-scale`, `--bar-width`, `--read-spacing`, and `--colors` values alongside the PNG so the animation step can pin them.

---

## Worked example: 16:4 wide-strip rendering at half PPT height

```bash
python scripts/KaryoScope_plot_reads.py \
  --bed "${BED_ARGS[@]}" \
  --colors KS_human_CHM13_v2.region_subtelomere_flat.colors.txt \
  --read-list all_rDNA_telogator_pass_reads.readlist.tsv \
  --output all_rDNA_telogator_pass_reads_16x4.png \
  --format png \
  --png-scale 2.0 \
  --bar-width 5 \
  --read-spacing 4 \
  --viewport-ratio none \
  --orient-satellite-top \
  --background black \
  --animate \
  --animate-zoom fixed \
  --animate-viewport 1920x480 \
  --animate-reads-per-second 5.15625
```

Notes:

- `--viewport-ratio none` is required when the source PNG's aspect ratio is unusual; otherwise the script will pad/crop to `16:9`.
- `--png-scale 2.0` produces a manageable source PNG; the animation engine receives `2.0` as `height_scale` and rescales `ratio` / `top_margin` to match.
- `--animate-zoom fixed` produces a constant-zoom pan. The 480-px viewport unmasks adaptive-zoom variation that the standard 1080-px viewport hides — `fixed` keeps the camera steady.
- `--animate-reads-per-second 5.15625` corresponds to 50 % of the default pacing.

---

## Troubleshooting

**The legend is missing a feature I expect to see.**
The legend only includes features present in the rendered reads. If the feature was filtered out by `--min-length`/`--max-length`/`--filter-group`, dropped during smoothing, or doesn't appear in any read at the chosen viewport, it won't be in the legend. Check the BEDs directly.

**The legend has two entries with the same color.**
Duplicate hex codes in the colors file. Audit and assign unique colors per feature, or merge the categories.

**The wrong subheading is shown.**
The subheading is taken verbatim from the `# ...` comment in the colors file directly above the feature. Edit the comment.

**Subheadings I expected aren't there at all.**
Either the colors file has no `#` comments, in which case `parse_sections` produces one `None`-headed section, or all the features under those headers were filtered out as not-in-use. The Section header itself is suppressed when the section is empty after filtering.

**The animation viewport zooms in and out as it pans.**
That's `--animate-zoom adaptive`. Switch to `--animate-zoom fixed` (the current default), or — at a tall viewport like 1080 px — the adaptive math is clamped by `max_zoom=1.0` and the variation isn't visible.

**Two scale bars visible at the start of the animation.**
You're on a build older than the `pan_start = max(left_margin, …)` fix in `create_panning_animation.py`. Update.
