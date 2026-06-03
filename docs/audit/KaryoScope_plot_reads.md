# Audit: KaryoScope_plot_reads.py (3620 lines)

## 1. Purpose
FACT: Standalone CLI tool that visualizes telomeric/satellite reads as stacked
colored bars, one feature-segment per BED interval. Each read becomes a vertical
bar (default) or horizontal bar (`--horizontal`); features within a read are
colored per a palette file. Module docstring lines 2-41. Supports three input
modes (KaryoScope results dir, direct BED files, cluster assignments), SVG/PNG
rendering, an optional categorical heatmap track, a two-tier sample/group label
system, an auto-filtered color legend, and panning-animation MP4 generation
(delegated to `create_panning_animation.py`).

ASSESSMENT: This is one of the two large "read raster" renderers in the repo
(the other being `KaryoScope_cluster_plot.py`). It is feature-dense and the
single `main()` (lines 3043-3616) carries a lot of orchestration weight.

## 2. CLI surface
FACT: `argparse` parser in `parse_args()` (lines 1655-2010). No subcommands;
all flags top-level. Groups:
- Input: `--samples`, `--results-dir`, `--bed`, `--clusters`, `--cluster-prefix`,
  `--colors` (required), `--read-list`, `--filter-group` (append),
  `--min-length`, `--max-length`, `--max-read-length` (note: two different
  length filters; see smells).
- Data path construction: `--analysis` (default `telogator`), `--database`
  (`KS_human_CHM13`), `--featureset` (`region`), `--smoothness` (`smoothed`).
- Output: `--output/-o` (required), `--format {svg,png,both}`, `--png-scale`
  (8.0), `--scale-bar-output`, `--batch-size`.
- Display: `--bar-width`, `--read-spacing`, `--sample-spacing`,
  `--subgroup-spacing`, `--ratio` (1/300), `--background {white,black,both}`,
  `--orient-telomere-top`, `--orient-chromosome-top`, `--orient-satellite-top`,
  `--markers`, `--marker-scale`, `--no-header`, `--no-scale-bar`,
  `--scale-bar-bp`, `--font-size` (11), `--viewport-ratio` (16:9), `--horizontal`,
  `--read-border`, `--legend`.
- Heatmap/labels: `--heatmap`, `--label-tier` (append, `COLUMN:DISPLAY`),
  `--heatmap-track` (append).
- Feature rendering: `--feature-mode {raw,transition,smooth}` (default smooth),
  `--min-feature-width` (0.5), `--min-width-exclude` (default `novel *arm* ct*`),
  `--oversample` (1).
- Animation: `--animate`, `--animate-direction`, `--animate-duration`,
  `--animate-reads-per-second` (20.625), `--animate-zoom`,
  `--animate-crop-ratio`, `--animate-fps`, `--animate-viewport`,
  `--animate-legend`, `--animate-include-header` (BooleanOptionalAction).

ASSESSMENT: ~50 flags on one flat parser. `--font-size` default is 11 in argparse
but many drawing functions default to 14 via `config.get("font_size", 14)`
(e.g. lines 929, 991, 1334, 1392) — the 14 fallbacks are dead since `font_size`
is always set in config (line 3390). `--marker-scale` help/usage is not in the
module docstring. `--animate-crop-ratio` and `--animate-zoom` are parsed but only
`--animate-zoom` is consulted indirectly; `--animate-crop-ratio` is never read
(dead flag — see smells).

## 3. Inputs & outputs
FACT inputs:
- Colors file via `load_color_mapping()` (line 2013), delegating to
  `karyoplot.core.colors.load_palette_file(parse_sections=True,
  suffix_both_ways=True, initial={"novel": "#ffffff"})`.
- BED features: `{results_dir}/{sample}/{analysis}/1/KaryoScope/{database}/
  {sample}.{analysis}.1.{database}.{featureset}.{smoothness}.features.bed`
  (`load_sample_bed_data` line 2198, builds path lines 2209-2215; tries `.gz`
  fallback lines 2218-2224). Direct BED via `load_bed_files_direct` (line 2231;
  `SAMPLE:PATH` or `PATH`). Cluster mode via `load_cluster_reads` (line 1478)
  reading `{prefix}.sequence_assignments.tsv` (uses pandas, imported locally
  line 1485).
- Read-list TSV via `load_read_list` (line 2028): col1 = read IDs, header row
  detected by first field in `{sequence,read,read_id}`, cols 2+ become metadata.
- Markers TSV (`--markers`): read_id, start, end, parsed in `main()` lines
  3344-3367.

FACT outputs:
- SVG via `drawsvg` (`draw_reads_svg` line 2785, `draw_reads_horizontal_svg`
  line 878).
- PNG via PIL+numpy (`draw_reads_png` line 219, `draw_reads_horizontal_png`
  line 615); `_save_png` (line 113) alpha-composites onto solid bg.
- Optional separate scale-bar SVG (`draw_scale_bar_svg` line 2713) which also
  shells out to `rsvg-convert` for a PNG.
- Optional separate legend PNG, scale-bar PNG, and MP4 for animation.
- Always writes a `.log` file next to `--output` (`setup_logging` line 62).
- `--background both` produces `_black`/`_white` variants; `--batch-size`
  produces `.batchNNN.*` files.

## 4. Pipeline / control flow (key functions + line numbers)
FACT: `main()` (3043):
1. Hardcode margins `top=30,left=30,bottom=15` (3048-3050); `setup_logging` (3053).
2. `--animate` forces `format=both` if svg (3062-3064).
3. Determine `input_mode` from flags (3066-3079).
4. `load_color_mapping` (3083).
5. Load reads by mode: `load_bed_files_direct` / `load_cluster_reads` /
   `load_sample_bed_data` (3087-3106).
6. Length filter (`--min/--max-length`) 3111-3117.
7. Read-list handling 3128-3244: `load_read_list`, parse `--label-tier`/
   `--heatmap-track`, build `read_groups`, `group_subgroup_order`,
   `read_metadata`/`metadata_columns`, `--filter-group`.
8. Sort `group_subgroup_order` by first-seen rank (3247-3255).
9. Heatmap validation/disable for horizontal (3258-3266).
10. Orientation: `orient_telomere_top` / `orient_chromosome_top` /
    `orient_satellite_top` (3274-3286) → `flipped_reads`.
11. `apply_read_list_grouping` (3293) relabels samples; builds `label_tiers`,
    `tier_display_names`, top-margin bump for extra tiers (3300-3312).
12. Heatmap color assign `assign_heatmap_colors` (3328); `--max-read-length`
    filter (3335-3341); markers load + flip (3344-3367).
13. `sort_reads` (3371) by sample order then descending length.
14. Build `config` dict (3374-3411).
15. Viewport-ratio scaling: recompute `ratio` only, reserve heatmap top-margin
    (3416-3466).
16. Per-bg-theme loop (3480-3612): batch or single; dispatch to the four
    `draw_reads_*` renderers; legend (composite or save); scale-bar PNG;
    `_run_animation` (once); separate scale-bar SVG.

FACT key helpers: `rasterize_features` (imported from karyoplot, used by all 4
renderers), `_build_colored_features` (155), `_draw_rect_rgba` (132),
`calculate_png_params`/`calculate_horizontal_png_params` (168/568),
`compute_legend_layout` (1210), `_filter_legend_features` (1165),
`compute_group_spans` (2126), `_orient_reads`/`_orient_reads_with_fallback`
(2358/2414), `draw_heatmap_grid_svg`/`_png` (2535/2580), `draw_scale_bar`/
`_svg`/`_png` (2669/2713/1528), `_run_animation` (1584).

## 5. Key design decisions (cite lines, include WHY if stated)
FACT:
- Reads are 4-tuples `(sample, read_id, read_length, features)` where
  `features=[(start,end,name),...]`; this shape threads through every function.
- Direct-to-PNG raster via numpy `_draw_rect_rgba` (132) "bypasses SVG element
  limits and produces images suitable for high-zoom panning animations"
  (docstring 222-224). WHY stated.
- PNG uses uniform height/width scaling so PNG is "proportionally faithful copy
  of the SVG, capped by MAX_PNG_DIMENSION" = 32000 (lines 110, 169-216).
  Comment line 110: "Stay under common image library limits."
- `feature_mode` default `smooth` (windowed majority-vote), with `transition`
  (min-width enforcement) and `raw` (integer pixel stacking) (parser help
  1916-1919). `min_width_exclude` default `novel *arm* ct*` exempts those from
  min-width inflation (1929-1932).
- Legend layout fills column-first, sections start new columns, and `num_rows`
  is grown to avoid section overlap (detailed comment 1267-1287) — WHY stated:
  "otherwise the later sections wrap back to col=0 and overwrite earlier content."
- Per-sample horizontal padding to prevent label overlap for small samples
  (comment 264-266, lines 271-282). WHY stated.
- Viewport-ratio only changes `ratio` (px/bp), keeping bar_width/spacing/heatmap
  constant (comment 3414-3415). WHY stated.
- `--filter-group` and group-order building deliberately scan `read_rows`
  (raw file order) not the `read_groups` dict because duplicate read IDs
  overwrite dict entries (comments 3166-3168, 3216-3217). WHY stated.
- Animation auto-direction: vertical reads → horizontal pan, horizontal reads →
  vertical pan (line 1602).

## 6. Assumptions (checkable statements)
FACT:
- Results-dir layout is exactly `{sample}/{analysis}/1/KaryoScope/{database}/...`
  with the `1` literal (line 2212). Hardcoded.
- BED is tab-separated with ≥4 cols (id,start,end,name); lines with <4 cols
  silently skipped (`_parse_bed_file` 2182).
- A read's length = `max(end for features)` (lines 2191, 1520); assumes features
  cover the read end.
- Read-list header detected only if first field lowercased ∈
  `{sequence,read,read_id}` (line 2059); otherwise cols never captured.
- Cluster assignments TSV has `cluster`,`sample`,`sequence` (or `read`) columns;
  cluster IDs are ints (`int(c)` line 3091; `int()` on `parts[1]/[2]` for
  markers 3352).
- Font "Basic Sans" availability assumed in SVG (hardcoded `font_family =
  "Basic Sans"` in every SVG draw, e.g. 932, 988, 1393, 2576, 2707, 2761);
  PNG uses `karyoplot.core.fonts.pil_font` which falls back if unavailable.
- `rsvg-convert` optional (try/except 2772-2782).
- Heatmap palette is fixed 6-color cycle + `#545454` missing (2512-2513).
- 25-entry hardcoded chromosome feature set incl. `_specific` variants
  (CHROMOSOME_FEATURES 2337-2347); satellite set 2350-2355; telomere set 2334.

## 7. Dependencies
FACT external libs: `drawsvg` (SVG), `PIL` (Image/ImageDraw/ImageFont),
`numpy`, `pandas` (local import in `load_cluster_reads` only, line 1485),
stdlib `argparse,fnmatch,gzip,logging,os,subprocess,sys,time,re,collections`.

FACT karyoplot usage (already delegated):
- `from karyoplot.core.fonts import pil_font as _load_font` (line 129).
- `from karyoplot.svg.reads import features_to_pixels_direct, rasterize_features,
  smooth_features_to_pixels` (lines 147-151) — used by all 4 renderers. The
  module docstring of `karyoplot/svg/reads.py` confirms these were ported from
  inline copies in BOTH `KaryoScope_cluster_plot.py` and this file (pure dedup).
- `load_color_mapping` (2013) delegates to `karyoplot.core.colors.
  load_palette_file`.

FACT NOT yet delegated (remains local):
- `hex_to_rgba`/`hex_to_rgb` (86/104) — duplicate of
  `karyoplot.core.colors.hex_to_rgba`/`hex_to_rgb` (colors.py:333/325).
- `_draw_rect_rgba` (132) — alpha-blend rect onto numpy array.
- Legend: `LegendItem`/`LegendLayout` namedtuples (1139-1142),
  `compute_legend_layout` (1210), `_filter_legend_features` (1165),
  `_build_heatmap_legend_items` (1145), `draw_legend_png`/`draw_legend_svg`
  (1077/1379), `composite_legend` (1435), `_estimate_*`/`estimate_*` height
  helpers (1329/1364). Note `karyoplot.svg.legend` already has
  `strip_label_suffixes` (legend.py:118) that duplicates the regex in
  `_filter_legend_features` (1177/1189/1198).
- Scale-bar drawing (PNG + SVG, 3 variants).
- Orientation logic (`_orient_reads`, fallback variant).
- Heatmap grid/legend drawing (4 functions).

FACT inter-script / shared-file deps:
- `_run_animation` (1584) does `sys.path.insert(0, script_dir)` then
  `from create_panning_animation import create_horizontal_panning,
  create_vertical_panning, create_adaptive_horizontal_panning` (1591-1594).
  This is a runtime import of a sibling script by path.
- Shares the read-rasterization stack with `KaryoScope_cluster_plot.py`
  (both import the same 3 functions from `karyoplot.svg.reads`).
- No `_feature_vocab` import found in this file.
- External tool: `rsvg-convert` via subprocess (2773).

## 8. Proposed home in new layout
ASSESSMENT:
- Subcommand: `karyoscope-analysis plot-reads`.
- `commands/plot_reads.py`: thin click command — define options, call into core.
- `core/plot_reads/` (or `core/reads_render.py`): the orchestration currently in
  `main()` (input-mode dispatch, read-list/grouping/heatmap setup, viewport
  scaling, per-theme loop) split into pure functions: `build_render_config()`,
  `prepare_reads()` (filter/orient/group/sort), `apply_viewport_ratio()`.
- `core/io/reads_bed.py`: `_parse_bed_file`, `load_sample_bed_data`,
  `load_bed_files_direct`, `load_cluster_reads`, `load_read_list`.
- `core/orientation.py`: `_orient_reads*`, `orient_*`, and the
  TELOMERE/CHROMOSOME/SATELLITE feature sets (these are domain vocab and may
  belong in a shared analysis constants module, possibly `core/_feature_vocab`).
- Render functions (`draw_reads_*`, scale-bar, heatmap) → see push-down below.

FACT/ASSESSMENT karyoplot push-down candidates + module:
- `hex_to_rgba`/`hex_to_rgb` → DELETE locally, import from
  `karyoplot.core.colors` (already exists).
- `_draw_rect_rgba` and the numpy raster bar-drawing loops → `karyoplot.svg.reads`
  or a new `karyoplot.mpl`/raster helper (it is the PIL/numpy counterpart to the
  already-extracted `rasterize_features`).
- Legend stack (`compute_legend_layout`, `LegendItem/Layout`,
  `_filter_legend_features`, `draw_legend_png/svg`, `composite_legend`,
  height estimators) → `karyoplot.svg.legend` / `karyoplot.mpl.legend`. This is
  the single biggest extraction; `karyoplot.svg.legend` already hosts grouped/
  column legend primitives and `strip_label_suffixes`.
- Scale-bar drawing → `karyoplot.svg` (SVG) and a PNG counterpart; per recent
  commits cluster_plot already routes scale-bar bp-pick/label through
  `karyoplot.core.coords.pick_round_scale_bp` + `core.text.format_genomic_distance`
  — this script does NOT yet (hardcodes "X Kbp" via `// 1000`).
- Heatmap grid/legend (SVG+PNG) → `karyoplot.svg.tracks` / `karyoplot.mpl.heatmap`.

## 9. Smells / risks / dead code / duplication (line-cited)
FACT/ASSESSMENT:
- Two overlapping length filters: `--max-length` (3111-3117) and
  `--max-read-length` (3335-3341). Confusing; `--max-read-length` is applied
  after grouping. Likely consolidate.
- `font_size` 14 fallbacks (`config.get("font_size", 14)`) are dead — config
  always sets it to `args.font_size` default 11 (line 3390). Lines 929, 991,
  1334, 1392, 2547, etc.
- `--animate-crop-ratio` parsed (1977-1982) but never read anywhere.
- Marker-arrowhead PNG loop unpacks with a trailing comma: `for mx, my, in
  pending_markers:` (line 412) — works but is a typo-style smell.
- Massive duplication between the 4 renderers (`draw_reads_svg`/`_png` and
  horizontal variants): identical feature-mode dispatch, `_build_colored_features`
  + `rasterize_features` blocks, read-border, scale-bar, and two-tier label
  logic appear 2-4× (e.g. label tiers SVG 2948-3023 vs PNG 490-561 vs horizontal
  762-838/996-1068). Each is ~150-300 lines of near-parallel code.
- `hex_to_rgba`/`hex_to_rgb` duplicate karyoplot (7).
- `_load_font(font_size)` is imported but PNG also recomputes/reloads fonts
  repeatedly inside loops (e.g. 268, 308, 448, 637, 2637, 2648) — perf and
  clarity smell.
- `draw_legend_svg` recomputes `compute_legend_layout` AND
  `estimate_featureset_legend_height` recomputes it again for the same data
  (1374 vs 1402) — double layout pass per render.
- Two legend-height estimators exist: `_estimate_heatmap_legend_height` (1329)
  appears UNUSED (only `estimate_featureset_legend_height` is called, line 2845;
  heatmap height folds into it via `extra_items`). Likely dead code.
- `draw_ctx = ImageDraw.Draw(img)` is reassigned but unused right after the numpy
  array build in `draw_reads_png` (304, 407) — first assignment at 304 is dead
  (overwritten at 407 after array round-trip).
- `subprocess` import only for `rsvg-convert` (2773); silent if missing.
- Sibling-script import via `sys.path` mutation (1590) is fragile under packaging.
- Hardcoded chromosome/satellite/telomere vocab (2334-2355) will silently fail
  on non-human / renamed features.
- `import sys as _sys` shadow inside `_run_animation` (1587) when `sys` already
  imported at module top.
- `compute_group_spans` and the `_label` composite-label helper are defined
  3× (2097, 2139, 3296) with identical `f"{g} — {s}"` logic.

## 10. Testability notes
FACT pure / easily unit-testable (no I/O):
- `hex_to_rgba`/`hex_to_rgb` (86/104), `_build_colored_features` (155),
  `compute_legend_layout` (1210), `_filter_legend_features` (1165),
  `_build_heatmap_legend_items` (1145), `compute_group_spans` (2126),
  `apply_read_list_grouping` (2081), `sort_reads` (2265),
  `assign_heatmap_colors` (2506), `_orient_reads`/`_orient_reads_with_fallback`
  (2358/2414), `calculate_png_params`/`calculate_horizontal_png_params`
  (168/568), `compute_viewport_params` (2284, currently UNUSED — see note).
- `compute_viewport_params` (2284) is defined but never called by `main()`
  (viewport logic is inlined 3416-3466) — dead function but a clean unit target
  if revived.

FACT I/O-bound (need fixtures / tmp files):
- `_parse_bed_file`, `load_sample_bed_data`, `load_bed_files_direct`,
  `load_cluster_reads` (pandas), `load_read_list`, `load_color_mapping`.
  Small synthetic BED/TSV/colors fixtures suffice.

ASSESSMENT integration-only (golden-image): the four `draw_reads_*` renderers,
heatmap/scale-bar/legend draws, `_run_animation`. Recent commits reference an
existing "fiberseq bench" / "18 byte-identical" SVG golden harness for
cluster_plot; this file should adopt the same byte-identical SVG golden approach
for `draw_reads_svg`/`draw_reads_horizontal_svg`, and image-diff (or dimension)
checks for PNG.

## 11. Open questions for the user
1. `--max-length` vs `--max-read-length`: intentional distinct semantics, or
   should they be merged into one `--max-length` in the new CLI?
2. `--animate-crop-ratio` is dead — drop it, or wire it through to
   `create_panning_animation`?
3. Should the hardcoded human chromosome/satellite/telomere feature vocab
   (2334-2355) move into a shared `core/_feature_vocab` (and become
   config/database-driven for non-human references)?
4. The legend layout engine here is independent of `karyoplot.svg.legend`'s
   grouped/column legends. Do you want a single unified legend module in
   karyoplot, or keep the analysis-specific section-aware grid layout separate?
5. Scale-bar labels here are hardcoded "X Kbp"; cluster_plot already routes
   through `karyoplot.core.coords.pick_round_scale_bp` + `core.text.
   format_genomic_distance`. Confirm we should converge this script onto those
   helpers (will change SVG output bytes / golden baselines).
6. `_run_animation` imports the sibling script by `sys.path`. In the package,
   should `create_panning_animation` become `core/animation.py` imported
   directly, with the animation primitives pushed to `karyoplot`?
7. Is `compute_viewport_params` (2284) intended to be used (currently dead while
   `main()` inlines its own viewport math)?
