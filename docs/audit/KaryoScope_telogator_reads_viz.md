# Audit: KaryoScope_telogator_reads_viz.py (856 lines)

## 1. Purpose
Visualizes telomeric reads (from Telogator output) as side-by-side vertical
colored bars. Each read is a thin vertical bar; coloured segments along the bar
represent region/satellite features (e.g. telomere, subtelomere, satellite
classes). Reads are grouped by sample and sorted longest→shortest within each
sample (module docstring lines 2-14; `sort_reads` lines 412-428). Produces an
SVG and a matching PNG (the PNG bypasses SVG element limits for very large read
sets), plus an optional standalone scale-bar SVG/PNG.

## 2. CLI surface
`argparse` parser in `parse_args` (lines 166-288). No subcommands; single
script entry via `main` (lines 756-852) guarded by `__main__` (855-856).

Input selection (mutually-exclusive in practice, validated lines 776-784):
- `--samples` (nargs+, line 172) + `--results-dir` (line 178) — auto-build BED
  paths.
- `--bed` (nargs+, line 183) — explicit `SAMPLE:PATH` or bare `PATH` specs.

Required: `--colors` (line 190), `--output/-o` (line 195).
Optional output: `--scale-bar-output` (line 198).

Data layout knobs: `--database` (default `KS_human_CHM13`, line 205),
`--featureset` (default `region`, 210), `--smoothness` (default `smoothed`, 215).

Display: `--bar-width` (3, 222), `--read-spacing` (3, 228), `--sample-spacing`
(20, 234), `--ratio` (1/300 bp→px, 240), `--background` {white,black} (black,
246), `--top-margin` (80, 252), `--left-margin` (60, 258), `--bottom-margin`
(40, 264), `--orient-telomere-top` (flag, 270), `--title` (276),
`--log-file/--no-log-file` (BooleanOptionalAction, default True, 280).

## 3. Inputs & outputs
Inputs:
- Region-feature BED files. Sample mode builds the path
  `{results_dir}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed` (`load_sample_bed_data` lines 316-322), trying `.gz` if the plain file is absent (325-331).
- BED columns used: col0=read_id, col1=start, col2=end, col3=feature
  (lines 342-345, 397-399). **Note the non-standard layout** — column 0 is the
  read id, not a chromosome (see Assumptions).
- A colors palette file (e.g. `KS_human_CHM13.region.colors.txt`).

Outputs:
- `--output` SVG (`draw_reads_svg`, line 749 save).
- Sibling PNG at `{output sans ext}.png` (`draw_reads_png`, always written,
  lines 838-840).
- `{output sans .svg}.log` tee'd console log when `--log-file` (lines 760-765).
- Optional `--scale-bar-output` SVG plus sibling `.png` via `rsvg-convert`
  (`draw_scale_bar_svg`, lines 526-593).

## 4. Pipeline / control flow (key functions + line numbers)
- `main` (756): logging setup → validate inputs (776-784) → load colors (788) →
  load BED (792-804) → optional orient (813-815) → sort (819) → SVG (835) →
  PNG (840) → optional scale-bar file (843-850).
- `load_color_mapping` (291-302): delegates to
  `karyoplot.core.colors.load_palette_file(suffix_both_ways=True, initial={"novel":"#ffffff"})`.
- `load_sample_bed_data` (305-355): builds path, gzip fallback, groups records
  by read id into `read_features`, read_length = max end (350).
- `load_bed_files_direct` (358-409): parses `SAMPLE:PATH`/`PATH` specs, returns
  `(all_reads, sample_order)`. **Body is a near-verbatim duplicate** of the
  parse/group/length loop in `load_sample_bed_data` (compare 334-353 vs 388-407).
- `sort_reads` (412-428): sort by sample rank then descending length.
- `orient_telomere_top` (435-480): flips reads whose telomere features
  (`TELOMERE_FEATURES`, line 432) sit in the second half so telomere is at top.
- `draw_reads_svg` (596-753): computes dimensions, draws bg, optional inline
  scale bar, per-read feature rectangles, per-sample top line + repeated labels
  (every 300 reads, line 688), inter-sample separators; returns image_height.
- `draw_reads_png` (55-163): Pillow re-implementation of the SVG renderer for
  large outputs.
- `draw_scale_bar` (483-523) inline SVG scale bar; `draw_scale_bar_svg`
  (526-593) standalone scale-bar SVG + rsvg PNG.
- Helpers: `TeeLogger` (30-46), `hex_to_rgb` (49-52).

## 5. Key design decisions (cite lines + WHY if stated)
- Dual SVG + PNG rendering. Stated WHY: PNG "handles large images better than
  SVG" / "bypasses SVG element limits" (lines 56, 837). Each read can be 1000s
  of rectangles, so a wide cohort blows past rsvg-convert/SVG limits.
- Separate scale-bar file option so the main image's left margin can shrink to
  20 px (`effective_left_margin`, lines 78, 633, 652) and the bar can be
  composited downstream (e.g. by the animation script).
- Labels repeated every 300 reads (`label_interval`, lines 131, 688) so the
  sample name stays visible while panning a very wide image. WHY tied to the
  animation use case (not stated in-file but implied by repetition logic).
- `--orient-telomere-top` normalises read orientation so the telomere end is
  always at y=0, making cross-read comparison meaningful (docstring 435-446).
- Colors loaded with `suffix_both_ways=True` and a `novel→#ffffff` seed
  (lines 298-302); unknown features fall back to white (`colors.get(feature,
  "#ffffff")`, lines 117, 675).
- Font: PNG path uses `karyoplot.core.fonts.pil_font(16)` (lines 94-95) with a
  Barthel-first fallback; SVG path hardcodes `font_family="Basic Sans"`
  (lines 689, 574) as a string only.

## 6. Assumptions (checkable statements; telogator output format)
- BED column 0 is a **read id**, not a chromosome; reads are grouped by col0
  (lines 342-346). This is a Telogator-specific feature BED, not a standard
  genome BED. Standard BED tools / `karyoplot.core.io.iter_bed_records` assume
  col0=chrom, so the column *semantics* differ even though the column
  *positions* match.
- Coordinates are per-read offsets starting near 0; read length = max feature
  end (lines 350, 404). Assumes features tile the read with no large trailing
  unannotated region beyond the last feature.
- The fixed sample-mode path template (lines 316-322) assumes the Telogator
  KaryoScope output tree `telogator/1/KaryoScope/{db}/...` exactly. Any layout
  change silently warns and skips the sample (lines 330, 384).
- Telomere features are exactly `{canonical_telomere, noncanonical_telomere}`
  (line 432); other telomere naming would defeat orientation.
- `--ratio` default 1/300 assumes downstream consumers use the same bp→px scale
  (the animation script's `--ratio` default also 1/300, kept in sync manually).
- `rsvg-convert` is on PATH for the optional scale-bar PNG (lines 586-593);
  failure is swallowed (592-593).

## 7. Dependencies
External libs: `drawsvg` (SVG, line 26), `Pillow` (`Image, ImageDraw,
ImageFont`, line 27 — `ImageFont` imported but unused here since fonts come from
karyoplot), `gzip`, `subprocess`, `argparse`, `os`, `sys`, `collections`.

External tools: `rsvg-convert` (scale-bar PNG only, lines 586-591).

karyoplot usage (already delegated):
- `karyoplot.core.fonts.pil_font` (lines 94-95).
- `karyoplot.core.colors.load_palette_file` (lines 297-298).

Inter-script / shared deps: this script's SVG/PNG output is the **input** to
`create_panning_animation.py`; the two share, by convention only, `--ratio`
(1/300), `--top-margin` (80), `--left-margin` (60), and the "labels repeated so
panning shows the sample name" idea. It is a sibling of
`KaryoScope_plot_reads.py` (which has its own, more elaborate read renderer in
`karyoplot.svg.reads`); this telogator viz does **not** use `karyoplot.svg.reads`
and reimplements its own simpler rectangle renderer twice (SVG + PNG).

## 8. Proposed home in new layout
- Subcommand: `karyoscope-analysis telogator-reads-viz` (thin
  `commands/telogator_reads_viz.py` → click command parsing args, delegating to
  `core/telogator_reads_viz.py`).
- Decomposition:
  - `core/io/telogator_bed.py`: BED loading/grouping (`load_sample_bed_data`,
    `load_bed_files_direct` collapsed into one loader + a path-builder).
  - `core/telogator_reads_viz.py`: `sort_reads`, `orient_telomere_top`, config
    assembly, render dispatch.
  - Renderers either stay local under `core/render/` or push to karyoplot.
- karyoplot push-down candidates + module:
  - Local `hex_to_rgb` (49-52) → already in `karyoplot.core.colors.hex_to_rgb`
    (delete the local copy).
  - BED parse loops → `karyoplot.core.io.iter_bed_records` (already streams
    `(c0, start, end, name, *rest)`); the per-read grouping is the only
    script-specific layer.
  - `TeeLogger` (30-46) → shared logging util (mirrors a pattern likely present
    in other analysis scripts; candidate for `karyoscope_analysis.core.logging`
    rather than karyoplot, since it's app-level not plot-level).
  - Vertical scale-bar drawing (`draw_scale_bar`, `draw_scale_bar_svg`) →
    `karyoplot.svg` (a `scale_bar` helper). The animation script has *three*
    more scale-bar variants (see other audit); consolidate all four.
  - rsvg PNG render (lines 586-591) → `karyoplot.svg.export.svg_to_png`
    (already exists; also gives `is_rsvg_convert_available`).
  - The bar/feature-rectangle renderer (SVG + PNG twins) → a `karyoplot.svg` /
    raster reads helper, ideally reconciled with `karyoplot.svg.reads`.

## 9. Smells / risks / dead code / duplication (line-cited)
- **Duplicated BED loader**: `load_sample_bed_data` (334-353) vs
  `load_bed_files_direct` (388-407) share the same parse/group/length body.
- **Duplicated renderer**: `draw_reads_svg` (596-753) and `draw_reads_png`
  (55-163) maintain two copies of layout math (dimensions 79-85 vs 634-640,
  positioning loops, label intervals). Easy to drift; e.g. PNG label loop
  (147-149) and SVG label loop (714-730) differ in y-offset (`-25` vs `-12`).
- **Duplicated scale-bar code**: `draw_scale_bar` (483-523) and the body of
  `draw_scale_bar_svg` (546-577) repeat the same bar/tick/label drawing.
- **Local `hex_to_rgb`** (49-52) duplicates `karyoplot.core.colors.hex_to_rgb`.
- **PIL pixel limit**: PNG path never sets `Image.MAX_IMAGE_PIXELS = None`
  (the animation script does, line 57). A very wide cohort can trip Pillow's
  `DecompressionBombError` on `Image.new`/`save` (lines 90, 160). Risk for the
  exact "large image" case this PNG path exists to serve.
- **Unbounded image width**: `image_width` (79-84) scales linearly with total
  reads with no cap or warning; combined with the pixel-limit gap above this can
  fail late after all parsing.
- **PNG separator alpha is a no-op**: `fill=(*text_color, 77)` on an RGB image
  (line 157, comment "0.3 opacity approximation") — RGB has no alpha channel, so
  the 77 is ignored and the line draws fully opaque (cosmetic only).
- **`TeeLogger.log` file never closed**: `sys.stdout` is reassigned (line 765)
  but `.close()` (45-46) is never called, and stdout is not restored.
- **Silent skips**: missing BED files only print a warning and `continue`
  (lines 330-331, 384); if *all* samples miss, the "No reads" error (808-810)
  returns 0 exit code (uses `return`, not `sys.exit(1)`), so failures look like
  success to a pipeline.
- **Unused params**: `left_margin`/`bottom_margin` pulled from config in
  `draw_reads_png` (60-61) but `left_margin` is overridden by the hardcoded
  `effective_left_margin = 20` (line 78); the configurable `--left-margin` has
  no effect on the PNG.
- `ImageFont` import (line 27) unused.

## 10. Testability notes
- Pure/easily unit-testable: `hex_to_rgb`, `sort_reads`, `orient_telomere_top`
  (deterministic transforms on tuples), path construction in
  `load_sample_bed_data`.
- BED loaders are testable with tiny fixture BEDs (incl. a `.gz`) — good golden
  candidates for the grouped `all_reads` structure.
- Rendering is hard to unit-test directly; recommend golden-image (PNG hash) or
  golden-SVG (string) tests on a 2-3 read fixture, matching the repo's existing
  "byte-identical bench" approach used in Phase 13 commits.
- `rsvg-convert` and font availability must be stubbed/skipped in CI; the
  Barthel fonts live under `~/Documents/...` and won't exist on CI runners
  (`pil_font` degrades gracefully, but SVG hardcodes "Basic Sans" — golden SVG
  is environment-independent there, golden PNG is font-dependent).
- `TeeLogger` + `sys.stdout` reassignment makes capturing output in tests
  awkward; moving logging behind the standard `logging` module (per
  gold-standard CLI) would fix this.

## 11. Open questions for the user
- Should the telogator feature-BED reader be unified with the generic
  `karyoplot.core.io` BED helpers despite the col0=read_id semantic difference,
  or kept as a distinct `telogator_bed` reader?
- Is the dual SVG+PNG output still required, or can we standardise on PNG (with
  `MAX_IMAGE_PIXELS=None`) for large cohorts and SVG only below a read-count
  threshold?
- Should this viz reuse `karyoplot.svg.reads` (the richer renderer used by
  `plot_reads`) instead of its own rectangle renderer, unifying the two read
  visualisations?
- Is the hardcoded Telogator path template (`telogator/1/KaryoScope/...`)
  stable, or should it become a configurable pattern?
- Should missing-input / no-reads conditions exit non-zero for pipeline safety?
