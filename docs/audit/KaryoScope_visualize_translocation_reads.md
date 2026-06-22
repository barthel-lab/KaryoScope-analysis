# Audit: KaryoScope_visualize_translocation_reads.py (734 lines)

## 1. Purpose

FACT: Renders per-group SVG/PNG figures of individual translocation reads as
stacked multi-featureset horizontal tracks (one row per read, one sub-track per
featureset), with a scale bar and a color legend. Two input modes: TSV batch
(from stage-1 output) and direct read specs (module docstring 2-21; `main`
519-619).

ASSESSMENT: This is the **visualize** (stage 3) step and the only one of the trio
that draws. It already participates in the Phase-6–13 karyoplot push-down (it
imports `karyoplot.core.text` and `karyoplot.core.colors`).

## 2. CLI surface

FACT (`argparse`, lines 520-562):

- Mutually exclusive, required (527-534): `--input-tsv` (Path) OR
  `--reads SPEC...` where each spec = `READ_ID:SAMPLE:DATA_TYPE:REPLICATE:TRANS_TYPE`.
- `--results-dir` (required, Path); `--output-dir` (required, Path);
  `--colors-dir` (required) — dir with `{database}.{featureset}.colors.txt`.
- `--database` (default `KS_human_CHM13`).
- `--title` (custom; else auto-generated).
- `--featuresets` (default `chromosome subtelomeric acrocentric repeat region`).
- `--max-reads-per-file` (int, default 50 → split into parts).
- `--filter` (substring on group key, TSV mode only).
- `--parts` (int; only first N parts per group, TSV mode).
- `--log-file` / `--no-log-file` (`BooleanOptionalAction`, default True →
  `visualize_translocation_reads.log`).

FACT: `argparse`, not click.

## 3. Inputs & outputs (formats passed between the three tools)

FACT — TSV-mode input (`load_translocation_reads`, 190-204): reads stage-1 TSV
columns `read_id, sample, data_type, replicate, read_length, translocation_type`
(via `csv.DictReader`). It does NOT read the `_bp`/`_pct` columns stage-1 writes.

FACT — BED inputs differ by mode and use TWO naming conventions:
- TSV mode → `load_all_features_batch` (124-154) reads **translocation** BEDs via
  `get_translocation_bed_path` (92-109). Chromosome featureset:
  `{s}.{dt}.{rep}.{db}.chromosome.smoothed.{trans}.translocations.bed.gz`;
  other featuresets:
  `{s}.{dt}.{rep}.{db}.{trans}.{featureset}.smoothed.translocations.bed.gz`.
  NOTE these are `.smoothed.` files under `{results}/{s}/{dt}/{rep}/KaryoScope/{db}/`,
  which differ from the `.presmoothed.` files that stages 1 & 2 discover.
- direct-reads mode → `load_features_for_reads` (157-187) reads **standard**
  (non-translocation) BEDs via `get_standard_bed_path` (112-121):
  `{s}.{dt}.{rep}.{db}.{featureset}.smoothed.features.bed.gz`, filtered to the
  requested read_ids.
- Color files: `{database}.{featureset}.colors.txt` under `--colors-dir`,
  loaded by `karyoplot.core.colors.load_featureset_palettes`
  (`load_color_files`, 80-89).

FACT — BED row parse (143-152, 174-185): tab cols col0=read_id, col1=start,
col2=stop, col3=feature; ≥4 cols required.

FACT — Outputs (`render_group`, 467-492): writes an SVG, converts to PNG via
`rsvg-convert --zoom=...`, then **deletes the SVG** (line 492). Filenames:
`{s}.{dt}.{rep}.{trans}.n{N}.translocation_reads.{svg,png}` or, when split,
`...n{total}.part{p}of{total}.translocation_reads.*`.

INTER-TOOL FACT: TSV-mode consumes stage-1's TSV directly. It re-derives read
length from the TSV (`read_length`, line 201) in TSV mode, but in direct mode
**recomputes** length as max feature stop (705-714). It does not consume stage-2
outputs.

## 4. Pipeline / control flow (key functions + line numbers)

- `main` (519-619): validate dirs, banner, `load_color_files` (611), dispatch to
  `_run_tsv_mode` (615) or `_run_reads_mode` (617).
- `_run_tsv_mode` (622-679): load TSV, group by `(sample,dt,rep,trans)`, per
  group batch-load features (647), sort by length desc (652), chunk by
  `max_reads_per_file` (654-673) calling `render_group` per part, else single.
- `_run_reads_mode` (682-730): parse specs (`_parse_read_spec` 499-516), group by
  `(sample,dt,rep)`, `load_features_for_reads`, compute lengths, drop zero-length,
  one combined `render_group`.
- `render_group` (325-492): layout math, drawsvg `Drawing`, title, per-read rows
  (label + per-featureset tracks 415-447), scale bar (456-460), legend (462-465),
  save + rsvg PNG conversion.
- Drawing helpers: `draw_scale_bar` (233-247), `draw_legend` (250-318).
- Misc: `natural_sort_key` (207-213), `compute_scale_bar_bp` (216-226).

## 5. Key design decisions (cite lines)

FACT: Reads sorted by length descending before rendering/chunking (340, 652) so
parts are deterministic and longest reads appear first.
FACT: Coordinate ratio is per-group: `ratio = drawable_width / max_length`
(351-352); all reads share one bp→px scale within a figure.
FACT: Label color is chosen by trans_type substring — red `#FF4444` if
`"chr2_chr13" in trans_type` else blue `#60A5FA` (line 403). ASSESSMENT:
hard-coded two-translocation assumption.
FACT: Scale bar bp picked as largest "nice" value ≤ 25% of the smallest read
(`compute_scale_bar_bp`, 216-226).
FACT: PNG zoom scales with read count: `min(2.0 + (n-5)//5*0.5, 4.0)` (484).
FACT: SVG is transient — created then unlinked, only PNG kept (481-492).
FACT: Detection of translocation reads is again upstream — this tool only
visualizes reads it is told about (TSV rows or specs); no breakpoint logic.

## 6. Assumptions (checkable)

- TSV has the stage-1 column names exactly (KeyError otherwise; 197-202).
- Results layout is exactly `{results}/{sample}/{data_type}/{replicate}/KaryoScope/{database}/...`
  (101-102, 118-119).
- Two distinct BED naming conventions are correct and present (see §3); chromosome
  featureset uses a transposed `chromosome.smoothed.{trans}` order vs others
  (103-109).
- BED rows: ≥4 tab cols, col0=read_id, ints in col1/col2 (143-152).
- Color file format parseable by `load_featureset_palettes`; colors returned as
  `tuple` (value_format="tuple", 85) → legend handles tuple-vs-str (289-292) and
  track draw unpacks `(color, opacity)` (438-441), defaulting `("#ffffff",1.0)`.
- `rsvg-convert` on PATH for PNG output; absence is silently tolerated (489-490)
  but then NO output file remains (SVG already deleted, 492).
- Direct-read spec is exactly 5 colon-separated fields (504-509) — read IDs
  containing `:` would break parsing.

## 7. Dependencies

FACT: stdlib `argparse, gzip, csv, subprocess, sys, collections.defaultdict,
pathlib` (23-29) + **drawsvg** (line 31). FACT: karyoplot —
`karyoplot.core.text.abbreviate_read_name` (77) and
`karyoplot.core.colors.load_featureset_palettes` (82). FACT: external tool
**rsvg-convert** via subprocess (485-487). No samtools/bedtools, no pandas.
INTER-SCRIPT: consumes stage-1 TSV; no calls to other scripts.

## 8. Proposed home in new layout

ASSESSMENT:
- Subcommand: `karyoscope-analysis visualize-translocation-reads`.
- Thin wrapper `commands/visualize_translocation_reads.py`; feature-loading +
  grouping in `core/translocation_viz.py`; rendering split to drawing module.
- karyoplot push-down candidates (high value):
  - `compute_scale_bar_bp` (216-226) → `karyoplot.core.coords.pick_round_scale_bp`
    (already exists; current local logic uses a "% of smallest read" target
    rather than px-window — reconcile).
  - `draw_scale_bar` (233-247) label `f"{bp//1000} kb"` → already covered by
    `karyoplot.core.text.format_genomic_distance(style="kb_short")`; the line/tick
    drawing → `karyoplot.svg` (cluster_plot's scale bar was already routed there
    per recent commits).
  - `draw_legend` (250-318) overlaps `karyoplot.svg.legend.draw_grouped_legend` /
    `merge_by_color` / `strip_label_suffixes` — strong consolidation candidate.
  - SVG→PNG (483-492) → `karyoplot.svg.export.svg_to_png(scale=..., quiet=True)`
    (handles missing rsvg gracefully and keeps the SVG).
  - feature rasterization to pixels overlaps `karyoplot.svg.reads.*`.
  - BED reads → `karyoplot.core.io.load_bed`.
- `TeeLogger` (37-53) → shared `core/logging.py`.

## 9. Smells / risks / dead code / duplication

- DUP: `TeeLogger` (37-53) identical across the trio.
- DUP: local `draw_legend`/`draw_scale_bar`/`compute_scale_bar_bp` reimplement
  karyoplot.svg/core helpers the library already exposes (see §8).
- DUP: `natural_sort_key` (207-213) likely overlaps chromosome ordering helpers
  in `karyoplot.core.chromosomes`.
- RISK: `rsvg-convert` missing → `except ... pass` (489-490) then SVG is unlinked
  in `finally` (492) → run "succeeds" but produces NO image. Should fall back to
  keeping the SVG (karyoplot `svg_to_png` does the right thing).
- RISK: read IDs containing `:` break `_parse_read_spec` (504-509) which hard-
  requires exactly 5 fields.
- BUG-ish: track feature color lookup `featureset_colors[fs].get(...)` (438) will
  `KeyError` if `fs` absent from the colors dict (e.g. a featureset with no color
  file) — `.get(fs, {})` is used in the legend (270) but not here.
- INCONSISTENCY: two BED naming conventions and two length sources (TSV
  `read_length` vs recomputed max-stop) between the two modes (§3).
- SMELL: many magic layout constants (59-69) and inline hex colors (403, 411,
  432, 447, 453) — no theme integration (`karyoplot.core.theme`).
- SMELL: `font_family="Basic Sans"` hard-coded in every Text (e.g. 247, 391) —
  no `--font-family` flag (cf. cluster_plot's configurable FONT_FAMILY).
- SMELL: `render_group` mutates its `reads` argument via in-place `.sort` (340).

## 10. Testability notes

ASSESSMENT: The IO/grouping layer is testable (`load_translocation_reads`,
`load_all_features_batch`, `load_features_for_reads`, `natural_sort_key`,
`compute_scale_bar_bp`, `_parse_read_spec` are pure-ish). Rendering is harder:
`render_group` writes files and shells to rsvg; best tested via golden SVG (assert
on the in-memory `draw.Drawing` element list before save, or snapshot the SVG
prior to deletion). Push-down to karyoplot would shrink the untested surface.
Blockers: global `sys.stdout` patch (TeeLogger), `sys.exit` on errors, SVG
deletion makes artifact inspection awkward.

## 11. Open questions for the user

1. Is silently producing no output when `rsvg-convert` is missing acceptable, or
   should the SVG be kept as fallback (switch to `karyoplot.svg.export.svg_to_png`)?
2. Should the local legend/scale-bar/SVG-export be replaced by the existing
   `karyoplot.svg` helpers now, or is byte-identical output a requirement (as with
   the cluster_plot bench)?
3. The red/blue label color is hard-coded to `chr2_chr13` vs other (line 403) —
   should label color come from `--colors-dir` / a translocation palette instead?
4. Why two BED naming conventions (`.smoothed.…translocations` for TSV mode vs
   `.smoothed.features` for direct mode), and is the chromosome-featureset
   transposed name (103-106) intentional?
5. Direct-mode read IDs cannot contain `:`. Is that an acceptable constraint?
6. Should font family be configurable (`--font-family`) to match other plotting
   tools?
