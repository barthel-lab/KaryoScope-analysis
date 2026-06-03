# Audit: KaryoScope_draw_legend.py (195 lines)

## 1. Purpose
FACT: Docstring (lines 2-23) states it generates "SVG legends from KaryoScope
color files" and describes itself as a "Thin CLI wrapper around
:func:`karyoplot.svg.legend.make_legend_drawing`" (line 4). It loads a TSV
`(feature, color)` colors file, optionally filters / merges-by-color / groups the
entries, builds a `Theme`, and calls `make_legend_drawing` to produce and save an
SVG (`main`, lines 139-191).

ASSESSMENT: The rendering is genuinely delegated to `karyoplot`, but the script is
NOT yet "thin": ~90 lines of CLI-input parsing/transform logic (`load_colors_file`,
`filter_items`, `parse_groups`, `group_items`, plus inline merge-override and theme
construction) live in the script and are not in the library (see section 8).

## 2. CLI surface
FACT: `argparse`-based (`parse_args`, lines 35-79), `RawTextHelpFormatter`. Options:
- `--colors` (required) — TSV colors file path (line 41).
- `--output` (default `legend.svg`) — output SVG path (line 43).
- `--columns` (int, default `None`) — column count, auto if unset (line 47).
- `--rows` (int, default `None`) — row count, auto if unset (line 49).
- `--include` / `--exclude` (default `None`) — comma-separated feature filters
  (lines 53-56).
- `--merge-same-color` (flag) — collapse features sharing a hex color (line 59).
- `--merge-label` (default `None`) — `feature=label,...` override map (line 61).
- `--groups` (default `None`) — `'Header:feat,feat;Header2:...'` category headers,
  prefix-matched (lines 65-67).
- Styling: `--swatch-size` (8), `--font-size` (12), `--background` (`#000000`),
  `--text-color` (`#FFFFFF`), `--stroke-color` (`#FFFFFF`), `--row-spacing` (14),
  `--col-spacing` (`None`), `--padding` (15) (lines 70-77).

ASSESSMENT: Surface is large (16 options) but maps cleanly onto
`make_legend_drawing`'s keyword parameters plus a handful of pre-processing toggles
— a good candidate for a click subcommand.

## 3. Inputs & outputs
FACT:
- Input: a TSV/whitespace-delimited colors file (`load_colors_file`, lines 82-91).
  Each line split on whitespace; lines with `< 2` fields or whose first field
  lowercased equals `"feature"` (header) are skipped (line 88); first two fields
  taken as `(feature, color)` (line 90). Example files exist in `scripts/`
  (`KS_human_CHM13.chromosome_acrocentric.colors.txt`).
- Output: an SVG written via `d.save_svg(args.output)` (line 189).
- stdout: progress prints — loading path (line 142), entry counts (lines 144, 148,
  159), header count (line 166), saved path + pixel dimensions (lines 190-191).
- No logging module; uses bare `print`.

## 4. Pipeline / control flow / public API
FACT (functions + line numbers):
- `parse_args()` — 35-79. Builds and returns the argparse namespace.
- `load_colors_file(filepath) -> list[tuple[str,str]]` — 82-91. Parses the colors
  TSV (see section 3).
- `filter_items(items, include=None, exclude=None)` — 94-101. Keeps items whose
  feature is in the include set (line 97) and/or not in the exclude set (line 100).
  Exact string match, not prefix.
- `parse_groups(groups_str) -> [(header, [prefixes])]` — 104-114. Splits on `;`
  then `:`; groups without a `:` are skipped (line 109).
- `group_items(items, groups_str)` — 117-136. Assigns each item to the first group
  whose any prefix (case-insensitive `startswith`, line 127) matches; emits a
  header row `(header, "", True)` before each non-empty bucket (line 131);
  unassigned items appended at the end as non-headers (lines 133-135).
- `main()` — 139-191. Orchestrates: load -> filter -> (merge-by-color OR
  strip-suffixes) -> optional group -> build `Theme` -> `make_legend_drawing` ->
  `save_svg` -> print dims.

Control-flow note (FACT): the merge branch (lines 150-159) parses `--merge-label`
into an `overrides` dict, calls `merge_by_color`, and re-wraps results as 3-tuples
`(label, color, False)`; the else branch (lines 160-161) applies
`strip_label_suffixes` to every feature and wraps as 3-tuples. If `--groups` is set
(lines 163-166), items are flattened back to 2-tuples (line 164) and re-expanded by
`group_items`.

## 5. Key design decisions
FACT:
- Rendering fully delegated to `karyoplot.svg.legend.make_legend_drawing` (import
  lines 28-32; call lines 177-188). `merge_by_color` and `strip_label_suffixes` are
  ALSO imported from the library (lines 30-31) rather than reimplemented — so the
  color-merge and suffix-strip *algorithms* already live in karyoplot.
- A `Theme` is constructed inline (lines 169-175) "to match the legacy CLI args
  (background/text/font)" (comment line 168): `background=args.background`,
  `text=args.text_color`, and `line`/`muted_line` both set to `args.text_color`
  (lines 172-174) — i.e. legend swatch strokes/lines follow the text color unless
  `--stroke-color` overrides (passed separately at line 187).
- The `(label, color, is_header)` 3-tuple shape is the contract `make_legend_drawing`
  expects (verified in library: `make_legend_drawing` items are
  `list[tuple[str,str,bool]]`, plotlib legend.py line 298).

ASSESSMENT: This is a backward-compat shim: it preserves the legacy argparse flag
surface and the dark default theme while routing all drawing through the shared
library. No v1/v2 vocab reconciliation is involved (that concern belongs to
`_feature_vocab`, not this file).

## 6. Assumptions (checkable statements)
1. The colors file is whitespace-delimited with feature in field 0 and a color
   string in field 1; any header row literally starts with `feature` (case-
   insensitive) (line 88). Checkable against the sample `.colors.txt` files.
2. Color strings are valid SVG fill values (hex like `#RRGGBB` or named) — the
   script never validates them; they pass straight to the swatch `fill` (library
   line ~404). `merge_by_color` upper-cases colors as the grouping key (plotlib
   line 284), so case-insensitive hex grouping is assumed.
3. `--merge-label` values are well-formed `feature=label` pairs; pairs without `=`
   are silently skipped (line 154).
4. `--groups` headers never legitimately contain `:` beyond the first (split with
   `maxsplit=1`, line 111) and feature lists never contain `,` within a name.
5. Group matching is by case-insensitive prefix (line 127), so a prefix can match
   multiple features and an item lands in the FIRST matching group only (`assigned`
   set, lines 125, 129).
6. `make_legend_drawing` and the `Theme` 4-field constructor signatures match what
   is called (lines 169-175, 177-188) — i.e. the installed `karyoplot` version is
   compatible.

## 7. Dependencies
FACT:
- External libs: `argparse` (stdlib, line 25) only.
- karyoplot usage (lines 27-32):
  - `karyoplot.core.theme.Theme` (line 27) — constructed at lines 169-175.
  - `karyoplot.svg.legend.make_legend_drawing` (line 28) — called lines 177-188.
  - `karyoplot.svg.legend.merge_by_color` (line 30) — called line 157.
  - `karyoplot.svg.legend.strip_label_suffixes` (line 31) — called line 161.
  - (Indirectly, `make_legend_drawing` returns a `drawsvg.Drawing`; `.save_svg`
    and `.width`/`.height` used at lines 189-191.)
- Who imports this: NOBODY. Grep of `scripts/` shows no other file imports
  `KaryoScope_draw_legend`; it is a standalone entrypoint (`if __name__ ==
  "__main__": main()`, lines 194-195).
- External tools: none.

## 8. Proposed home in new layout
RECOMMENDATION: becomes the click subcommand `karyoscope-analysis draw-legend`.

Is it "thin enough to become a trivial subcommand"? PARTIALLY — not yet.
- The DRAWING is already delegated (good). But the script still owns four
  input-shaping helpers that are NOT in karyoplot:
  `load_colors_file` (82-91), `filter_items` (94-101), `parse_groups` (104-114),
  `group_items` (117-136).
- These are reusable, side-effect-free transforms over `(feature, color)` lists.
  Two clean options:
  1. Push them into `karyoplot.svg.legend` (or a `karyoplot.mpl/svg`-shared helper)
     alongside the already-present `merge_by_color`/`strip_label_suffixes` — they
     are the same family of legend pre-processing utilities. Then the analysis
     subcommand truly becomes a trivial click->library shim.
  2. Keep them in `karyoscope_analysis/core/legend.py` (analysis-local) if the
     team prefers not to grow karyoplot's surface. The colors-file *format* is a
     KaryoScope-engine artifact, so a thin `core/io/colors.py` reader for
     `load_colors_file` is also defensible (matches the gold-standard `core/io/`
     convention).
- Suggested split for the new layout:
  - `commands/draw_legend.py` — click command: option definitions + orchestration
    (mirrors `karyoscope/commands/annotate.py` structure: thin, delegates to core).
  - `core/legend.py` (or pushed to karyoplot) — `load_colors_file`, `filter_items`,
    `parse_groups`, `group_items`, merge-override parsing, theme construction.
- The `print(...)` progress lines (142-191) should become `click.echo`/logging to
  match the gold-standard engine convention (annotate.py uses `click.echo` + a
  module logger).

ASSESSMENT: After extracting the four helpers, the command body is genuinely
trivial (load -> transform -> `make_legend_drawing` -> save), so the answer is
"trivial after a small refactor, not as-is."

## 9. Smells / risks / dead code / duplication
- `load_colors_file` (lines 85-91) opens the file without `encoding=` and does not
  guard against a missing/garbled file beyond Python's default exception — no
  user-friendly error (gold-standard uses `click.ClickException`).
- No validation of `--columns`/`--rows` (could be negative/zero), color strings, or
  that `--merge-label` features actually exist in the input.
- 3-tuple / 2-tuple juggling: items are wrapped to 3-tuples in both merge branches
  (lines 158, 161) and then flattened back to 2-tuples for grouping (line 164) and
  re-expanded (line 165). Mildly error-prone; a single consistent shape would be
  cleaner.
- `merge_by_color` returns `(label, color)` 2-tuples (plotlib line 293) but the
  script must re-add the `is_header` flag (line 158) — a contract seam worth a test.
- Interaction gap: when `--merge-same-color` is used, suffix-stripping is skipped
  (else branch line 161 not taken) — `merge_by_color` internally strips suffixes
  only for the chosen shortest label (plotlib lines 290-292), so non-merged labels
  in the merge path may retain suffixes inconsistently. Checkable.
- `print`-based output (not logging) — inconsistent with the engine's
  `click.echo`/logger convention.
- No tests, no `__all__`, no type hints on most helpers (only `load_colors_file`
  is annotated, line 82).
- DEAD/UNUSED IMPORT RISK: none found — `Theme`, `make_legend_drawing`,
  `merge_by_color`, `strip_label_suffixes` are all used.

## 10. Testability notes
- The four helpers are pure functions over in-memory lists/strings and are directly
  unit-testable: `filter_items` (include/exclude precedence, line 95-100),
  `parse_groups` (malformed groups skipped, line 109), `group_items` (first-match
  assignment, header insertion, trailing unassigned items, lines 122-135).
- `load_colors_file` needs a tiny temp TSV fixture (header skip, `<2` fields skip).
- `main()` is currently hard to test in isolation: it reads argparse + writes a
  file + prints. Splitting orchestration from a pure `build_legend(items, opts) ->
  Drawing` function (no I/O) would allow asserting on the returned `Drawing`
  (dimensions, element count) without touching disk — and enables a golden-SVG test.
- Integration/golden test: run against a committed sample `.colors.txt` and compare
  the emitted SVG (byte- or structure-level) — consistent with the repo's existing
  "byte-identical bench" practice noted in recent commits.
- No external tools needed; only `karyoplot` must be importable.

## 11. Open questions for the user
1. Should `load_colors_file`/`filter_items`/`parse_groups`/`group_items` move INTO
   `karyoplot.svg.legend` (next to `merge_by_color`/`strip_label_suffixes`), or stay
   analysis-local in `karyoscope_analysis/core/`? (Determines how thin the
   subcommand becomes.)
2. Is the `.colors.txt` format owned by the core `karyoscope` engine? If so, should
   the reader live in `core/io/` and be shared rather than redefined here?
3. Should the legacy argparse defaults (dark theme `#000000`/white text) be the
   command defaults, or should the command adopt karyoplot's `DEFAULT_THEME`
   (which is also dark — plotlib line 331) and drop the inline `Theme` construction?
4. Should progress output move from `print` to `click.echo`/logging to match the
   engine convention?
5. Should the suffix-stripping inconsistency between the merge and non-merge paths
   (section 9) be normalized?
6. Are `--include`/`--exclude` intended to be exact-match only (current behavior,
   lines 97/100), while `--groups` is prefix-match (line 127)? The asymmetry may
   surprise users.
