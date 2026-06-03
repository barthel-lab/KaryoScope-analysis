# Audit: KaryoScope_select_representatives.py (519 lines)

## 1. Purpose
Selects the best representative reads for each cluster so they can be passed to
`KaryoScope_cluster_plot.py` via `--reads-file`. Per the module docstring (lines
3-9), selection balances (1) feature matching — reads must contain the cluster's
defining features — and (2) read length, preferring longer reads for
visualization. Three selection modes coexist (see §5): length-only,
cluster `_top`-feature scoring, and named feature-group scoring.

## 2. CLI surface
Argparse, defined in `parse_args()` (lines 55-71). FACT — flags:
- `--cluster-analysis` (required) — path to `cluster_analysis.tsv`. NOTE: this
  arg is declared but never read anywhere in the program (see §9).
- `--read-assignments` (required) — `sequence_assignments.tsv`.
- `--cluster-labels` (default None) — labels file with `_top` columns (TSV or Excel).
- `--feature-groups` (nargs='+', default None) — named groups e.g.
  `telomeric=canonical_telomere,noncanonical_telomere`.
- `--bed-prefix` (required) — base dir for BED files.
- `--database` (default `KS_human_CHM13`), `--smoothness` (default `smoothed`).
- `--n-per-cluster` (int, default 5), `--feature-min-pct` (float, default 10.0).
- `--clusters` (comma-separated IDs; default all).
- `--preferred-min-length` (int, default 20000), `--preferred-max-length` (int, default 30000).
- `--output` (required) — output TSV.

No click; no subcommands; `main()` is the entry point (lines 357-515, run at 518-519).

## 3. Inputs & outputs
Inputs (FACT):
- `--read-assignments` TSV read with `pandas.read_csv(sep='\t')` (line 374);
  column `read` renamed to `sequence` if needed (lines 375-376); requires
  `cluster`, `sample`, `sequence`, `centroid_distance` columns and optionally
  `read_span`/`read_length`.
- `--cluster-labels` TSV or Excel (`.xlsx`/`.xls` → `read_excel`, else `read_csv`),
  lines 382-385; expected `cluster_id` column (line 434) and `*_top` columns
  (`region_top`, `subtelomeric_top`, `repeat_top`; lines 77, 390).
- BED feature files discovered by glob-like fixed path patterns (lines 107-110),
  supporting `.bed.gz` and `.bed`. Featuresets hardcoded to
  `['region', 'subtelomeric', 'repeat']` (line 412).

Outputs (FACT):
- `--output` TSV via `to_csv(sep='\t')` (line 482); columns vary by mode (lines
  462-479).
- A sibling `*.reads.txt` file (path derived via `.replace('.tsv', '.reads.txt')`,
  line 486) containing one sequence ID per line (lines 487-490).
- Progress/summary printed to stdout (lines 360, 493-515).

## 4. Pipeline / control flow (key functions + line numbers)
- `main()` (357): parse args → optionally parse feature groups (`parse_feature_groups`,
  345) → load read assignments (374) → optionally load labels (381-394) →
  determine cluster IDs (396-401) → optionally pre-load BED features
  (`load_all_features`, 91; gated on labels/groups, 408-420) → per-cluster loop
  calling `select_representatives` (426-445) → `normalize_by_rank` (453) → build
  rows + write outputs (456-490) → summary (492-515).
- `parse_top_features()` (74): regex-parses `name (NN.N%)` tokens from `_top`
  columns, keeping features ≥ `min_pct`. Returns `(feat, pct, featureset)` tuples.
- `load_all_features()` (91): reads BED files into a nested cache
  `{(read_id, sample): {featureset: {feature: bp}}}`, summing `end - start` per feature.
- `check_read_features_cached()` (133): returns `(has_top_feature, n_matches)`.
- `score_multi_feature_groups()` (150): per-group total bp + count of groups present.
- `select_representatives()` (173): tier-based selection (see §5).
- `normalize_by_rank()` (298): reorders reads so the read at rank N has similar
  length across clusters, using per-rank median target lengths.
- Helpers `_get_length` (35), `_create_rep_dict` (42).

## 5. Key design decisions (cite lines)
- Length tiering in `select_representatives` (lines 200-218): candidates split
  into `in_range` (within preferred min/max, sorted by distance to midpoint),
  `below_range` (≥ `MIN_READ_LENGTH`=10000 and < min, sorted by closeness to min),
  `above_range` (> max, sorted by closeness to max), `very_short` (< 10000, sorted
  desc by length). Rationale (docstring 180-182): prefer preferred-range reads
  that carry defining features; for out-of-range reads prefer those nearest the
  boundaries. This is a heuristic, NOT centroid-based selection.
- `centroid_distance` is CARRIED THROUGH from the input file into outputs
  (lines 49, 233, 260, 469) but is NOT used as a ranking/selection criterion
  anywhere — selection is purely length + feature driven. ASSESSMENT: the column
  name suggests centroid-distance selection but the code does not implement it.
- Three mutually-exclusive scoring modes (lines 220-295): (a) length-only when
  neither features nor groups given (221-234); (b) feature-group mode — sort by
  `(-n_groups, -total_bp, dist_to_midpoint)` (237-267); (c) `_top`-feature mode —
  take tiered reads having `has_top or n_matches>0`, then backfill (270-295).
- Rank normalization (298-342): per-rank median target lengths so plots line up
  visually across clusters. Uses `read_span` via `_get_length` default.
- Constants `MIN_READ_LENGTH=10000`, `DEFAULT_PREFERRED_MIN=20000`,
  `DEFAULT_PREFERRED_MAX=30000` (lines 30-32).

## 6. Assumptions (checkable statements)
- BED layout exactly matches the hardcoded pattern
  `{prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed[.gz]`
  (lines 107-110); otherwise features silently absent.
- BED columns: col0=read_id, col1=start(int), col2=end(int), col3=feature; ≥4
  fields required (lines 119-124).
- Read-assignments has `cluster`, `sample`, and `sequence`/`read`,
  `centroid_distance`; `read_span` or `read_length` present (lines 195-196, 374-376, 469).
- Labels file has `cluster_id` and at least one `_top` column when in `_top` mode
  (lines 390-393, 434).
- `--clusters` values are integers (`int(...)`, line 398); `cluster` IDs are ints.
- Output path ends in `.tsv` for the `.reads.txt` sibling to be named sensibly
  (line 486).

## 7. Dependencies
External libs (FACT): `pandas` (line 27); stdlib `argparse`, `gzip`, `os`, `re`,
`collections.defaultdict`. `read_excel` implies an Excel engine (e.g. openpyxl)
is needed when `.xlsx` labels are used — NOT a declared import.
No sklearn, no numpy, no matplotlib.
karyoplot usage: NONE currently. Inter-script: produces output consumed by
`KaryoScope_cluster_plot.py` (`--reads-file`, docstring line 9).
External tools: none invoked directly; depends on upstream telogator/KaryoScope
BED outputs existing on disk.

## 8. Proposed home in new layout
- Subcommand: `karyoscope-analysis select-representatives`.
- Decomposition: thin `commands/select_representatives.py` (click wrapper +
  arg validation) → `core/representatives.py` holding the pure functions
  (`parse_top_features`, `score_multi_feature_groups`, `check_read_features`,
  `select_representatives`, `normalize_by_rank`, `parse_feature_groups`,
  `_get_length`).
- karyoplot / io push-down candidates:
  - BED feature loading (`load_all_features`, 91) and the hardcoded telogator BED
    path pattern (107-110) belong in shared IO. Candidate: `karyoplot.core.io`
    (or a new `karyoscope_analysis/core/io/features.py`). NOTE the gold-standard
    engine already has `karyoscope/core/io/features.py` and `telo.py` — the path
    convention likely should live next to those readers, not be re-hardcoded here.
  - `_get_length` read-length/read_span fallback duplicates logic likely present
    in other analysis scripts — candidate for a shared `core` helper.
- This script does NOT overlap with `karyoplot.mpl.comparison` (that module is
  feature-rate condition comparison, unrelated to read selection).

## 9. Smells / risks / dead code / duplication (line-cited)
- BUG (argument-order) at lines 47-48: `_get_length(row, 'read_length', length_column)`
  and `_get_length(row, 'read_span', length_column)`. `_get_length(row,
  primary_col, fallback_col, default=0)` (35) takes the SECOND positional as
  `primary_col`, so here `read_length`/`read_span` is treated as primary and
  `length_column` as the fallback — workable when `length_column` is the other
  real column, but the intent (primary=length_column with fixed fallback) is
  inverted vs. the simple-mode usage at 230-231 / 258-259 which pass the same.
  ASSESSMENT: confusing and fragile; the helper's signature does not match how
  it is called.
- Unused required arg `--cluster-analysis` (line 57) — declared, never referenced.
  Misleading required input.
- Silent failure: `load_all_features` wraps file parsing in bare
  `except Exception: pass` (lines 125-126) — malformed/missing BEDs vanish with
  no warning; combined with the silent path-pattern assumption this can yield
  zero features and a confusing "length-only-like" result.
- `normalize_by_rank` lambda closure over `target` (line 337) is fine here (called
  immediately) but `min(..., key=lambda i: ...available[i]...)` mutates
  `available` via `.pop` inside the loop — correct but easy to break.
- `_create_rep_dict` emits `has_top_feature`/`n_feature_matches` even though
  the length-only and feature-group paths build their own dicts — three parallel
  dict shapes (lines 227-233, 255-265, 44-52) that must stay in sync with the
  output-row builder (462-479). Duplication risk.
- `score_multi_feature_groups` iterates all featuresets × all group features for
  every read (lines 162-169) — O(reads × groups × feats); fine at current scale,
  noted for large inputs.
- Two unused module constants are shadowed: `DEFAULT_PREFERRED_MIN/MAX` (31-32)
  are defaults for `select_representatives` but `main` always passes
  `args.preferred_*`, so the constant defaults are effectively dead in practice.

## 10. Testability notes
Prime pure-function unit targets (no I/O):
- `parse_top_features` (74) — regex + threshold; easy table-driven tests.
- `parse_feature_groups` (345) — including the malformed (`no '='`) warning path.
- `score_multi_feature_groups` (150) — deterministic bp sums.
- `check_read_features_cached` (133) — has_top / n_matches given a cache dict.
- `select_representatives` (173) — feed a small DataFrame + fake cache; assert
  tier ordering and mode selection (the three branches).
- `normalize_by_rank` (298) — assert per-rank median targeting with crafted lengths.
`load_all_features` (91) needs a tmp BED fixture (integration-level).

## 11. Open questions for the user
1. Is `centroid_distance` meant to influence selection (the filename/columns imply
   centroid-distance selection) or is it intentionally only metadata for output?
2. Should `--cluster-analysis` be removed (currently unused/required) or wired in
   (e.g., to source per-cluster centroid info)?
3. Are the hardcoded featuresets `['region','subtelomeric','repeat']` (412) and the
   telogator BED path convention (107-110) stable, or should they be configurable /
   sourced from the shared engine's IO layer?
4. Should silent BED-parse failures (125-126) become warnings/errors?
5. Is the `read_length` vs `read_span` distinction settled — which is canonical
   for "length" in the new package?
