# Audit: KaryoScope_cluster_annotate.py (1719 lines)

## 1. Purpose

FACT: Step 2 of a two-step annotation pipeline (docstring lines 8-10). It consumes
per-read sequence annotations (a TSV produced by `KaryoScope_sequence_annotate.py`) plus
clustering outputs and produces a per-cluster annotation TSV. For each cluster it
aggregates pre-computed per-read feature statistics into cluster-level columns
(`{fs}_readpct__{feat}`, `{fs}_bppct__{feat}`, `{fs}_dmax/dmin/dmedian/dfirst/dlast__{feat}`,
plus `dterminal`, `dterminal_min`, `max_block_bp` and `interspersion_*`), optionally
auto-labels clusters via a structural decision tree (`auto_label_cluster`, lines 727-899),
optionally selects representative reads per cluster (lines 1591-1622), and optionally runs a
feature-importance analysis with a 6-panel PDF (lines 1654-1715).

ASSESSMENT: It is an aggregation + heuristic-classification + optional-reporting tool. It does
NOT recompute anything from BED files — all per-read stats are read from the upstream TSV
(see lines 194-233 docstring "Reads pre-computed density columns ... from sequence annotations").

## 2. CLI surface

FACT: `argparse` CLI in `main()` (parser created line 1266). Arguments:

- `--prefix` (required, 1271): auto-finds `{prefix}.read_assignments.tsv` /
  `.sequence_assignments.tsv` (1336-1338) and `{prefix}.cluster_analysis.tsv` (1339); also
  `{prefix}.feature_matrix.npz` for SVD (1677).
- `--sequence-annotations` (required, 1273): per-read annotations TSV.
- `--adaptive-thresholds` (default None, 1275): auto-derived from the sequence-annotations
  path as `{sa_prefix}.adaptive_thresholds.tsv` (1449-1451) if omitted.
- `--featuresets` (default `region,subtelomeric,chromosome,acrocentric,repeat,gene`, 1277).
- `--output/-o` (required, 1279).
- `--top-n` (int, default 3, 1281).
- `--clusters` (CSV of int IDs, default all, 1283).
- `--min-size` (int, default 1, 1285).
- `--exclude-features` (default `*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized,canonical_telomere*`, 1287).
- `--log-file` (BooleanOptionalAction, default True, 1290).
- `--auto-label` (store_true, 1293).
- `--feature-importance` (store_true, 1296).
- `--alt-samples` (CSV sample prefixes, 1299) + `--alt-threshold` (float, default 80, 1301).
- `--select-representatives N` (int, default None, 1303).
- `--rep-strategy` (`annotation`|`centroid`, default `annotation`, 1307).

## 3. Inputs & outputs

FACT — Inputs:
- `{prefix}.read_assignments.tsv` OR `.sequence_assignments.tsv` (1336-1338). Must have
  `cluster`, `sample`, and a read-id column `sequence` (or `read`, renamed at 1352-1353).
  Optionally `read_span`/`read_length` and `centroid_distance` (used by centroid strategy).
- `{prefix}.cluster_analysis.tsv` (1339; warns if missing, 1422). Read for `enrichment`,
  `p_value`, `q_value`, `odds_ratio`, and auto-detected entity stat columns with suffixes
  `_count/_pct/_pval/_odds` (1376-1400).
- `--sequence-annotations` TSV (1432) — primary data source; featuresets/features
  auto-detected from `{fs}_frac__{feat}` columns (load_sequence_annotations, 57-89).
- `{sa_prefix}.adaptive_thresholds.tsv` (optional; columns `featureset`, `feature`,
  `threshold`, 92-106). If missing, thresholds computed from `_frac__` columns (1463-1474).
- `{prefix}.feature_matrix.npz` (optional, for `--feature-importance`; needs `svd_components`,
  `svd_feature_names`, `svd_explained_variance_ratio`, 986-991).

FACT — Outputs:
- `--output` TSV: one row per cluster (1625).
- `{output_stem}.log` when `--log-file` (1322-1327, via `TeeLogger`).
- With `--feature-importance`: `.feature_importance_annotations.tsv` (1669),
  `.feature_importance_svd.tsv` (1690), `.feature_importance_svd_layers.tsv` (1697),
  `.feature_importance.pdf` (1708).

## 4. Pipeline / control flow (key functions + line numbers)

FACT:
- `load_sequence_annotations` (57-89): reads TSV, detects featuresets/features from
  `_frac__` columns.
- `load_adaptive_thresholds` (92-106): TSV → `{fs: {feat: threshold}}`.
- `matches_any_pattern` (109-114): fnmatch wildcard helper for exclude patterns.
- `summarize_featureset` (117-151): top-N feature string by summed `{fs}_bp__` across cluster
  reads; applies exclude patterns (display only).
- `score_cluster_features` (154-169): `{fs}_readpct__` = % reads with `_frac__` > threshold.
- `compute_cluster_bp_scores` (172-191): `{fs}_bppct__` = feature bp / `{fs}_total_bp` ×100.
- `compute_cluster_window_densities` (194-233): median across cluster reads of pre-computed
  density columns (dmax/dmin/dmedian/dfirst/dlast/dterminal/dterminal_min/max_block_bp),
  zero-padded for reads missing from annotations; ×100 except max_block.
- `compute_cluster_interspersion` (236-257): median of `interspersion_*` columns.
- `score_read_against_annotation` (260-388): per-read representativeness score (weighted sum,
  weights adapt for ECTR/satellite labels, 290-308).
- `compute_length_score` (391-400): log2 length-proximity score.
- `normalize_representatives_by_length` (403-451): rank-assign candidates by combined
  feature+length score.
- `select_centroid_representatives` (454-502): pick smallest `centroid_distance`.
- `select_annotation_representatives` (505-647): two-pass annotation-aware selection.
- `TeeLogger` (650-673): stdout+file tee.
- `_print_params_and_command` (680-724): param table.
- `auto_label_cluster` (727-899): decision-tree labeler (see §5).
- `analyze_annotation_importance` (904-968): per-column CV/sparsity + Spearman corr.
- `analyze_svd_loadings` (971-1057): per-feature SVD importance + layer decomposition.
- `plot_feature_importance` (1060-1262): matplotlib 6-/4-panel PDF.
- `main` (1265-1715): orchestration — load inputs, per-cluster loop (1488-1583), build/sort
  DataFrame (1586-1589), optional rep selection (1592-1622), save (1625), summaries, optional
  feature-importance block (1655-1715).

## 5. Key design decisions (cite lines)

FACT — Auto-labeling (`auto_label_cluster`, 727-899): structural decision tree keyed on
terminal telomere density. Thresholds are hard-coded constants (746-756): `CAN_ECTR=70`,
`NCAN_ECTR=10`, `CAN_SUB=15`, `NCAN_SUB=5`, `DMAX_HIGH=35`, `ENRICH_DMAX=35`,
`NCAN_VARIANT=25`, `SAT_DOMINANT=80`, `CT_ENRICH=20`, `ALT_BLOCK_BP=6000`, `ARM_PRESENT=30`.
Stated WHY (docstring 728-734): dfirst/dlast classify telomere-at-both-ends (ECTR),
one-end (Subtelomere), or internal (Interstitial); Subtelomeres with a ≥6 kb contiguous
canonical block become "Type II ALT subtelomere"; enrichment qualifiers appended. Tree order
(859-899): (1) ECTR if both ends telomeric OR one end without arm/SegDup (862); (2) Subtelomere
/ Type II ALT (867-870); (3) Interstitial telomere (873); (4) Interstitial ITS/TAR1 (877);
(5) satellite-dominant by readpct ≥80 (895); (6) unlabeled "" (899).

FACT — Enrichment uses `dmax` for qualifiers but `readpct` for satellite *dominance* (rule 5),
commented at 798. v2 satellite types explicitly enumerated for forward-compat (799-821).

FACT — Type I ALT relabeling (1573-1581): if sum of `{sample}_pct` over `--alt-samples`
exceeds `--alt-threshold`, prepend "Type I ALT " (lowercasing first char), excluding ECTR /
Type II ALT labels. Sample-specific business logic baked into a generic tool.

FACT — Representative selection has two strategies (1592-1606): `centroid` uses
`centroid_distance` from assignments (454-502); `annotation` (default) scores reads vs the
cluster annotation profile, with a two-pass length normalization that re-targets to the
longest Pass-1 representative (619-647). WHY (docstring 403-419, 454-465): normalize lengths
across clusters; centroid is embedding-space, heuristic-free.

FACT — Featureset prefix selection prefers v2 `region_subtelomere_flat`, falls back to legacy
`telomere_region`, for both auto-label (1565-1568) and reps (1600-1604).

## 6. Assumptions (checkable statements)

FACT/ASSESSMENT:
- Featuresets/features are derivable from `{fs}_frac__{feat}` column names (71-79). If the
  upstream TSV omits `_frac__` columns, NO featuresets are detected.
- `{fs}_bp__{feat}` and `{fs}_total_bp` columns exist for bp scoring (125, 178-181). If
  `{fs}_total_bp` missing, bppct silently returns `{}` → all `_bppct__` become 0 (181-183).
- Per-read density columns are 0-1 raw fractions; cluster output scales ×100 (docstring 199;
  226-228). `max_block_bp` is NOT scaled (225-226).
- Read-id column is named `sequence` (or `read`) in BOTH assignments and annotations (1352,
  119, 156). Annotations TSV is assumed to use `sequence` (not renamed there).
- Auto-label requires `region_subtelomere_flat` OR `telomere_region` in `--featuresets`, else
  hard exit (1567-1569).
- rDNA is read from a DIFFERENT featureset prefix: hard-coded `acrocentric_dmax__rDNA`
  (785) — not `{pfx}`.
- Feature-vocab assumption: satellite names follow v2 canonical with v1 fallback via
  `lookup_satellite_col` (see §7); `_feature_vocab.py` documents the v1/v2 alias semantics.
- Color files `KS_human_CHM13.chromosome_acrocentric*.colors.txt` exist in `scripts/` but are
  NOT referenced by this script (no color loading here).
- Cluster IDs are integers (`--clusters` parses with `int`, 1479).

## 7. Dependencies

FACT — External libs: `argparse`, `fnmatch`, `math`, `os`, `sys` (41-45); `pandas` (50);
lazily imported inside functions: `numpy`, `scipy.stats` (910-911), `numpy` (981, 1065, 1656),
`matplotlib`+`PdfPages` (1066-1069), `scipy.cluster.hierarchy`/`scipy.spatial.distance`
(1154-1155), `matplotlib.patches.Patch` (1123).

FACT — karyoplot usage: `from karyoplot.core.fonts import register_fonts, resolve_family`
(1072), called at 1073-1078 to register brand fonts and set `font.family`. This is the ONLY
karyoplot import; confirmed `karyoplot/core/fonts.py` exposes `register_fonts` and
`resolve_family` (plus `set_default_font`, `is_available`, `pil_font`).

FACT — Inter-script / shared deps: `_feature_vocab.py` (imported via `sys.path.insert` at
53-54: `from _feature_vocab import lookup_satellite_col`). Uses `lookup_satellite_col` at
368, 822, 827. `_feature_vocab.py` (109 lines) defines v1/v2 satellite vocab sets, alias
tables (`SATELLITE_V1_TO_V2`, `SATELLITE_V2_TO_V1`), arm/ct/telomere feature sets, and helpers
`is_satellite` / `lookup_satellite_col`. CONFIRMED `_feature_vocab` is shared: also imported by
`KaryoScope_sequence_annotate.py` (the upstream step-1 script). The Barthel color palette and
featureset→color map are hard-coded locally (1082-1096), NOT from karyoplot.core.colors.

FACT — External tools: none invoked (no subprocess). Upstream dependency on
`KaryoScope_sequence_annotate.py` and `KaryoScope_cluster_analysis.py` outputs only.

## 8. Proposed home in new layout

ASSESSMENT — Subcommand name: `karyoscope-analysis annotate-clusters` (or `cluster-annotate`,
matching the file stem; pairs with a step-1 `annotate-sequences`).

Decomposition:
- `commands/cluster_annotate.py` — thin click subcommand: option definitions (mirroring §2),
  logging setup, orchestration calling into core. Replaces `main()` (1265-1715).
- `core/cluster_annotate.py` — pure aggregation/labeling logic:
  `summarize_featureset`, `score_cluster_features`, `compute_cluster_bp_scores`,
  `compute_cluster_window_densities`, `compute_cluster_interspersion`, `auto_label_cluster`,
  the Type I ALT relabel helper, and the per-cluster row builder (1488-1583).
- `core/representatives.py` — `score_read_against_annotation`, `compute_length_score`,
  `normalize_representatives_by_length`, `select_centroid_representatives`,
  `select_annotation_representatives`. NOTE this overlaps with standalone
  `KaryoScope_select_representatives.py` (which has its own `select_representatives`,
  `normalize_by_rank`, `_get_length`) — consolidate both into one core module to remove
  duplication.
- `core/feature_importance.py` — `analyze_annotation_importance`, `analyze_svd_loadings`
  (numeric analysis, no plotting).
- `core/io/` — readers for assignments / cluster_analysis / sequence_annotations /
  adaptive_thresholds / feature_matrix.npz (currently inline `pd.read_csv` at 1349, 1368,
  63, 98, and `np.load` at 983).

karyoplot push-down candidates + module:
1. `TeeLogger` (650-673) → `karyoplot.core.io` (or analysis-pkg `core/io/logging.py`) — generic
   stdout/file tee, duplicated across analysis scripts.
2. `plot_feature_importance` (1060-1262) → `karyoplot.mpl` (new `feature_importance.py`, or
   reuse `karyoplot.mpl.heatmap` for panels C/D and `karyoplot.mpl.style` for the palette).
   The hard-coded Barthel `COLORS`/`FEATURESET_COLORS` (1082-1096) belong in
   `karyoplot.core.colors` / `karyoplot.mpl.style`.
3. The font/rcParams setup (1072-1079) is already a karyoplot delegate; a small
   `karyoplot.mpl.style.apply_brand_rcparams()` would remove the inline `plt.rcParams.update`.
   `analyze_annotation_importance`'s CV/sparsity/Spearman block is a candidate for
   `karyoplot.mpl.statistics`.

`_feature_vocab` placement: it is shared across analysis scripts (sequence_annotate +
cluster_annotate) and encodes domain feature vocabulary, NOT plotting. ASSESSMENT: put it in
the analysis package as `karyoscope_analysis/core/feature_vocab.py` (importable, no sys.path
hack). If the same v1/v2 vocab is also needed by the core `karyoscope` engine or plotlib,
promote to `karyoplot.core` (e.g. `karyoplot.core.feature_vocab`) — but only if a cross-repo
consumer actually exists; default to analysis-pkg core.

## 9. Smells / risks / dead code / duplication (line-cited)

FACT/ASSESSMENT:
- `sys.path.insert(0, ...)` import hack (53-54) — fragile; resolved by packaging.
- Module-level side effect: `_original_command = ' '.join(sys.argv)` at 48, before pandas
  import (50) — works but unusual ordering.
- `score_read_against_annotation` `s_bppct` (349-358) approximates read bppct from median
  density (comment 355 "Approximate") — heuristic, possibly inaccurate.
- Hard-coded `acrocentric_dmax__rDNA` (785) couples auto-label to a specific featureset name
  while everything else uses `{pfx}` — silent 0 if that featureset is absent.
- Type I ALT relabel (1573-1581) is sample-name-specific business logic embedded in a generic
  tool; label-mutation via string surgery (`label[0].lower() + label[1:]`, 1581) is brittle.
- `summarize_featureset` exclude filter (141-142) is applied only to the `{fs}_top` display
  string, NOT to the numeric `_readpct/_bppct/_dmax...` columns — easy to misread as global.
- `score_cluster_features` (162-166): comment says zero-fill for missing reads, but
  `cluster_ann[col].fillna(0)` only fills reads PRESENT in annotations; reads absent from
  annotations are correctly handled by dividing by `n_reads` (the full cluster size). The
  `.fillna(0)` and the divisor interact subtly — verify intent.
- `_print_params_and_command` `_fmt` (684-694) default-detection logic is convoluted
  (None-vs-default branching) and only some params routed through it (697-699 bypass it).
- `compute_length_score` median selection `sorted(...)[len//2]` (427, 644) is a non-interpolated
  median (off-by-one for even n) vs `pd.Series.median()` used elsewhere — inconsistent.
- Lazy in-function imports of numpy/scipy/matplotlib (910, 981, 1065-1069) scattered — fine for
  optional deps but should be consolidated.
- No `dead code` found per se; `entity_columns` summary loop (1647-1652) recomputes `seen`
  already available as `seen_entities`.
- Duplicated read-span column detection (`read_span`/`read_length`) appears 3× (491, 515-525,
  555) — extract a helper.

## 10. Testability notes

FACT/ASSESSMENT:
- Pure/easily unit-testable (dict/Series in, dict/str out, no I/O): `matches_any_pattern`
  (109), `score_cluster_features` (154), `compute_cluster_bp_scores` (172),
  `compute_cluster_window_densities` (194), `compute_cluster_interspersion` (236),
  `score_read_against_annotation` (260), `compute_length_score` (391),
  `normalize_representatives_by_length` (403), `auto_label_cluster` (727),
  `analyze_annotation_importance` (904), `analyze_svd_loadings` (971). `auto_label_cluster` is
  the highest-value target — its decision tree (862-899) deserves table-driven tests per branch
  (ECTR / Subtelomere / Type II ALT / Interstitial telomere / ITS-TAR1 / satellite / unlabeled).
- Fixture-needing (read DataFrames/TSVs): `load_sequence_annotations`,
  `load_adaptive_thresholds`, `select_centroid_representatives`,
  `select_annotation_representatives` — small synthetic TSVs suffice.
- Integration-only: `plot_feature_importance` (1060, PDF output — smoke test it runs without
  error), and `main` end-to-end (golden TSV from a fixed input set).
- `TeeLogger` (650) reassigns global `sys.stdout` (1327) and never restores it — complicates
  test isolation; tests should patch/restore stdout.

## 11. Open questions for the user

1. Should `_feature_vocab` live in the analysis package's `core/` (default proposal) or be
   promoted into `karyoplot.core` for cross-repo reuse? Is any non-analysis repo a consumer?
2. The auto-label thresholds (746-756) are hard-coded magic numbers tuned for human CHM13.
   Should these be configurable (CLI/config file) or remain fixed in code?
3. Type I ALT relabeling (1573-1581) hard-codes ALT sample logic. Keep it in the generic
   subcommand, or move to a separate post-processing step / config?
4. `KaryoScope_select_representatives.py` is a standalone script with overlapping
   representative-selection logic. Should the two be unified into one `core/representatives.py`
   (and which strategy set is authoritative)?
5. The Barthel palette + `plot_feature_importance` are local. Push the palette to
   `karyoplot.core.colors` and the plotting to `karyoplot.mpl`, or keep feature-importance
   plotting analysis-local since it is bespoke?
6. Confirm the canonical featureset prefix going forward: `region_subtelomere_flat` (v2) vs
   `telomere_region` (legacy) — can legacy support be dropped?
