# Audit: KaryoScope_compare_clusterings.py (577 lines)

## 1. Purpose
Compares two clustering results produced from different featuresets (e.g. region
vs repeat) on the same set of reads, to assess concordance and identify
complementary signal (docstring lines 3-6). Emits a text report, a multi-page
PDF of concordance plots, a cluster-to-cluster overlap matrix, and an auto-label
cross-tabulation.

## 2. CLI surface
Argparse, in `main()` (lines 362-381). FACT — flags:
- `--clustering1` (required) — first `sequence_assignments.tsv`.
- `--clustering2` (required) — second `sequence_assignments.tsv`.
- `--label1` (default `clustering1`), `--label2` (default `clustering2`) — display labels.
- `--output-prefix` (required) — prefix for the four output files.
- `--dark-mode` (store_true, default False) — dark plot styling.

No click; no subcommands. Entry point `main()` (run at lines 576-577).

## 3. Inputs & outputs
Inputs (FACT):
- Two assignment TSVs via `read_csv(sep='\t')` (line 41); `sequence` renamed to
  `read` if needed (44-45). Requires `read`/`sequence` and `cluster` columns.
- Auto-discovered sidecar files via string `.replace` on the assignments path
  (lines 49, 59): `*.cluster_analysis.tsv` (for the `enrichment` map; columns
  `cluster_id`, `enrichment`, line 54) and `*.cluster_annotations.tsv` (for the
  `cluster_name` map; columns `cluster_id`, `cluster_name`, lines 62-64). Both
  optional — absence → `'unknown'` / `'(unlabeled)'` fills (56, 66-67).

Outputs (FACT, all prefixed by `--output-prefix`):
- `.comparison_report.txt` — `generate_report` (235-359).
- `.comparison_plots.pdf` — up to 5 pages via `PdfPages` (443-521).
- `.comparison_matrix.tsv` — cluster×cluster crosstab (525-529).
- `.comparison_labels.tsv` — auto-label crosstab with pct_of_a/pct_of_b
  (536-549), only when both label columns present.
- Stdout summary (552-573).

## 4. Pipeline / control flow (key functions + line numbers)
- `load_clustering()` (39): loads assignments, builds enrichment & cluster_name
  maps; returns `(df, enrichment_map, cluster_name_map)`.
- `main()` (362): parse args → load both clusterings (388-398) → build label-suffixed
  subframes and inner-merge on `read` (401-408; exits if empty merge, 410-412) →
  `generate_report` (417) → PDF pages (443-521) → matrix (525-529) → label crosstab
  (535-549) → summary (552-573).
- Metric/compute helpers: `compute_cluster_overlap_matrix` (72),
  `compute_enrichment_concordance` (77), `compute_adjusted_rand_index` (82),
  `compute_normalized_mutual_info` (88), `compute_label_retention` (94),
  `compute_label_purity` (115).
- Plot helpers: `plot_enrichment_sankey` (136, actually a heatmap),
  `plot_cluster_size_comparison` (155), `plot_top_cluster_mappings` (171),
  `plot_concordance_by_cluster` (205).
- Report writer: `generate_report` (235).

## 5. Key design decisions (cite lines; metrics + WHY)
- Reads are matched across the two clusterings by an INNER merge on `read`
  (line 407); only reads present in both are compared (exit if none, 410-412).
- Concordance metrics in the report (258-281):
  - Raw enrichment-label agreement = fraction of reads whose `*_enrich` labels are
    equal (261-262). Drives the qualitative LOW/MODERATE/HIGH interpretation at
    <0.5 / <0.75 / else (343-356).
  - Adjusted Rand Index (ARI) via `compute_adjusted_rand_index` (82) and
    Normalized Mutual Information (NMI) via `compute_normalized_mutual_info` (88).
    WHY (stated only implicitly via docstring "assess concordance"): ARI/NMI are
    standard label-invariant cluster-agreement measures — appropriate because
    cluster IDs are arbitrary between the two runs, so a permutation-invariant
    metric is required. No explicit rationale comment in code.
- Label flow analysis: `compute_label_retention` (94, forward: % of each A-label
  retaining same label in B) and `compute_label_purity` (115, reverse). These
  operate on auto-label STRINGS, treating a match as label equality (lines 100,
  121) — only meaningful when both runs share a label vocabulary.
- Plots favor enrichment-flow heatmaps (row-normalized %, 141, 498) and
  Pre/Post color coding in `plot_top_cluster_mappings` (189-196) — domain-specific
  (pre/post BIR enrichment categories).
- `pdf.fonttype`/`svg.fonttype` set to keep text editable (34-36) and re-applied
  after `style.use` resets them (439-441) — deliberate.

## 6. Assumptions (checkable statements)
- Both assignment files share a common `read` ID space (merge on `read`, 407).
- `cluster_analysis.tsv` has `cluster_id` + `enrichment`; `cluster_annotations.tsv`
  has `cluster_id` + `cluster_name` (54, 62-64).
- Sidecar files are co-located and named by the exact `.replace` substring
  `'.sequence_assignments.tsv'` (49, 59) — fragile to path variations.
- Enrichment category strings contain `'Pre'`/`'Post'` substrings for the
  color logic to fire (189-196).
- sklearn is importable for ARI/NMI — but see §9 (the ARI import name is wrong).

## 7. Dependencies
External libs (FACT): `matplotlib` (27, `PdfPages` 31), `numpy` (28, imported but
NOT used — see §9), `pandas` (29), `seaborn` (30), `scipy.stats` (32, imported
but NOT used). `sklearn.metrics` imported lazily inside the two metric helpers
(84, 90). stdlib `argparse`, `os`, `sys`.
karyoplot usage: NONE.
Overlap with `karyoplot.mpl.comparison`: NO functional overlap. That module
(volcano/dot-strip/lollipop, log2FC, Fisher/Mann-Whitney via `mpl.statistics`)
compares per-sample FEATURE RATES between biological CONDITIONS. This script
compares two CLUSTERINGS of the same reads (ARI/NMI, crosstabs, label flow). The
name collision is superficial. Shared style primitives DO overlap (dark-mode,
fonttype, save) — see §8.
Inter-script: consumes `*.sequence_assignments.tsv` / `*.cluster_analysis.tsv` /
`*.cluster_annotations.tsv` produced by `KaryoScope_cluster_analysis.py` /
`KaryoScope_cluster_annotate.py`. External tools: none.

## 8. Proposed home in new layout
- Subcommand: `karyoscope-analysis compare-clusterings`.
- Decomposition: thin `commands/compare_clusterings.py` (click) → pure metric
  layer `core/clustering_comparison.py` (the `compute_*` functions) →
  report writer `core/clustering_comparison_report.py` (or keep in core) →
  plotting `core/clustering_comparison_plots.py` (or push to karyoplot, below).
- karyoplot push-down candidates:
  - Style/setup boilerplate (dark-mode palette, `pdf.fonttype`/`svg.fonttype`,
    figure facecolor) at 34-36, 426-441 duplicates `karyoplot.mpl.style`
    (`apply_default_style`, `fg_color`, `save_fig`) — replace with those.
  - Generic plotting helpers (`plot_enrichment_sankey` heatmap 136,
    `plot_cluster_size_comparison` 155, stacked-bar concordance 205) are
    reusable; a NEW `karyoplot.mpl` module — e.g. `karyoplot.mpl.clustering` or
    `clustering_concordance` — is the right home (do NOT fold into
    `mpl.comparison`, which has different semantics; see §7).
  - The ARI/NMI metric functions are pure and belong in
    `karyoplot.mpl.statistics` (or a `core` stats module) alongside the existing
    Fisher/Mann-Whitney helpers.

## 9. Smells / risks / dead code / duplication (line-cited)
- BUG (silent, high impact): line 84 imports `adjusted_rand_index` from
  `sklearn.metrics`, but the real name is `adjusted_rand_score`. The import
  raises `ImportError`, which is caught at 271-272 and 565-566, so ARI is
  NEVER computed and the report just prints "(sklearn not available)" — even
  when sklearn IS installed. Verified: `adjusted_rand_index` does not exist in
  sklearn; `adjusted_rand_score` does. NMI (line 90) uses the correct name.
- The broad `except ImportError` (271, 280, 565) also masks the wrong-name bug:
  a genuinely-missing-sklearn message hides a typo. Should narrow or surface.
- Dead imports: `numpy as np` (28) and `from scipy import stats` (32) are never
  used.
- Dead/unused params: `generate_report` accepts `df1, df2, name1_map, name2_map`
  (235-236) but never uses `df1`/`df2`/`name1_map`/`name2_map` in the body.
- Dead local: `annot_matrix = ct_pct.copy()` (509) is assigned and never used.
- Misleading name: `plot_enrichment_sankey` (136) draws a heatmap, not a Sankey.
- Sidecar discovery by `str.replace` (49, 59) is fragile — breaks if the literal
  substring is absent or appears elsewhere in the path.
- `load_clustering` signature has `analysis_file=None` param (39) that `main`
  never passes — auto-discovery is the only path used; param is effectively dead
  from the CLI's perspective.
- Title says "row-normalized %, cells <1% hidden" (516) but the mask threshold is
  `< 0.5` (513) — off-by-comment (hides <0.5%, not <1%).
- `compute_label_retention`/`purity` compare label strings for equality across
  the two clusterings (100, 121); if the two runs use disjoint label vocabularies,
  "same_label_pct" is trivially ~0 and misleading. ASSESSMENT: only meaningful
  with a shared vocabulary — undocumented precondition.

## 10. Testability notes
Prime pure-function unit targets (no I/O, deterministic):
- `compute_cluster_overlap_matrix` (72) and `compute_enrichment_concordance` (77)
  — crosstab shape/values on tiny frames.
- `compute_adjusted_rand_index` (82) / `compute_normalized_mutual_info` (88) —
  known-answer cases (identical labels → ARI=1, NMI=1; random → ~0). A test here
  would immediately have caught the `adjusted_rand_index` typo.
- `compute_label_retention` (94) and `compute_label_purity` (115) — forward/reverse
  percentages and top-N formatting on a small merged frame.
- `load_clustering` (39) — needs tmp TSV fixtures (integration-level); good for
  testing the sidecar auto-discovery and the `'unknown'`/`'(unlabeled)'` fills.
Plot functions: smoke-test against a `matplotlib` Agg axis (no assertions on pixels).

## 11. Open questions for the user
1. Is ARI silently missing from existing reports a known issue, or have all prior
   runs been emitting "(sklearn not available)"? (The `adjusted_rand_index` typo
   means ARI has likely never been produced.)
2. Should label retention/purity require a shared label vocabulary between the two
   clusterings, and should the tool warn when vocabularies are disjoint?
3. Is the Pre/Post enrichment color coding (189-196) general, or specific to one
   study? Should it move into a shared theme/palette in karyoplot?
4. Should the comparison plots be reimplemented on `karyoplot.mpl` style helpers
   (and a new `mpl.clustering` module), or stay self-contained in the analysis
   package?
5. Are the sidecar naming conventions (`.cluster_analysis.tsv`,
   `.cluster_annotations.tsv`) guaranteed, or should paths be explicit CLI args?
