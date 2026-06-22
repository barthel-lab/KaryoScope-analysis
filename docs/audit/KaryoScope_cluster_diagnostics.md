# Audit: KaryoScope_cluster_diagnostics.py (1731 lines)

## 1. Purpose
FACT: A post-hoc, read-level diagnostic/visualization script. It consumes an
already-clustered, *annotated reads* TSV (one row per read, with per-read
alignment metrics + a `cluster` column) and produces (a) exploratory diagnostic
PDFs comparing per-read metric distributions across clusters / samples /
enrichment categories, and (b) "publication-quality" figures + a statistics
TSV comparing user-selected target clusters vs all others. Module docstring
(lines 2-40) states it "Creates diagnostic plots comparing clusters across
various metrics from annotated sequence data" and supports "both exploratory
analysis and publication-quality figure generation."

ASSESSMENT: This is a *downstream consumer* of clustering output, NOT a
clustering or k-selection tool. Despite "diagnostics" in the name, it does NOT
overlap with cluster_analysis's k-selection diagnostics (silhouette/Davies-
Bouldin/knee-elbow) — see section 8/9. Its real domain is per-cluster metric
comparison + significance testing.

## 2. CLI surface
FACT: argparse (`main`, line 1224). Single entry point; no subcommands.
Arguments:
- `--annotated` (required, line 1230): annotated reads TSV.
- `--cluster-analysis` (line 1232): cluster_analysis TSV for enrichment/q-value.
- `--output-prefix` (required, line 1234).
- `--dark-mode` (line 1236).
- `--format` {pdf,svg,png} (line 1238): default pdf; only applied to pub figures.
- `--compare-clusters "82,88"` (line 1243): box comparison of listed clusters vs others.
- `--compare-significant` (line 1245): each FDR-sig cluster vs pooled non-sig.
- Publication group (lines 1249-1271): `--pub-comparison`, `--pub-counts`,
  `--pub-heatmap`, `--pub-stats`, `--pub-all`, `--clusters`, `--cluster-labels`
  ("82:α-sat,88:rDNA"), `--heatmap-order` {enrichment,size,hierarchical},
  `--min-cluster-size` (default 5), `--no-exploratory`.

ASSESSMENT: `--no-exploratory` (line 1270) is declared but **never read**
anywhere in `main` — exploratory plots always run. Dead flag (section 9).

## 3. Inputs & outputs
FACT inputs:
- `--annotated` TSV (`load_annotated_data`, line 187): must contain `cluster`
  (line 193). Optional columns consumed: `sample`, `group`, `enrichment`,
  `sequencing_approach`, and numeric metrics `read_length`, `centroid_distance`,
  `primary_mapq`, `primary_de`, `primary_align_len`, `primary_align_fraction`,
  `total_align_fraction`, `n_alignments`, `n_secondary`, `n_supplementary`,
  plus legacy `mapq/de/align_len/align_fraction` (line 378).
- cluster_analysis TSV: read at lines 1311-1312 and AGAIN at 1547-1550;
  requires `cluster_id` + `enrichment` (+ optional `q_value`).
- Auto-discovery (lines 1316-1327): if `--cluster-analysis` not given, replaces
  `.read_assignments.annotated.tsv` -> `.cluster_analysis.tsv` in the annotated
  path and loads if present.

FACT outputs (prefix-based): `{prefix}.cluster_metrics.pdf`,
`.cluster_composition.pdf`, `.cluster_summary.tsv` (always);
`.cluster_comparison.pdf` (`--compare-clusters`);
`.significant_clusters_comparison.pdf` (`--compare-significant`);
`.pub_comparison.{fmt}`, `.pub_counts.{fmt}`, `.pub_heatmap.{fmt}`,
`.pub_statistics.tsv` (pub flags). No SVG/PNG twin output; no return codes.

## 4. Pipeline / control flow (key functions + line numbers)
FACT:
- Stat/format helpers: `format_pvalue` (101), `format_pvalue_stars` (113),
  `compute_effect_size` (127, rank-biserial from Mann-Whitney U),
  `format_metric_label` (155), `format_metric_title` (171).
- IO: `load_annotated_data` (187).
- Exploratory per-axis plotters (seaborn box+strip + Kruskal-Wallis title):
  `plot_metric_by_cluster` (197), `plot_metric_by_enrichment` (228),
  `plot_metric_by_sample` (254), `plot_metric_by_cluster_and_sample` (285),
  `plot_cluster_composition` (315), `plot_cluster_sizes` (339),
  `plot_correlation_heatmap` (363).
- Summary: `compute_cluster_summary` (376) -> per-cluster mean/median/std table.
- Comparison figures: `plot_significant_vs_nonsig` (412, box for continuous +
  stacked-bar proportions for `n_*` count metrics), `plot_cluster_vs_others`
  (570, Mann-Whitney + effect size) **[never called]**.
- Publication: `pub_comparison_figure` (684, boxplots + significance brackets,
  Mann-Whitney + rank-biserial), `pub_count_distribution_figure` (872, stacked
  bars + Chi-square), `pub_summary_heatmap` (1046, z-scored median matrix,
  RdBu_r imshow, enrichment-colored row labels, q-value stars, orderings
  enrichment/size/hierarchical-ward), `pub_statistics_table` (1172).
- `main` (1224): parse -> set style (1276-1299) -> load (1306) -> dual enrichment
  load (1309-1327) -> exploratory PDFs (1343-1523) -> summary TSV (1527) ->
  re-load enrichment/q (1544-1550) -> compare flags (1553-1616) -> pub block
  (1622-1716, each wrapped in broad `try/except Exception`) -> summary print.

## 5. Key design decisions (cite lines; diagnostic metrics + WHY if stated)
FACT (tests, with stated rationale where present):
- Kruskal-Wallis for multi-group exploratory metric panels (lines 216, 242, 273,
  302). Non-parametric; chosen for >=2 groups. No WHY stated.
- Mann-Whitney U + rank-biserial effect size for cluster-vs-others 2-group
  comparisons (`compute_effect_size` 127-152; pub_comparison 820-821; stats
  table 1195-1196; plot_cluster_vs_others 657-662). Rationale implied by
  "rank-biserial correlation (effect size for Mann-Whitney U)" docstring (129).
- Chi-square contingency for discrete count metrics `n_*` (line 1011); counts
  binned with a `max_val+ "N+"` overflow bucket capped at 10 (lines 471-472,
  929-932); only tested when each side has >=5 obs and >=2 nonzero columns
  (1004-1009) — WHY: avoid invalid chi-square on sparse tables.
- Significance stars 4-tier incl. `****` p<0.0001 (`format_pvalue_stars` 113-124).
- Heatmap uses MEDIAN per cluster then z-scores across clusters (1076-1083),
  clamps colorbar to vmin/vmax -3/3 (1109); constant columns std=0 -> 1.0 to
  avoid NaN (1081-1082) — explicit guard. Hierarchical ordering uses Ward
  linkage on z-matrix (1096).
- Count metrics split out from continuous by `m.startswith('n_')` prefix
  convention (lines 420-421, 465).
- "Significant" = enrichment != 'mixed' AND q_value < 0.05 (lines 1587-1590).
- pdf.fonttype=42 / svg.fonttype='none' for editable text (82-83, 1298-1299),
  re-applied after `plt.style.use` resets them (comment line 1297).

## 6. Assumptions (checkable statements)
FACT:
- `df['cluster']` exists and is integer-castable (`int(c.strip())` at 1554,
  1632); `--clusters`/`--compare-clusters` are comma-separated ints.
- Annotated TSV is one row per read with the listed metric columns; missing
  columns are tolerated (each plotter guards `if metric not in df.columns`).
- cluster_analysis TSV keys on `cluster_id` and has `enrichment` (hard-coded
  at 1312, 1325, 1550) and optionally `q_value` (1548).
- Enrichment strings contain substrings `'E6E7'` / `'primary'` for color/sort
  routing (lines 448-451, 752-755, 1089, 1125-1127) — hard-coded to a specific
  E6E7-vs-primary experiment.
- Fonts 'Arial'/'Helvetica' available (PUBLICATION_RCPARAMS line 65).
- Auto-discovery assumes the `.read_assignments.annotated.tsv` naming convention
  (line 1320).

## 7. Dependencies (external libs; karyoplot usage; inter-script/shared deps; external tools)
FACT external libs: `argparse`, `os`, `sys` (unused), `typing`, `matplotlib`
(pyplot, patches, backend_pdf.PdfPages), `numpy`, `pandas`, `seaborn`,
`scipy.stats`, `scipy.cluster.hierarchy` (lazy import line 1094).
FACT karyoplot usage: **NONE.** Does not import karyoplot at all (imports
lines 42-53). All style/stat/format helpers are local reimplementations.
FACT inter-script deps: consumes the *annotated reads* TSV (produced upstream by
KaryoScope_sequence_annotate.py / a `read_assignments.annotated` producer) and
the cluster_analysis TSV produced by KaryoScope_cluster_analysis.py. No imports
from sibling scripts.
FACT external tools: none.

## 8. Proposed home in new layout
ASSESSMENT:
- **Subcommand: `karyoscope-analysis cluster-diagnostics`.** Thin
  `commands/cluster_diagnostics.py` (click; map argparse flags 1230-1271) →
  `core/cluster_diagnostics.py` for orchestration.
- **Decomposition:** `core/io/annotated_reads.py` (`load_annotated_data` +
  enrichment/q-value join + auto-discovery 1309-1327); `core/cluster_metrics.py`
  (pure metric/stat computations: `compute_effect_size`, `compute_cluster_summary`,
  chi-square count binning, z-scored median matrix builder); plotters stay in
  `core` but call shared theme/stat helpers.
- **karyoplot push-down candidates (highest value):**
  1. Stat/format helpers → already exist in `karyoplot.mpl.style.sig_label`
     and `karyoplot.mpl.statistics.apply_fdr`. Replace local `format_pvalue_stars`
     (113) / FDR-style q-value logic; push `format_pvalue` (101) and
     `compute_effect_size` (127, rank-biserial) DOWN into `karyoplot.mpl.statistics`.
  2. Style/rcParams + dark-mode block → `karyoplot.mpl.style.apply_default_style`
     / `fg_color` already cover lines 705-716, 886-895, 1276-1299; delete local
     `PUBLICATION_RCPARAMS` (63) + dark dicts.
  3. The z-scored median heatmap (`pub_summary_heatmap` 1046) and the box-with-
     significance-bracket comparison (`pub_comparison_figure` 684) overlap with
     `karyoplot.mpl.heatmap.plot_heatmap` and `karyoplot.mpl.comparison`
     (plot_dot_strip/lollipop) — push a generic "metric matrix heatmap" and
     "grouped boxplot + sig bracket" primitive DOWN. Also adopt
     `karyoplot.mpl.style.save_fig` (SVG+PNG) to replace ad-hoc `fig.savefig`.
- **Overlap with cluster_analysis / shared core module:** ASSESSMENT: The two
  scripts do NOT share clustering logic. cluster_analysis OWNS clustering +
  k-selection (silhouette/calinski/davies-bouldin/composite-knee, lines 34/83-91)
  and emits `enrichment`/`q_value` + the read_assignments table (read_length,
  centroid_distance per read, lines 2176-2199). diagnostics CONSUMES those.
  The genuine shared surface is the *statistics/plotting* layer (Mann-Whitney,
  rank-biserial, FDR, sig-stars, dark-mode style, heatmap/box primitives), which
  should live in `karyoplot`, not a private analysis core module. Recommend:
  share via karyoplot, keep diagnostics a distinct subcommand. They are NOT a
  consolidation-into-one-command candidate.

## 9. Smells / risks / dead code / duplication (line-cited; esp. vs cluster_analysis)
FACT:
- **Dead function:** `plot_cluster_vs_others` (570) defined, never called.
- **Dead flag:** `--no-exploratory` (1270) parsed but never consulted.
- **Unused imports:** `sys` (44), `Optional` (45).
- **Duplicate enrichment load:** cluster_analysis TSV `pd.read_csv` twice
  (1311-1312 and 1547-1550) building `enrichment_map` redundantly; the first
  pass also doesn't capture q_value, the second does.
- **Format flag half-honored:** `--format` only affects pub outputs (1661, 1679,
  1699); exploratory + comparison outputs hard-coded `.pdf` via PdfPages
  (1341, 1496, 1565, 1602).
- **In-place df mutation:** adds `supplementary_contribution` column to the
  shared df (1448) and `enrichment` (1313/1326); page 6 relies on page 5 having
  run (line 1471 checks the derived column) — ordering coupling.
- **Broad `except Exception` swallowing** around all pub figures (1666, 1684,
  1704, 1715) prints `Error:` and continues — silent partial output.
- **Hard-coded experiment semantics:** `'E6E7'`/`'primary'` substrings and
  fixed hex colors (86-98, 448-453, 752-757, 1125-1130) — not generalizable.
- **Seaborn `palette=` without `hue=`** (e.g. 207-210, 235, 264) is deprecated
  in modern seaborn and will warn/break.
- **`np.random.uniform` jitter without seed** (803) → non-reproducible figures.
- **Duplication vs cluster_analysis:** rank-biserial / Mann-Whitney logic appears
  in BOTH (diagnostics 127-152, 657-662; cluster_analysis imports mannwhitneyu at
  1912 for per-group enrichment). Sig-star + FDR logic duplicated vs
  `karyoplot.mpl.style.sig_label` / `statistics.apply_fdr`. Dark-mode rcParams
  duplicated vs `karyoplot.mpl.style`. `compute_cluster_summary` (376) overlaps
  the per-cluster `*_count/*_pct` records cluster_analysis builds at 1906-1909.
- **`format_pvalue_stars` vs karyoplot `sig_label` divergence:** diagnostics adds
  a 4th tier `****` (p<0.0001, line 115) that karyoplot's `sig_label` lacks —
  behavior must be reconciled on push-down.

## 10. Testability notes
ASSESSMENT:
- **Prime pure-function unit targets:** `format_pvalue` (101), `format_pvalue_stars`
  (113), `compute_effect_size` (127) — deterministic, table-driven; assert
  boundary p-values and the rank-biserial sign/magnitude on known arrays.
  `compute_cluster_summary` (376) — assert column set + values on a tiny df.
- **Refactor-then-test:** the chi-square count binning (929-1018) and z-scored
  median-matrix construction (1076-1083, incl. std=0 guard) should be extracted
  to pure functions returning arrays/DataFrames, then unit-tested independently
  of matplotlib.
- **Integration-only:** the figure builders (`pub_*`, `plot_*`) return Figures /
  write PDFs; test via `Agg` backend smoke tests + golden checks on the emitted
  `.cluster_summary.tsv` / `.pub_statistics.tsv` rather than pixels.
- **Fixtures:** small synthetic annotated TSV (2-3 clusters, `sample`/`group`/
  `enrichment` cols, all metric columns) + a matching cluster_analysis TSV with
  `cluster_id`/`enrichment`/`q_value`. Add a "missing optional columns" fixture
  to exercise the `if metric not in df.columns` guards.

## 11. Open questions for the user
1. Is the `'E6E7'`/`'primary'` enrichment coloring/sorting specific to one
   experiment, or should enrichment categories + colors come from config /
   the cluster_analysis output so this generalizes?
2. Should `--format` also control exploratory/comparison outputs (currently
   PDF-only), and should we emit SVG+PNG via `karyoplot.mpl.style.save_fig`?
3. Reconcile the 4-tier `****` sig stars with karyoplot's 3-tier `sig_label`
   — which threshold scheme is canonical?
4. Is `plot_cluster_vs_others` (570) and `--no-exploratory` (1270) intended
   functionality that was never wired up, or safe to delete?
5. Confirm the upstream producer/naming of `.read_assignments.annotated.tsv`
   so the IO loader + auto-discovery can be made robust.
6. Should the per-cluster summary (`compute_cluster_summary`) be unified with the
   per-cluster records cluster_analysis already emits, to avoid two definitions?
