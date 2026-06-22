# Audit: KaryoScope_cluster_plot.py (7213 lines)

## 1. Purpose (1–3 sentences)
Renders a publication-quality SVG (optionally PNG) of representative reads from each cluster produced by `KaryoScope_cluster_analysis.py`, drawing per-read feature bars (chromosome/subtelomeric/region/repeat/fiberseq tracks) annotated with cluster dendrograms, enrichment bubbles/grids, sample×cluster count matrices, and color legends. It supports both a "horizontal" layout (reads as columns, dendrogram on top) and a "vertical" layout (reads as rows, dendrogram on left, matrix sidebars), plus a separate "structural" multi-panel per-chromosome mode auto-detected from the feature matrix `mode` flag.

## 2. CLI surface
FACT. `argparse` (`parse_args`, line 77) with `RawTextHelpFormatter`. No subcommands; single script. `parse_args` also snapshots all non-None defaults into the module global `_argparse_defaults` (lines 271–276) for the parameter-table printer.

Required:
- `--cluster-analysis-prefix` (dest `cluster_prefix`, line 84) — prefix that auto-discovers `.sequence_assignments.tsv`, `.feature_matrix.npz`, `.sample_metadata.tsv`, `.cluster_analysis.tsv`.
- `--output` (line 86) — output SVG path.
- `--colors` (dest `colors_dir`, line 97) — directory of `{database}.{featureset}.colors.txt` files.

Data sources (one of `--bed` or `--input-bed-prefix` effectively required, enforced at lines 5762–5799):
- `--bed` (nargs='+', line 90), `--input-bed-prefix` (line 92), `--database` (line 95), `--featuresets` (default `chromosome,subtelomeric,region`, line 99), `--custom-beds` (`NAME:PATH`, line 101), `--density-featuresets` (line 104), `--density-bin-size` (default 300, line 107), `--density-line-plot` (`fsA:fsB`, line 109), `--rect-plot` (line 113), `--fiberseq` (dir auto-discovery, line 119).

Display/layout: `--background {white,black,both}` (default black, line 125), `--png` (line 128), `--dendro-cut` (line 130), `--bar-width` (8), `--bar-spacing` (0), `--read-spacing` (12), `--cluster-spacing` (30), `--ratio` (1/300, line 141), `--oversample` (1), `--smoothness` (default `smoothed`, line 147), `--feature-mode {smooth,transition}` (default transition, line 149), `--min-feature-width` (0.5, line 153; help text says default 1.0 — mismatch), `--min-width-exclude` (default `['novel','*arm*','ct*']`, line 155), `--target-width`/`--target-height` (auto-ratio), `--font-family` (line 268).

Mode flags: `--show-dendrogram`, `--hide-brackets`, `--no-reorder`, `--hide-dendrogram`, `--full-dendrogram`, `--fresh-dendrogram`, `--fresh-cluster-dendrogram`, `--dendro-linkage {ward,average,complete,weighted,single}` (default average), `--vertical`, `--show-matrix`/`--sample-count-matrix`, `--show-bar-plots`, `--column-tracks`, `--show-group-matrix`, `--show-group-enrichment`, `--enrichment-grid`/`--sample-enrichment-grid`, `--enrichment-normalization {raw,telomeric,total}` (default raw; see §9), `--total-reads-file`, `--orient-telomere-top` (see §9), `--show-read-indices`, `--show-clade-id`, `--show-clade-count`, `--show-cluster-numbers`, `--hide-read-labels`, `--show-threshold`, `--structural-threshold`/`--st` (0.25).

Read selection: `--reads-file` (line 161), `--curated-reps` (line 208), `--cluster-labels` + `--label-column` (default `curated_annotation`, lines 190/193), `--n-per-cluster` (dest `max_reps`), `--use-centroids` (line 212), `--priority-samples` (line 231).
Filtering: `--max-qvalue` (line 167), `--filter-enrichment` (line 169).
Logging: `--log-file/--no-log-file` (`BooleanOptionalAction`, default True, line 233).

## 3. Inputs & outputs
FACT — Inputs (all keyed off `--cluster-analysis-prefix`, lines 5736–5739):
- `{prefix}.sequence_assignments.tsv` (required, line 5757) — TSV with columns `sequence` (or `read`, renamed at lines 896/2790), `sample`, `cluster`, `group`, `enrichment` (derived), `read_length`, `rank`, `centroid_distance`; in structural mode also `chromosome`, `cluster_type`, `raw_divergence`, `norm_divergence` (lines 5451–5524).
- `{prefix}.cluster_analysis.tsv` (optional) — per-cluster stats; columns consumed: `cluster_id`, `enrichment`, `*_pct`, `*_pval`, `*_odds`, `*_count`, `p_value`, `q_value`, `odds_ratio`, `size`, `centroid_read`, plus `{group}_pval`/`{group}_odds` (lines 469–528, 3286–3297). Also `region_top`/`subtelomeric_top`/`repeat_top` parsed by `parse_top_features` (line 540) — though that function's output appears unused by the live selection path (see §9).
- `{prefix}.feature_matrix.npz` (optional) — keys consumed: `mode`, `cluster_linkage`, `cluster_ids_ordered`, `cluster_centroids`, `above_cut_linkage`, `above_cut_cluster_order`, `adj_matrix`, `read_names` (lines 1377–1503, 1614–1615, 6039–6045).
- `{prefix}.sample_metadata.tsv` (required when `--input-bed-prefix`, line 5776; also read directly at lines 5780, 6276) — columns `sample`, `group`, plus color/display columns consumed by `karyoplot.core.sample_metadata`.
- Per-sample feature BED files at `{prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed[.gz]` (lines 588–589, 1141); columns `read_id/scaffold, start, stop, feature` (4-col BED). Custom/fiberseq BEDs same 4-col format (line 1219).
- `--cluster-labels`, `--curated-reps`, `--total-reads-file`, `--reads-file`, color `.colors.txt` files.

FACT — Outputs:
- One or two SVGs (`black`/`white` themes; `{out_base}_{bg}.svg` when `both`, else `--output`) saved via `d.save_svg`.
- Optional PNG per SVG via `_svg_to_png` (`karyoplot.svg.export.svg_to_png`, line 59) when `--png`.
- `{output%.svg}.log` when `--log-file` (TeeLogger, lines 5283/5732).
- Structural mode: per-chromosome `{out_base}{theme}.{chrom}.svg` plus a combined `{out_base}{theme}.all_chromosomes.svg` grid (lines 5673–5714); may create a `*.FIRE_LINKER.bed` file on disk in the fiberseq dir (lines 5840–5849).

## 4. Pipeline / control flow
FACT — `main()` (line 5720):
1. Parse args, set module `FONT_FAMILY` (5721–5723); set up `TeeLogger` (5726–5732).
2. Derive the four prefix paths (5736–5739). Probe `feature_matrix.npz` for `mode == 'structure'`; if so call `plot_structural_mode(args, fm)` and exit (5742–5754).
3. Resolve `sample_bed_paths` + `database` from `--bed` (`parse_bed_paths`, line 1047) or `--input-bed-prefix` (reads metadata for sample list, builds `{...}/telogator/1/KaryoScope/{db}` dirs) (5762–5803).
4. Parse featuresets; merge custom-bed names; process `--fiberseq` dir (auto-discovers FIRE/LINKER/m6A/5mC, builds combined FIRE_LINKER bed and density-line config) (5811–5863). Pull density/rect/density-line placeholder tracks out of the main featureset list and append synthetic track names (`density_line`, `rect_plot`) (5867–5897).
5. Load data: `load_sample_metadata` (374), `load_cluster_analysis` (424), `load_cluster_labels` (390), `load_representative_reads` (874), `load_feature_matrix` (1031), `load_color_files` (1094), `load_bed_data` (1128), `load_custom_bed_files` (1191) (5922–5984).
6. Optionally compute fresh read/cluster dendrograms (`compute_fresh_dendrogram` 1659, `compute_fresh_cluster_dendrogram` 1742) (5986–6011).
7. Generate sample/enrichment colors (`generate_sample_colors` 1323, `get_enrichment_colors` 1241) (6013–6028).
8. Compute cluster ordering: prefer `above_cut_linkage` subtree extraction inline (6037–6153) else `compute_cluster_dendrogram_order` (1354, which round-trips scipy linkage → Newick via Bio.Phylo to prune to displayed clusters). Apply `--dendro-cut` via `fcluster` (6160–6183). Honor `--hide-dendrogram`/`--full-dendrogram` (6185–6200).
9. Branch: `if args.vertical:` → vertical drawing path (6205–6748, ends with `return`); else horizontal path (6750–7209).
   - Vertical: compute left-margin budget for dendrogram + enrichment bubbles/grid + sample/group/enrichment matrices (6216–6306); compute y-positions per read/cluster (6308–6353); rasterize features via `rasterize_features` (6419); per-theme loop draws dendrogram, scale bar, feature bars (row or `draw_feature_bars_column_mode`), matrices, enrichment bubbles/grid, vertical cluster labels, right-side vertical legends; save (6489–6731).
   - Horizontal: compute group_width, top_margin (dendrogram+grid+bracket), x-positions per read (6765–6832); rasterize features incl. density-line/rect-plot synthetic tracks (6868–7023); compute `read_heights` for borders (7026–7043); per-theme loop draws dendrogram header, horizontal enrichment grid, cluster brackets, annotation bars, feature bars, read labels, per-read featureset labels, top legends + bottom color legends; save (7065–7193).
10. Emit uncolored-feature warnings + summary + `_print_params_and_command` (279) (7195–7209).

Key helper clusters: data loading (374–1234), color/enrichment helpers (1241–1347), dendrogram compute (1354–1830), dendrogram/bracket drawing (1833–2333), vertical-mode drawing (2340–3754), enrichment bubbles/grids/legends (3782–4552), annotation/feature-bar/read-label/color-legend drawing (4588–5274), density helpers (5309–5394), structural mode (5400–5717).

## 5. Key design decisions
FACT (with stated rationale where present):
- Representative selection is delegated to `KaryoScope_select_representatives.py`; this script only *loads* selected reads (`load_representative_reads`, docstring lines 877–878). Multiple selection sources with a precedence chain: `--use-centroids` > `--curated-reps` > `--reads-file` > `representative_read_*` columns in `--cluster-labels` (lines 922–979).
- Cluster ordering tiers (`get_enrichment_tier`, lines 492–506): tier 0 = 100% enriched, tier 1 = ≥80%, tier 2 = other; sorted by (tier, p_value). Comment-stated as a quality heuristic.
- Per-sample mode detection filters `*_pval` candidate names against `known_samples` so `{group}_pval` columns aren't mistaken for samples (lines 472–486; rationale in docstring 432–436).
- Above-cut linkage is used directly when all clusters displayed, else an inline subtree-extraction collapses hidden leaves preserving original merge distances (6055–6142). The older `compute_cluster_dendrogram_order` path uses Bio.Phylo Newick round-trip + `optimal_leaf_ordering` (1414–1564) — explicitly to "preserve the original tree structure for displayed clusters" (1412–1413).
- Matrix heatmap uses log1p intensity scaling (`math.log1p(count)/math.log1p(max_count)`), with separate dark→yellow (black bg) / light→dark-blue (white bg) gradients (lines 2974–2996). Rationale stated: log helps visualize low values (2967–2969).
- Sample columns in the matrix are clustered within group via Ward linkage on proportion-normalized count profiles + `optimal_leaf_ordering` (2855–2923). Rationale: compares relative distributions not absolute counts (2863–2864).
- Bubble encoding (`_compute_bubble_style`, 3874): size ∝ %, color = one-sided white→red on log2(OR) clamped to ±4, alpha = binary FDR threshold OR a `-log10(p)/10` gradient when no FDR passed. The group-enrichment matrix deliberately omits `fdr=` to get per-cell within-row opacity differentiation (extensive rationale 3277–3334; mirrored in `draw_enrichment_grid` comment 4002–4006).
- `draw_enrichment_bubbles` (3782) uses a *two-sided* diverging red/blue colormap (Tumor/Normal) and binary alpha at q<0.1 (3862–3863) — a different encoding from `_compute_bubble_style`.
- Scale bar bp/label picking delegated to `karyoplot.core.coords.pick_round_scale_bp` + `format_genomic_distance` (Phase 13.2, lines 2665–2675).
- `--target-width`/`--target-height` auto-derive `ratio` after subtracting hardcoded margin estimates (400px vertical at 6376, 350px horizontal at 6862).
- Fresh dendrograms build an edge-only feature matrix (adjacent-feature co-occurrence, excluding `novel`), `log1p` + `StandardScaler` z-score, Ward linkage (1679–1719). Stated to mirror `KS_allchr_dendrogram.py`.

## 6. Assumptions (as checkable statements)
- BED files live at exactly `{input_bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed[.gz]` (lines 588, 1141, 5787). The `telogator/1/KaryoScope` path segment is hardcoded.
- BED filenames encode `{sample}.telogator.1.{database}...`, i.e. database is the 4th dot-delimited token (line 1076).
- BED files are 4-column `read/scaffold start stop feature`; lines with <4 cols or non-int coords are skipped (1172–1178).
- `sequence_assignments.tsv` has columns `sequence`/`read`, `sample`, `cluster`, `group`; structural mode additionally needs `chromosome`, `cluster_type`, `raw_divergence`, `norm_divergence`.
- `sample_metadata.tsv` has a `sample` column and a `group` column; matrix sorting specifically treats `Normal`/`Control` groups as the "first" group (lines 2798, 2878).
- Hardcoded enrichment-label conventions: labels end in `-enriched` (1259, 6026); special colors for `post`/`pre`/`Normal`/`Tumor` substrings (1282–1287, 3459–3463).
- `feature_matrix.npz` is a structure-mode file iff it contains key `mode` equal to `'structure'` (5747).
- Color file naming `{database}.{featureset}.colors.txt` (1107); missing file is fatal (`sys.exit(1)`, 1110).
- `pick_round_scale_bp` options assume read lengths roughly in the 1–20 kb range (line 2671).
- Environment provides `drawsvg`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `biopython`, and `rsvg-convert` on PATH (only when `--png`); `openpyxl` is implicitly needed for `.xlsx` cluster-labels.
- Reads are assumed ≤~20 kb for the length-tier heuristics in `_select_by_strategy` (16–20 kb "ideal", lines 760–762).

## 7. Dependencies
**External libraries (FACT):** `argparse`, `fnmatch`, `glob`, `gzip`, `os`, `subprocess` (imported line 43 but never used — see §9), `sys`, `collections`, `math`; `drawsvg`, `numpy`, `pandas`; lazy imports of `scipy.cluster.hierarchy`/`scipy.spatial.distance`/`scipy.stats.fisher_exact`, `sklearn.preprocessing.StandardScaler`, `Bio.Phylo`, `re`, `math`, `io.StringIO`.

**`karyoplot` usage — ALREADY delegated (FACT):**
- `karyoplot.core.colors.TAB10, TAB20` (line 58); `load_palette_file` (line 1100, in `load_color_files`).
- `karyoplot.svg.export.svg_to_png` re-exported as `_svg_to_png` (line 59).
- `karyoplot.svg.reads.smooth_features_to_pixels, features_to_pixels_direct, rasterize_features` (lines 60–64; comment at 2594 notes Phase 13.1 move).
- `karyoplot.core.sample_metadata.load_sample_metadata` (line 380, wrapped to legacy 4-tuple).
- `karyoplot.core.coords.pick_round_scale_bp` + `karyoplot.core.text.format_genomic_distance` (lines 2665–2666, Phase 13 P2.c).
- `karyoplot.core.text.abbreviate_read_name` (line 4823, Phase 13.2).

**NOT yet delegated (local logic that overlaps karyoplot):** all dendrogram-drawing functions; all matrix/bubble/grid/legend drawing (`draw_color_legends*`, `draw_grid_legend*`, `draw_bubble_legend*`); density binning (`compute_density_features`, `compute_density_line`); `_compute_bubble_style`; `_build_legend_items`/`_find_common_label`; `TeeLogger`.

**Inter-script / shared-file deps (FACT):** Documented to consume outputs of `KaryoScope_cluster_analysis.py` and `KaryoScope_select_representatives.py` (header lines 6–25). No imports from sibling scripts; no use of `scripts/_feature_vocab.py`; the two `scripts/*.colors.txt` files are reference color files (not imported — passed via `--colors`). No `data/` directory reads. Does not import `karyoplot.mpl.statistics` (re-implements Fisher's exact / FDR locally).

**External tools / subprocess (FACT):** `--png` calls `karyoplot.svg.export.svg_to_png` which (per the import comment) shells out to `rsvg-convert`. The top-level `import subprocess` (line 43) is unused locally.

## 8. Proposed home in new layout
**Target CLI subcommand:** `karyoscope-analysis cluster-plot` (thin wrapper in `commands/cluster_plot.py`). The auto-detected structural mode should become an explicit sibling, e.g. `cluster-plot --structural` or a separate `structural-plot` subcommand, rather than sniffing the npz `mode` key.

**Decomposition (this file MUST be split):**
- `commands/cluster_plot.py` — click command mirroring `parse_args`; calls into core; owns `TeeLogger` setup and `_print_params_and_command` (or move the latter to a shared `core/logging.py`).
- `core/cluster_plot/io.py` (`core/io/`) — `load_sample_metadata` wrapper, `load_cluster_labels`, `load_cluster_analysis`, `load_representative_reads`, `load_curated_representatives`, `load_feature_matrix`, `parse_bed_paths`, `load_color_files`, `load_bed_data`, `load_custom_bed_files`, `get_read_features`.
- `core/cluster_plot/selection.py` — `parse_top_features`, `score_read_features`, `compute_balanced_score`, `select_fallback_read`, `_select_by_centroid`, `_select_by_strategy`, `_select_fallback` (note: these appear largely dead, §9 — confirm before keeping).
- `core/cluster_plot/dendrogram.py` — `compute_cluster_dendrogram_order` (incl. Bio.Phylo Newick round-trip), `compute_full_dendrogram`, `compute_fresh_dendrogram`, `compute_fresh_cluster_dendrogram`, and the inline above-cut subtree extraction currently embedded in `main()` (6055–6142).
- `core/cluster_plot/colors.py` — `get_enrichment_colors`, `get_cluster_colors`, `get_primary_color`, `generate_sample_colors`.
- `core/cluster_plot/density.py` — `compute_density_features`, `compute_density_line`.
- `core/cluster_plot/layout.py` + `render_horizontal.py` + `render_vertical.py` + `render_structural.py` — split the 1500-line `main()` and the structural path; each render module owns its drawing orchestration.
- `core/cluster_plot/draw/` — the many `draw_*` helpers grouped (dendrogram, matrix, bubbles/grids, legends, feature bars), thin over karyoplot.

**Push DOWN into `karyoplot` (and which module):**
- SVG dendrogram drawing (`draw_dendrogram`, `draw_cluster_dendrogram`, `draw_cluster_dendrogram_vertical`, `draw_full_dendrogram[_header]`, `draw_dendrogram_scale_axis`, `draw_mini_dendrogram`) → new `karyoplot.svg` module (e.g. `dendrogram.py`); the leaf-ordering/linkage-subtree math overlaps `karyoplot.mpl.heatmap.fix_leaf_ordering`/`push_leaves_to_edge`/`cluster_and_reorder` — fold the matrix manipulation there into `karyoplot.core`.
- Enrichment bubble/grid + `_compute_bubble_style` + legends → `karyoplot.svg.legend` already hosts `draw_grouped_legend`/`merge_by_color`; the bubble-style colormap and `_find_common_label`/`_build_legend_items` belong there.
- Heatmap matrix (`draw_sample_matrix`/`draw_group_matrix`) log1p color gradient → `karyoplot.core.colors` (gradient helper) + a `karyoplot.svg` matrix drawer.
- Fisher's exact + FDR enrichment math → `karyoplot.mpl.statistics` (already has `compare_two_conditions`, `apply_fdr`).
- `TeeLogger` → `karyoplot.core` or a shared analysis `core/logging.py`.

## 9. Smells / risks / dead code / duplication
- **`main()` is ~1490 lines (5720–7209)** with two near-duplicate layout engines (vertical 6205–6748, horizontal 6750–7209). Read-position loops, ratio auto-calc, uncolored-feature warning, and summary printing are duplicated across the two branches (e.g. 6308–6353 vs 6795–6832; 6734–6746 vs 7196–7209).
- **Likely dead selection code:** `parse_top_features` (540), `score_read_features` (609), `get_read_features` (569), `compute_balanced_score` (655), `select_fallback_read` (679), `_select_by_strategy` (730), `_select_fallback` (791) — none are called from `main()`; live selection is `_select_by_centroid` (called at 1012) plus file-based filters. Confirm before deleting.
- **Unused import:** `import subprocess` (line 43) never used in this file.
- **Default/help mismatch:** `--min-feature-width` default is `0.5` but help says "default 1.0" (lines 153–154).
- **`'cluster_dendro_data' not in dir()` (line 6032)** is a fragile way to detect a prior local assignment; brittle and hard to follow. Replace with explicit init.
- **Three different bubble colormaps / alpha conventions** (`draw_enrichment_bubbles` diverging+binary 3847–3863; `_compute_bubble_style` one-sided+FDR/gradient 3918–3927; matrix gradients) risk inconsistent semantics across panels.
- **Advertised-but-unimplemented flags:** `--enrichment-normalization {telomeric,total}` and `--total-reads-file` are parsed (252–260) but never consumed (the normalization is never applied). `--orient-telomere-top` (261) parsed but never used. `--show-clade-id`/`--show-clade-count` (236–239) parsed but never used. `--show-threshold`/`--structural-threshold` only affect structural mode. These are documented behaviors that don't happen.
- **`_lev` (Levenshtein) defined inside the `chrom` loop** in `plot_structural_mode` (5560) — recompiled per chromosome, O(n²) pairwise on encoded strings; fine for small n but a perf/clarity smell.
- **Broad `except Exception: pass`** swallowing in many loaders (603, 919, etc.) hides malformed-input bugs; several others `traceback.print_exc()` then fall back silently.
- **`plot_structural_mode` writes a derived `*.FIRE_LINKER.bed` into the user's fiberseq dir** (5840–5849) — a side-effecting file creation hidden in a plotting command.
- **Hardcoded magic numbers everywhere** (margins 400/350/200, panel_width 1200, n_cols 5, color stops) — should become named constants/params.
- **`_svg_to_png` import comment says "re-exported for legacy call sites"** (line 59) but the only caller is local; the legacy alias is unnecessary.
- **Reading large TSVs twice:** `load_representative_reads` reads `representatives_file`, then `draw_sample_matrix` re-reads the same file (2789) and `--input-bed-prefix` path re-reads metadata (5780, 6276).
- **Bio.Phylo Newick round-trip** (1429–1553) to prune a scipy tree is heavyweight and fragile (non-binary node handling, `.6f` branch-length formatting); the inline above-cut extraction (6094–6142) is a cleaner reimplementation of the same idea — two code paths solving one problem.

## 10. Testability notes
**Pure functions easily unit-testable (no I/O, deterministic):**
- `parse_top_features` (540), `compute_balanced_score` (655), `_compute_bubble_style` (3874), `_find_common_label` (5007), `_build_legend_items` (5074), `compute_density_features` (5309), `compute_density_line` (5358), `get_enrichment_colors` (1241), `get_cluster_colors` (1292), `get_primary_color` (1304), `generate_sample_colors` (1323), the inner `_lev` (5560).
- Dendrogram compute functions (`compute_full_dendrogram`, `compute_fresh_dendrogram`, `compute_fresh_cluster_dendrogram`, `compute_cluster_dendrogram_order`) are testable with small synthetic feature matrices / cluster_reads dicts — verify ordering + linkage shape (they print/return rather than draw).

**Needs fixtures:** loaders (`load_cluster_analysis`, `load_representative_reads`, `load_curated_representatives`, `load_color_files`, `load_bed_data`, `parse_bed_paths`) need small TSV/BED/npz/colors fixtures and the exact directory layout.

**Integration / golden-image only:** every `draw_*` function and the two `main()` render branches + `plot_structural_mode` — these append `drawsvg` primitives and `save_svg`. Best covered by golden SVG snapshot tests (string-stable since the repo already does "byte-identical" benchmarks per the git log) across representative flag combinations (vertical+matrix, horizontal+grid, full-dendrogram, structural, fiberseq). `TeeLogger` and `_print_params_and_command` are side-effecting (stdout/file).

## 11. Open questions for the user
- Should the auto-detected **structural mode** become its own subcommand (`structural-plot`) or a `--structural` flag? It has a completely separate code path, inputs, and outputs.
- Are the **feature-scoring selection functions** (`parse_top_features`, `score_read_features`, `_select_by_strategy`, etc.) truly dead, or invoked by an external caller / planned re-enable? Confirm before dropping them in the new layout.
- The **`--enrichment-normalization telomeric/total`, `--total-reads-file`, `--orient-telomere-top`, `--show-clade-id/-count`** flags are parsed but never used. Implement, or remove from the CLI?
- Is the hardcoded **`telogator/1/KaryoScope/{database}` BED path layout** stable across all datasets, or should it be configurable?
- Do you want the **two layout engines (vertical/horizontal) unified** behind one renderer, or kept separate? They share substantial logic but differ enough that unification is a real effort.
- How much of the drawing should move into **`karyoplot.svg`** vs. staying analysis-specific? Specifically: dendrogram drawing, bubble/grid drawing, and the heatmap-matrix gradient are the strongest shared-library candidates — confirm scope.
- Confirm the **`--min-feature-width` intended default** (code 0.5 vs help 1.0).
- Is **writing the derived `*.FIRE_LINKER.bed`** into the input directory acceptable, or should it go to a temp/output location?
- Should **enrichment statistics (Fisher/FDR)** be unified on `karyoplot.mpl.statistics` rather than the local re-implementations, to keep p-values consistent with `cluster_analysis`?
