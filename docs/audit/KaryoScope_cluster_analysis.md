# Audit: KaryoScope_cluster_analysis.py (3103 lines)

## 1. Purpose
FACT: Per the module docstring (lines 1-3), the script "Analyzes hierarchical clustering results to identify biologically interesting clusters with sequence assignments sorted by centroid distance for visualization."

FACT: It loads one or more KaryoScope feature BED files, builds a per-read feature/edge adjacency matrix, optionally reduces dimensions with truncated SVD, performs hierarchical (Ward by default) clustering, auto-selects the number of clusters `k`, tests each cluster for group/sample enrichment (Fisher's exact + FDR), and emits TSVs, an NPZ matrix bundle (consumed downstream by `KaryoScope_cluster_plot.py`, line 3042-3044), and a battery of diagnostic plots (k-selection, knee, dendrogram, circular dendrogram, UMAP, enrichment bubble/heatmap).

FACT: A second, additive analysis mode (`structure`, lines 221-417) instead clusters reads per-chromosome by Levenshtein distance of feature sequences to flag structural outliers vs a "Major" consensus, then `sys.exit(0)` (line 417).

ASSESSMENT: This is the statistical/clustering core of the analysis repo and the producer of the artifacts most other scripts consume.

## 2. CLI surface
FACT: Uses `argparse` (top-level module code, lines 57-188), `RawTextHelpFormatter`. Arguments:

- `--bed` (required, nargs='+', line 61): input BED file(s), gzip ok, concatenated.
- `--output-prefix` (required, line 64).
- `--sample-metadata` (line 66): TSV `sample, group, color`.
- `--comparison-mode` {two-group, per-sample} default two-group (line 70).
- `--control-group` (line 75): reference group, auto-detect default.
- `--n-clusters` (int, line 77): fixed k; default None = auto.
- `--min-k` (40), `--max-k` (300) (lines 79, 81).
- `--k-selection` {composite, silhouette, cosine-silhouette, calinski, davies-bouldin, composite-knee} default `composite-knee` (line 83).
- `--min-cluster-size` (3, line 92).
- `--min-sequence-length` (10000), `--max-sequence-length` (50000) (lines 94, 96).
- `--sequence-list` (line 98): whitelist of read names.
- `--exclude-features` default `"novel,canonical_telomere*"` (line 100), comma-list, wildcard, colon-component match.
- `--linkage-method` default ward (line 102).
- `--matrix-type` {binary, count, length_weighted, count_log1p, count_log1p_zscore, count_log1p_zscore_blockweight} default `count_log1p_zscore_blockweight` (line 104).
- `--edges` {directional, symmetric} default symmetric (line 114).
- `--matrix-mode` {layered, combined} default combined (line 119).
- `--include-edges` / `--no-include-edges` default True (BooleanOptionalAction, line 124).
- `--abundance` / `--no-abundance` default True (line 127).
- `--umap` / `--no-umap` default True (line 130); `--umap-neighbors` (25), `--umap-min-dist` (0.2), `--umap-html` (False) (lines 136-142).
- `--circular-dendrogram` default True (line 133).
- `--perfect-threshold` (1.0), `--strong-threshold` (0.80) (lines 143, 145).
- `--early-stopping` (150 iterations, line 147).
- `--silhouette-sample-size` (2000, line 149).
- `--enrichment-normalization` {raw, telomeric, total} default raw (line 151); `--total-reads-file` (line 157).
- `--reduce-dims` (int, 500; 0 disables, line 160).
- `--background` {white, black, both} default white (line 164).
- `--log-file` / `--no-log-file` default True (line 170).
- `--fdr-threshold` (0.05), `--fdr-method` {bh, by} default bh (lines 173, 175).
- `--analysis-mode` {enrichment, structure} default enrichment (line 180); `--structural-threshold`/`--st` (0.25, line 185).

FACT: Cross-arg validation: at least one of `--edges`/`--abundance` required (lines 209-210).

ASSESSMENT: 40+ flags on a single top-level `argparse` parser. There is no subcommand structure; `analysis-mode=structure` is effectively a hidden second program bolted on via a branch.

## 3. Inputs & outputs
FACT inputs:
- BED file(s) parsed by `load_bed_file` (lines 474-499): columns `read/start/end/feature[/chromosome]`; chromosome defaults to "unknown" if <5 cols. Sample label from filename prefix via `extract_sample_label` (lines 501-508).
- `--sample-metadata` TSV loaded via `load_sample_metadata` (lines 597-609) which delegates to `karyoplot.core.sample_metadata.load_sample_metadata`.
- `--total-reads-file` TSV (`sample`, `total_reads`) read directly with pandas (lines 952-958).
- `--sequence-list` newline file (lines 853-854).

FACT outputs (enrichment mode):
- `{prefix}.log` (TeeLogger, lines 191-214).
- `{prefix}.cluster_k_analysis.tsv` (line 1586).
- `{prefix}.cluster_analysis.tsv` (line 2165).
- `{prefix}.sequence_assignments.tsv` (line 2201).
- `{prefix}.feature_matrix.npz` (line 2207-2232): adj_matrix, seq_names, read_names (dup), cluster_labels, cluster_centroids, cluster_linkage, above_cut_linkage, optional SVD components/feature names.
- `{prefix}.sample_metadata.tsv` (line 2238).
- `{prefix}.umap_coordinates.tsv` (line 2706) and `{prefix}.umap.html` (line 2820).
- PDFs: `.svd_scree`, `.k_selection`, `.composite_knee_diagnostic`, `.cluster_analysis`, `.circular_dendrogram`, `.umap`, `.enrichment_bubble`, `.{group|sample}_percentage` (with optional `_dark` suffix).

FACT outputs (structure mode, lines 404-416): `{prefix}.sequence_assignments.tsv` (different schema: read, chromosome, cluster, cluster_type, divergence metrics) and `{prefix}.feature_matrix.npz` (mode=structure, per-chrom linkages).

## 4. Pipeline / control flow (key functions + line numbers)
FACT: The script is mostly **top-level module code**, not functions. Control flow:
1. Imports + font registration (lines 21-54); argparse build+parse (57-188); logging setup (191-214).
2. Helper/`def`s declared before the procedural body: `run_structure_mode` (221-417), `get_plot_style`/`apply_plot_style`/`get_backgrounds_to_generate` (420-464), `load_bed_file` (474-499), `extract_sample_label` (501-508), `get_edges` (511-536), `detect_feature_layers` (539-556), `split_feature_into_layers` (559-577), `get_layer_features` (580-595), `load_sample_metadata` (597-609), `generate_group_colors` (612-642), `calculate_enrichment_two_group` (645-712), `calculate_enrichment_per_sample` (715-824).
3. Procedural body begins line 827: load+concat BEDs (832-847); sequence-list filter (850-862); length/span calc (865-871); span filter (876-891); excluded-feature filter (894-921); annotated length (925-926).
4. **Branch** (lines 928-932): if `analysis_mode == 'structure'` -> `run_structure_mode` -> exit.
5. Metadata + group/color setup (934-1005).
6. Matrix build: `get_weighted_edges` (1024-1039), `_build_matrix_from_features` (1042-1184, the shared builder), `build_layer_matrix` (1187-1213), `build_combined_matrix` (1216-1231); dispatch by matrix-mode (1234-1296).
7. Optional truncated-SVD reduction + scree plot (1309-1397).
8. `pdist`(euclidean) + `linkage` (1401-1402).
9. k-selection: `fast_enrichment_check` (1412-1450); per-k loop computing silhouette/cosine-silhouette/CH/DB + enrichment counts + composite (1472-1582, with early stopping); k-selection plot (1591-1697); knee detection (1699-1726); knee plot (1729-1790); selected_k chosen by metric (1802-1832).
10. `fcluster` cut (1839); per-cluster enrichment + centroid (1850-1964).
11. Cluster-level dendrogram + above-cut linkage remap incl. recursive `get_cluster_for_node` (1966-2086).
12. FDR via `false_discovery_control` + `update_enrichment_label` (2088-2122); console summary table (2125-2162).
13. Save TSV/NPZ/metadata (2164-2246).
14. Plots: `get_cluster_style` (2258-2264); cluster_analysis 2x2/3x2 (2281-2448); circular dendrogram (2451-2576); UMAP static + Plotly (2579-2829); bubble + percentage heatmap (2831-3015).
15. Summary, parameters table (`fmt_param`, 3053-3097), echo command (3017-3103).

## 5. Key design decisions (cite lines; clustering/linkage/normalization/k-selection)
FACT — Linkage: Default `ward` (line 102); distance is `pdist(adj_matrix, metric='euclidean')` (line 1401), then `linkage(dist_matrix, method=args.linkage_method)` (line 1402). ASSESSMENT: Ward+Euclidean is consistent (Ward assumes Euclidean), but note Ward is run on the *SVD-reduced* matrix when `--reduce-dims>0`.

FACT — Matrix construction (`_build_matrix_from_features`, 1042-1184): two blocks, an **edge** block (feature transitions; symmetric mode sorts the pair alphabetically so A->B == B->A, lines 1065-1073) and an **abundance** block (per-feature summed bp, lines 1125-1132). matrix-type transforms:
- binary: 0/1 presence (1082-1087, 1120-1123).
- count: raw transition counts / raw bp totals.
- length_weighted: edge weight = avg adjacent feature length / read_len (1096-1101), abundance = total_len/read_len (1134-1140).
- count_log1p: `np.log1p` on both blocks (1103-1104, 1142-1143).
- count_log1p_zscore: per-column `StandardScaler` on edge and abundance blocks (1145-1153).
- count_log1p_zscore_blockweight (DEFAULT): additionally rescales the abundance block by `sqrt(n_edge/n_abund)` so the two blocks contribute equal total variance (lines 1155-1162). Stated WHY: "for equal variance contribution" / "equal edge/abundance contribution" (lines 113, 1162).

FACT — Dimensionality reduction: TruncatedSVD with `n_components=min(reduce_dims, n_samples-1, n_features)`, random_state=42 (lines 1313-1323). Stated WHY: "Recommended for merged BED files which can create very high-dimensional matrices" (lines 161-162).

FACT — k-search range: `range(min_k, min(max_k+1, n_reads//10))` (lines 1456-1457). Stated WHY: "need at least 10 reads per cluster on average" (line 1455).

FACT — Per-k metrics (lines 1477-1485): silhouette + cosine-silhouette (subsampled to `silhouette_sample_size`, random_state=42, lines 1478-1483), Calinski-Harabasz, Davies-Bouldin.

FACT — Composite score (lines 1541-1546): `0.5*silhouette_norm + 0.1*enriched_ratio + 0.4*perfect_ratio`, where silhouette is mapped [-1,1]->[0,1] (line 1543). Stated WHY: "Weights favor biologically pure clusters over moderate enrichment" (line 1542).

FACT — Composite-knee (DEFAULT k-selection, lines 1699-1726): Kneedle-style. Normalize k and composite score to observed [0,1] range; `knee_distance = score_norm - k_norm`; smooth with a centered rolling mean of window `max(3, int(0.2*k_range))` (lines 1717-1720); `selected_k = argmax(smoothed)` (line 1721, used at 1816). Raw (unsmoothed) knee also computed for the diagnostic (1724-1726). Comment acknowledges knee "may vary slightly with max-k due to normalization" (line 1705).

FACT — Enrichment statistics:
- k-loop fast check (`fast_enrichment_check`, 1412-1450) uses an **odds-ratio with 0.5 pseudocount**, threshold >1.5 to call enrichment (lines 1444-1446); purity = max group fraction (1427).
- two-group (`calculate_enrichment_two_group`, 645-712): Fisher's exact two-sided on a control-vs-treatment 2x2; direction from odds ratio; significance at raw p<0.05 (lines 677-696).
- per-sample (`calculate_enrichment_per_sample`, 715-824): one-sided (`alternative='greater'`) Fisher per sample-vs-rest; for telomeric/total normalization it computes size factors (`denominator/median_denominator`) and Fisher's on **rounded scaled counts** (lines 770-795). Stated WHY: "so samples are comparable" (lines 771-772).
- FDR: `scipy.stats.false_discovery_control(method=bh|by)` over cluster p-values (line 2094); labels demoted to 'mixed' when q>=threshold (`update_enrichment_label`, 2101-2110).

FACT — Centroid: per-cluster centroid = mean vector; representative read = argmin Euclidean distance to mean (lines 1870-1875); assignments sorted by centroid distance and ranked (lines 2191-2194).

FACT — Structure mode distance: normalized Levenshtein on feature-string encodings (`raw/max_len`), Ward linkage, `fcluster(criterion='distance', t=structural_threshold)` (lines 304-320); reports also raw, binary (symmetric set diff), and length-weighted divergence vs the major consensus structure (lines 379-401).

## 6. Assumptions (checkable statements)
- FACT: BED is tab-separated; col0=read, col1/2=int start/end, col3=feature, optional col4=chromosome (lines 482-489). Lines with <4 fields are silently dropped (line 483).
- FACT: Sample label = first dot-delimited token of the basename (lines 504-507). Two files sharing that prefix would collide.
- FACT: One read maps to exactly one sample (`groupby('sequence')['sample'].first()`, line 844) — duplicate read IDs across samples take the first.
- FACT: Read span is `end.max()-start.min()` computed BEFORE feature exclusion (lines 869-871, comment lines 867-868, 874-875).
- FACT: A read needs >=2 features to contribute any edge (`get_edges` returns [] for len<=1, lines 524-525).
- FACT: Control group defaults to the alphabetically-first group when 2+ groups and not specified (lines 973-975).
- FACT: Default length filter keeps only reads with span in [10000, 50000] bp (lines 94-96, 881-882) — telomere-length-specific.
- FACT: composite-knee, the default, requires `score_max>score_min` and `k_max_obs>k_min_obs` else falls back to zeros arrays (lines 1710-1711).
- ASSESSMENT: With one BED file and no metadata, every sample is its own group and two-group mode degenerates (the `len(groups)<2` branch at lines 651-659 returns p=1.0 'mixed').

## 7. Dependencies
FACT — External libs: `numpy`, `pandas`, `argparse`, `gzip`, `fnmatch`, `collections` (lines 21-30). `scipy`: `cluster.hierarchy` (linkage, dendrogram, fcluster, set_link_color_palette), `spatial.distance` (pdist, squareform), `stats` (fisher_exact, false_discovery_control, mannwhitneyu imported inline at 1912) (lines 31-33). `sklearn`: `metrics` (silhouette_score, calinski_harabasz_score, davies_bouldin_score), `decomposition.TruncatedSVD`, `preprocessing.StandardScaler` (lines 34-36). `matplotlib` (Agg) + `matplotlib.colors` (lines 37-39). Optional at runtime: `umap` (line 2582), `plotly` (line 2713).

FACT — karyoplot usage (what's already delegated):
- `karyoplot.core.fonts.register_fonts` called twice (lines 46-48) — fonts delegated.
- `karyoplot.core.sample_metadata.load_sample_metadata` wrapped by local `load_sample_metadata` (lines 602-609), returning the legacy 3-tuple.

FACT — karyoplot facilities NOT yet used but overlapping (push-down candidates):
- `karyoplot.mpl.style.apply_default_style/fg_color/sig_label/save_fig` overlap local `get_plot_style`/`apply_plot_style`/`get_backgrounds_to_generate` (420-464) and the inline `***/**/*` logic at lines 2987-2996.
- `karyoplot.mpl.statistics.apply_fdr` (different, BH-only manual impl) overlaps the FDR step (line 2094); `compare_two_conditions` overlaps Fisher logic.
- `karyoplot.core.io.load_bed` overlaps local `load_bed_file` (474-499).
- `karyoplot.core.colors.qualitative_palette` overlaps `generate_group_colors` tab10 logic (612-642).

FACT — Inter-script / shared-file deps: Output NPZ + sequence_assignments.tsv + sample_metadata.tsv are consumed by `KaryoScope_cluster_plot.py` (lines 3042-3044). Upstream merge helper `KaryoScope_merge_beds.py` referenced in usage (lines 62-63). This script does NOT import `_feature_vocab` (confirmed: only `KaryoScope_sequence_annotate.py` and `KaryoScope_cluster_annotate.py` do). No custom color/font files read directly (delegated to karyoplot fonts).

FACT — External tools: none (no subprocess calls).

## 8. Proposed home in new layout
ASSESSMENT — Subcommand: `karyoscope-analysis cluster` (the enrichment mode). The `structure` mode is a distinct workflow and should become its own subcommand `karyoscope-analysis cluster-structure` (or `structure`) rather than a `--analysis-mode` branch.

ASSESSMENT — Decomposition:
- `commands/cluster.py` — click command: arg definitions, logging setup, orchestration; thin.
- `commands/cluster_structure.py` — the per-chromosome structural workflow.
- `core/matrix.py` — `get_edges`, `get_weighted_edges`, `detect_feature_layers`, `split_feature_into_layers`, `get_layer_features`, `_build_matrix_from_features`, `build_layer_matrix`, `build_combined_matrix`, SVD reduction (pure, NumPy-only; prime unit-test target).
- `core/clustering.py` — linkage/cut, `fast_enrichment_check`, the per-k metric loop, composite score, knee detection, `selected_k` logic, above-cut linkage remap + `get_cluster_for_node`.
- `core/enrichment.py` — `calculate_enrichment_two_group`, `calculate_enrichment_per_sample`, size-factor normalization, FDR demotion (`update_enrichment_label`).
- `core/structure.py` — Levenshtein + per-chrom outlier logic.
- `core/io/bed.py` — BED loading / sample-label extraction / sequence-list & feature filters.
- `core/io/results.py` — NPZ + TSV writers (cluster_analysis, sequence_assignments, sample_metadata, umap_coordinates).
- `plots/` (or `core/plots/`) — scree, k-selection, knee, dendrogram, circular dendrogram, UMAP, bubble, percentage figures.

ASSESSMENT — karyoplot push-down candidates (module):
1. Plot style/theme + significance stars + dual-background save -> `karyoplot.mpl.style` (extend `apply_default_style`/`save_fig`; `sig_label` already exists).
2. BED loader / sample-label extraction -> `karyoplot.core.io`.
3. Group/sample color generation -> `karyoplot.core.colors` (`qualitative_palette`).
4. Generic Fisher-exact + size-factor normalization + FDR helpers -> `karyoplot.mpl.statistics` (or `karyoplot.core` if mpl-free).
ASSESSMENT: The clustering/k-selection math is analysis-specific and should stay in `karyoscope_analysis.core`, NOT pushed to karyoplot.

## 9. Smells / risks / dead code / duplication (line-cited)
- RISK (significance threshold inconsistency): k-loop calls enrichment at odds-ratio>1.5 (lines 1444-1446); two-group/per-sample use raw p<0.05 (lines 688, 798); final labels use FDR q<threshold (lines 2103-2108). Three different significance definitions across the pipeline.
- RISK (FDR on the SAME p-values for both metrics): in per-sample mode `cluster_df['p_value']` is the *min* sample p-value (line 819) and `false_discovery_control` is applied over clusters (line 2094), not over the sample x cluster grid — multiple-testing across samples within a cluster is not corrected.
- RISK (one-sided vs two-sided mismatch): per-sample Fisher is one-sided 'greater' (lines 766, 793) but two-group Fisher is two-sided (line 678); direction handled differently.
- RISK (rounded size-factor scaling): scaling counts then `round()` (lines 779-780, 1937-1940) can zero out small samples and is clamped with `max(0,...)` (lines 790-791, 1942-1943) — distorts Fisher inputs; "telomeric" vs "total" both fall through to the same code path (lines 747-753).
- SMELL (mostly top-level script): ~2200 lines of procedural module-level code with no `main()` / `if __name__`; not importable or testable as-is.
- SMELL (duplication): `cluster_to_enrichment` rebuilt 3x (lines 2279, 2596, 2722); `enrichment_colors`/`get_enrichment_color` rebuilt per plotting block (2314-2324, 2419-2423, 2472-2476, 2604-2608); `cluster_color_map` reassigned for UMAP (2613) shadowing the earlier global (2267).
- SMELL (broad bare excepts): `except:` swallows all errors in per-group Fisher (lines 1927, 1946).
- SMELL (pandas `groupby().apply(lambda)`): lines 869, 1012 — slow and deprecation-prone on newer pandas; also `read_feature_lengths` uses `.unstack` keyed on a `read` column that doesn't exist in `chrom_data` (it's `sequence`) — see next.
- POSSIBLE BUG (structure mode): line 355 groups by `['read','feature']` but `load_bed_file` names the column `sequence`, not `read` (line 492). This would raise `KeyError` unless every chromosome has <2 reads. ASSESSMENT: needs runtime verification — structure mode may be effectively unexercised/broken for multi-read chromosomes.
- RISK (knee depends on max-k): observed-range normalization makes selected_k sensitive to `--max-k` and early-stopping cutoff; acknowledged in comment (lines 1700, 1705) but still the default selection metric.
- SMELL (`level`/`n_total` sparsity, `mannwhitneyu` import at 1912 is unused after import).
- DEAD-ish: `node_to_cluster_leaf`/`leaf_counter`/`below_cut_threshold` partly unused scaffolding in above-cut remap (lines 2021-2026); `_palette`-style `read_names` duplicates `seq_names` in NPZ (lines 2210-2211) for back-compat.
- RISK (`sys.exit` inside helper): `run_structure_mode` exits the process (line 417), and several error paths `sys.exit(1)` (1463) — hard to wrap in a CLI/test harness.
- RISK (TeeLogger never closed; `random_state=42` fixed everywhere — reproducible but SVD/UMAP/silhouette all hard-coded).

## 10. Testability notes
ASSESSMENT — Prime pure-function unit targets (currently module-level, must be extracted first):
- `_build_matrix_from_features` / `build_combined_matrix` / `build_layer_matrix` and the layer helpers (`detect_feature_layers`, `split_feature_into_layers`, `get_layer_features`, `get_edges`, `get_weighted_edges`): deterministic NumPy in/out — assert matrix shapes, z-score, and block-reweight factor `sqrt(n_edge/n_abund)`.
- Knee detection (lines 1699-1726): feed synthetic composite curves, assert selected k at known elbows; test max-k sensitivity.
- Composite score formula (1541-1546): trivial arithmetic test.
- `calculate_enrichment_two_group` / `calculate_enrichment_per_sample` (645-824): build small contingency fixtures, assert odds/p and the size-factor scaling path; test the single-group edge case (651-659) and empty-sample case (728-738).
- `fast_enrichment_check` (1412-1450): assert odds-ratio + purity on hand counts.
- `update_enrichment_label` (2101-2108): table-driven.
- Levenshtein (lines 241-255) and normalized-distance matrix (304-312): classic unit tests.

ASSESSMENT — Fixtures: a tiny 2-sample BED with a handful of reads/features; a 2-row sample-metadata TSV; a total-reads TSV. Integration-only: SVD/clustering full pipeline, all matplotlib/plotly/UMAP outputs (UMAP requires `umap-learn`; assert files exist + NPZ keys). Golden tests: hash `cluster_analysis.tsv` / `sequence_assignments.tsv` under fixed `--n-clusters` and `random_state` seeds. The numerous side-effecting plot blocks should be isolated behind functions so the numeric core can be tested without rendering.

## 11. Open questions for the user
1. Should `--analysis-mode structure` become its own subcommand, and is structure mode actually used/working? (Suspected `read` vs `sequence` column bug at line 355 — see Section 9.)
2. Is `composite-knee` truly the desired default given its documented `--max-k` sensitivity, or should silhouette/Calinski be the default for reproducibility?
3. For per-sample enrichment, do you want FDR applied over the full sample x cluster p-value grid rather than only over per-cluster min-p (current behavior)?
4. Is the rounded size-factor normalization (telomeric/total) the intended compositional correction, or would CLR/proportion-test be preferred? Note "telomeric" and "total" currently share one code path.
5. Should the local FDR (scipy `false_discovery_control`, supports BH+BY) stay, or migrate to `karyoplot.mpl.statistics.apply_fdr` (BH-only)? They differ.
6. Which plotting belongs in shared `karyoplot.mpl` vs analysis-local? (Style/save/sig-stars vs clustering-specific dendrogram/knee figures.)
7. Default length window [10000, 50000] bp — keep as global default or make telomere-specific/per-subcommand?
8. Backward-compat: downstream `cluster_plot` reads the exact NPZ keys (incl. duplicate `read_names`); must the new writer preserve them byte-for-byte?
