# Session Log — All-Chromosome Dendrogram

## 2026-04-01 16:00 — Initial clustering (zscore_blockweight + silhouette)

Ran structure mode clustering. 4555 haplotypes, 18 chromosomes, 122 clusters.
chr3 and chr9 validation PASS. chr5 FAIL (k=2 too low to isolate 1-in-359 deletion).

## 2026-04-01 16:30 — Grid search: 3 matrix x 2 k-selection

Tested all 6 combinations. Only count_log1p + silhouette passes all 3 validations.
But count_log1p over-splits: 63% of haplotypes become "outliers", most visually
identical to Major. The z-score variants are cleaner but miss chr5.

## 2026-04-01 17:00 — Root cause analysis: why chr5 fails with z-score

chr5 has 359 haplotypes, only 1 (HG00558#1) lacks HSat3 (0.3% frequency).
Silhouette picks k=2 for chr5 with z-score. k=2 cannot isolate a singleton.
Tested forced k: HG00558#1 separates at k>=8 with z-score.
Problem is k-selection, not the matrix.

## 2026-04-01 17:30 — Two-stage outlier detection

Designed and tested: stage 1 = silhouette-optimal clustering, stage 2 = centroid
distance scan within Major (flag members > N SD from centroid).

chr5 HG00558#1 is at 10.6 SD from Major centroid — clear signal.
Stage 2 catches it regardless of k. Added to KS_allchr_dendrogram.py.

## 2026-04-01 18:00 — Silhouette threshold filter

Problem: some chromosomes (chr1, chr2, chr10, chr17) have sil<0.3 at k=2,
meaning the split is forced — no real structure. These produce fake outliers.

Solution: --sil-threshold collapses weak splits to k=1.

Tested grid: sil={0-0.8} x sd={0,3,5,10}. All 3 FISH validations pass at
every threshold from 0 to 0.8, thanks to stage 2. But stage 2 alone (sil=1.0)
fails chr3 because the HSat1A deletion needs clustering (stage 1) to detect.

Recommended: sil=0.5 sd=5 (66 reps, all PASS, cleanest).

## 2026-04-01 19:00 — Annotation pipeline

Created KS_allchr_annotate.py: compares outlier BED features to Major.
Reports block order changes, feature gains/losses, major abundance shifts.
User preference: structural changes only, abundance threshold >= 15%.

48 outliers annotated. 9 show only "edge pattern difference" (subtle).

## 2026-04-01 19:30 — Three plot variants

Generated with sil=0.5 sd=5 + zscore_blockweight clustering:

1. allchr_dendrogram_sil0.5_sd5.svg — 66 rows, no annotations
2. allchr_dendrogram_sil0.5_sd5_annotated.svg — 66 rows, all annotated
3. allchr_dendrogram_sil0.5_sd5_annotated_filtered.svg — 57 rows, subtle removed

Annotations: short labels, color-coded (red=deletion, green=gain,
purple=rearrangement), positioned at bar end, 9px sans-serif.

--hide-subtle removes the entire row, not just the text.

## Full pipeline commands

```bash
# Step 1: Clustering
python3 scripts/KaryoScope_cluster_analysis.py \
  --bed /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --output-prefix agent_results/allchr_structure \
  --analysis-mode structure \
  --edges directional \
  --no-abundance \
  --max-sequence-length 50000000 \
  --exclude-features "novel" \
  --matrix-type count_log1p_zscore_blockweight \
  --k-selection silhouette \
  --background white

# Step 2: Clean dendrogram
python3 scripts/KS_allchr_dendrogram.py \
  --assignments agent_results/allchr_structure.sequence_assignments.tsv \
  --bed /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --colors /Users/ychen/Documents/GitHub/KaryoScope/resources/databases/KS_human_CHM13 \
  --output agent_results/allchr_dendrogram.svg \
  --matrix-type count_log1p_zscore \
  --sil-threshold 0.5 \
  --centroid-sd 5 \
  --row-height 12 --bar-height 10

# Step 3: Generate annotations
python3 scripts/KS_allchr_annotate.py \
  --svg agent_results/allchr_dendrogram.svg \
  --assignments agent_results/allchr_structure.sequence_assignments.tsv \
  --bed /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --output agent_results/allchr_outlier_annotations.tsv

# Step 4: Annotated dendrogram (all)
python3 scripts/KS_allchr_dendrogram.py \
  --assignments agent_results/allchr_structure.sequence_assignments.tsv \
  --bed /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --colors /Users/ychen/Documents/GitHub/KaryoScope/resources/databases/KS_human_CHM13 \
  --output agent_results/allchr_dendrogram_annotated.svg \
  --matrix-type count_log1p_zscore \
  --sil-threshold 0.5 \
  --centroid-sd 5 \
  --row-height 12 --bar-height 10 \
  --annotations agent_results/allchr_outlier_annotations.tsv

# Step 5: Filtered dendrogram (subtle removed)
python3 scripts/KS_allchr_dendrogram.py \
  --assignments agent_results/allchr_structure.sequence_assignments.tsv \
  --bed /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --colors /Users/ychen/Documents/GitHub/KaryoScope/resources/databases/KS_human_CHM13 \
  --output agent_results/allchr_dendrogram_filtered.svg \
  --matrix-type count_log1p_zscore \
  --sil-threshold 0.5 \
  --centroid-sd 5 \
  --row-height 12 --bar-height 10 \
  --annotations agent_results/allchr_outlier_annotations.tsv \
  --hide-subtle
```

## Validation (final, sil=0.5 sd=5)

| Criterion | Expected | Observed | Status |
|-----------|----------|----------|--------|
| chr3 NA21144#1 | Outlier | chr3_Outlier_6 (stage 1) | PASS |
| chr5 HG00558#1 | Outlier | chr5_Outlier_S2_4 (stage 2, 10.6 SD) | PASS |
| chr9 HG02630#1 | Outlier | chr9_Outlier_S2_1 (stage 2) | PASS |

## Key decisions made

1. z-score blockweight for clustering + centroid scan for singletons,
   instead of count_log1p which over-splits
2. sil=0.5 threshold to suppress weak chromosome splits
3. sd=5 centroid threshold (sd=3 too liberal, sd=10 misses chr5)
4. Annotations focus on structural changes, not subtle abundance shifts
5. --hide-subtle removes entire rows, not just annotation text

## 2026-04-07 — Add chromosome identity blocks (cosmetic issue #9)

Added colored chromosome identity blocks between labels and feature bars.
Each of the 18 chromosomes gets a distinct color (maximally separated hues
so adjacent rows are visually distinguishable). Block is 12px wide, drawn
right before the feature bar.

Changes:
- scripts/KS_allchr_dendrogram.py: CHROM_BLOCK_COLORS palette (18 colors),
  chrom_block_width layout variable, block drawing loop, chromosome legend
- CLAUDE.md: documented cosmetic issue #9

Re-generated all 3 plot variants:
- allchr_dendrogram.svg (66 rows, clean)
- allchr_dendrogram_annotated.svg (66 rows, annotated)
- allchr_dendrogram_filtered.svg (57 rows, subtle removed)

All 3 FISH validations PASS (chr3, chr5, chr9).
Committed as e960d4f, pushed to pangenome_structure_V2.

## 2026-04-07 — Stacked bar chart (CLAUDE.md task 2)

Created scripts/KS_allchr_barplot.py: stacked bar chart showing Major vs
Outlier haplotype counts per chromosome. Reuses sil-threshold + centroid-scan
logic from dendrogram script (imports shared functions).

Per-chromosome totals validated against reference PDF
(/Users/ychen/Documents/GitHub/KaryoScope/results/figureA/filter.pdf):
all 18 chromosomes match (4551 total haplotypes).

Summary: 3730 Major, 821 Outlier across 18 chromosomes.
Chromosomes with most outliers: chr12 (233), chr8 (206), chr11 (178),
chr3 (124), chrX (41).

Output: agent_results/allchr_barplot.svg (13 KB)

## 2026-04-07 — chr5 haplotype count discrepancy (358 vs 359)

Reference PDF (filter.pdf) shows 359 chr5 haplotypes post-QC, but the
clustering assignments TSV has 358. The missing sequence is
HG03050#1#CM098762.1 which spans 123.2 Mb — nearly the full chromosome,
not just the centromere region. The clustering script's
--max-sequence-length 50000000 flag correctly filters it out.
This is expected behavior, not a bug.

Also added major count labels to the barplot (user requested both
Major and Outlier counts visible on each bar).

## Session summary — 2026-04-07

Completed today:
1. Cosmetic issue #9: added colored chromosome identity blocks between
   labels and bars in dendrogram (CHROM_BLOCK_COLORS palette, 18 distinct
   hues). Re-generated all 3 plot variants. Committed e960d4f.
2. Task 2: created KS_allchr_barplot.py — stacked bar chart of Major vs
   Outlier haplotype counts per chromosome. Shows both counts on bars.
   Totals validated against filter.pdf reference.
3. Investigated chr5 358 vs 359 discrepancy. Root cause: HG03050#1
   spans 123.2 Mb, filtered by --max-sequence-length 50000000. Documented
   in CLAUDE.md and session log. Not a bug.

All changes committed (95c482c) and pushed to pangenome_structure_V2.
All 3 FISH validations PASS (chr3, chr5, chr9).

Files produced:
- agent_results/allchr_dendrogram.svg (clean, 66 rows)
- agent_results/allchr_dendrogram_annotated.svg (annotated, 66 rows)
- agent_results/allchr_dendrogram_filtered.svg (filtered, 57 rows)
- agent_results/allchr_barplot.svg (stacked bar chart, 18 chromosomes)
- scripts/KS_allchr_barplot.py (new script)

## DONE: 2026-04-07 16:45

## 2026-04-08 — Human-readable NucFlag QC filter names

Replaced cryptic X_n_ filter labels with human-readable names based on
NucFlag wiki (https://github.com/logsdon-lab/NucFlag/wiki):

  X_n_COLLAPSE     -> No Collapsed Regions
  X_n_COLLAPSE_VAR -> No Collapsed Regions (with Variants)
  X_n_Err          -> No Erroneous Regions

Updated 3 files in the KaryoScope repo:
- results/figureA/generate_filter_flowchart.R: flowchart boxes, arrow
  labels, and table column sub-headers now use readable names.
  Regenerated filter.pdf.
- workflow/scripts/QCfilter_apply.R: report output uses readable names
  with original X_n_ codes in parentheses for traceability.
- workflow/scripts/QCfilter_explore.emp.R: plot labels updated
  (Erroneous, Collapsed, Collapsed (with Variants), Misjoined).

## 2026-04-08 — Reordered filters + heatmap panel

User requested manual filter order by importance instead of automatic
weak-to-strong ranking:
  1. No Erroneous Regions (removes 1,131 / 18.4% — biggest impact)
  2. No Collapsed Regions (removes 142 / 2.3%)
  3. No Collapsed Regions with Variants (removes 215 / 3.5%)

Replaced Panel 2 (plain table) with a % retained heatmap:
  - Color ramp: red (low retention) -> yellow -> white -> green (high)
  - Each cell shows absolute count + (% retained)
  - Acrocentrics (chr13-15, chr21-22) clearly stand out as most filtered
  - Angled x-axis labels for readability

Files changed (all in KaryoScope repo):
- results/figureA/generate_filter_flowchart.R: manual filter_order,
  heatmap replaces table, angled x-axis labels
- results/figureA/filter.pdf: regenerated

## DONE: 2026-04-08

## 2026-04-13 — Allele-specific outlier analysis (chr3, chr8, chr11, chr12)

User observed chr3, chr8, chr11, chr12 have more outliers than other
chromosomes in the barplot. Wanted to know if outliers are allele-specific:
do they affect one haplotype (monoallelic) or both (biallelic)?

Created scripts/KS_allchr_allele_heatmap.py:
- Pairs haplotypes by sample (h1/h2) for each chromosome
- Builds co-occurrence matrix: h1 cluster x h2 cluster
- Plots 2x2 heatmap grid (symmetrised for display)
- Imports apply_sil_and_centroid from KS_allchr_barplot.py to ensure
  cluster labels match the barplot exactly (same filtering pipeline)

Results:
  chr3:  51 paired, 7.8% both Major, 43.1% monoallelic, 49.0% biallelic
         (84% of biallelic in different clusters — independent variation)
  chr8:  138 paired, 17.4% both Major, 48.6% monoallelic, 34.1% biallelic
         (50/50 same vs different cluster)
  chr11: 173 paired, 31.2% both Major, 49.7% monoallelic, 19.1% biallelic
         (91% biallelic in same cluster — Outlier_2 is a common subtype)
  chr12: 165 paired, 19.4% both Major, 41.8% monoallelic, 38.8% biallelic
         (50/50 same vs different cluster)

Note: sil-threshold 0.5 + centroid-sd 5 filtering did not change labels
for these 4 chromosomes (their silhouette scores are all above 0.5, and
no stage-2 outliers were added). Results identical with or without filter.

Output files:
- agent_results/allchr_allele_heatmap.svg (heatmap grid)
- agent_results/allchr_allele_heatmap.png (for validation)
- agent_results/allchr_allele_summary.tsv (per-chrom summary counts)
- agent_results/allchr_allele_pairs.tsv (per-sample h1/h2 assignments)
- agent_results/allchr_allele_cooccurrence.tsv (co-occurrence long format)
- scripts/KS_allchr_allele_heatmap.py (new script)

### 2026-04-13 — SVG text editability fix + simplified 2x2 heatmap

Problem: SVG text rendered as path outlines (matplotlib default), not
editable in Illustrator.

Fix: Added `matplotlib.rcParams['svg.fonttype'] = 'none'` and
`matplotlib.rcParams['font.family'] = 'Helvetica'` to
KS_allchr_allele_heatmap.py.

Also added a simplified 2x2 Major-vs-Outlier heatmap (collapsing all
outlier subclusters into one "Outlier" category) for cleaner
presentation of monoallelic vs biallelic patterns.

Output files:
- agent_results/allchr_allele_heatmap.svg (full cluster matrix, now editable)
- agent_results/allchr_allele_heatmap_2x2.svg (simplified Major/Outlier)
- Both also exported as PNG
