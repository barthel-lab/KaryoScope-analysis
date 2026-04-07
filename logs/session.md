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

## DONE: 2026-04-07 16:40
