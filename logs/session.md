# Session Log — All-Chromosome Dendrogram

## 2026-04-01 19:30 — Clustering analysis (structure mode)

**Command:**
```
python3 scripts/KaryoScope_cluster_analysis.py \
  --bed pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
  --output-prefix agent_results/allchr_structure \
  --analysis-mode structure --edges directional --no-abundance \
  --max-sequence-length 500000000 --exclude-features "novel" \
  --matrix-type count_log1p_zscore_blockweight --k-selection silhouette
```

**Output:** 4555 haplotypes across 18 included chromosomes, 122 total clusters.

**Validation:**
- chr3 NA21144#1#CM094092.1 → `chr3_Outlier_6` — PASS
- chr5 HG00558#1#CM088494.1 → `chr5_Major` — **FAIL** (known issue, raw_div=19.35 but only 2 clusters found)
- chr9 HG02630#1#CM091811.1 → `chr9_Outlier_9` — PASS

**Decision:** Also tested `--matrix-type count_log1p_zscore --k-selection calinski` for chr5 — same result (2 clusters, HG00558#1 in Major). This is a pre-existing per-chromosome clustering limitation; chr5 only finds 2 clusters because the outlier's edge signature is not sufficiently distinct from major. Flagged and proceeding.

## 2026-04-01 19:35 — New script: KS_allchr_dendrogram.py

**Created:** `scripts/KS_allchr_dendrogram.py`

Features:
- Loads structure mode assignments + raw BED
- Excludes chr13/14/15/21/22/chrY
- Selects 1 representative per cluster (or `--all-haplotypes` for all)
- Builds GLOBAL feature matrix (directional edges, log1p+zscore) across ALL chromosomes
- Computes GLOBAL pairwise distance and Ward linkage (not reusing per-chrom distances)
- Draws: dendrogram | chromosome labels | feature bars | outlier labels
- Per-row chromosome coloring (handles interleaved rows from dendrogram ordering)
- Contiguous-run chromosome labels (only labels blocks >= 2 rows)
- Scale bar, legend, Major/Outlier dot indicators

## 2026-04-01 19:39 — Generated plots

| File | Description | Size |
|------|-------------|------|
| `agent_results/allchr_dendrogram.svg` | 122 cluster reps, 12px rows | 267 KB |
| `agent_results/allchr_dendrogram_allhaps.svg` | 4555 haplotypes, 8px rows | 8.3 MB |

Canvas dimensions:
- Representatives: 1070 x 1744 px
- All haplotypes: 1070 x 36720 px

Cosmetic fixes applied (vs screenshot1.png):
1. Bar height: 10px (representative) / 7px (all-haps) — meets 8px target
2. Zero gap between dendrogram and bar panel — dendrogram tips flush with bars
3. Global distance recalculation — not reusing per-chromosome distances

## Validation summary

| Criterion | Expected | Observed | Status |
|-----------|----------|----------|--------|
| chr3 NA21144#1 | Outlier cluster | chr3_Outlier_6 | PASS |
| chr5 HG00558#1 | Distinct outlier | chr5_Major (raw_div=19.35) | FAIL (known) |
| chr9 HG02630#1 | Outlier cluster | chr9_Outlier_9 | PASS |

## 2026-04-01 19:45 — chr5 fix: switch to count_log1p

**Problem:** `count_log1p_zscore_blockweight` and `count_log1p_zscore` both find only 2 clusters for chr5, placing HG00558#1 (HSat3 deletion) in Major. Z-score normalization flattens the magnitude differences that distinguish the deletion.

**Investigation:** HG00558#1 has 0% hsat3 vs normal ~5-7%. With edges-only + z-score, the missing `hor->hsat3` transition is one column among hundreds — too subtle.

**Solution:** `--matrix-type count_log1p` (no z-score) preserves raw magnitude differences. Result: chr5 now has 9 clusters, HG00558#1 → `chr5_Outlier_5`. Trade-off: more clusters overall (some over-splitting), but all 3 validation criteria now pass.

**Final validation (count_log1p + silhouette):**

| Criterion | Expected | Observed | Status |
|-----------|----------|----------|--------|
| chr3 NA21144#1 | Outlier | chr3_Outlier_8 | PASS |
| chr5 HG00558#1 | Outlier | chr5_Outlier_5 | PASS |
| chr9 HG02630#1 | Outlier | chr9_Outlier_4 | PASS |

**Files produced:**
- `agent_results/allchr_structure_log1p.sequence_assignments.tsv` (512 KB)
- `agent_results/allchr_dendrogram_log1p.svg` (310 KB) — 139 cluster reps

## Recommended next steps

1. Consider adding sample-name labels to the representative plot
2. PNG export needs `rsvg-convert` or native cairo — currently SVG only
3. Could add interactive HTML version with hover tooltips
4. May want to tune k-range for count_log1p to reduce over-splitting on some chromosomes
