# Quick Start Guide

This guide walks through the full KaryoScope analysis pipeline on the Core-3 example dataset: clustering, sequence annotation, cluster annotation, and visualization.

## Prerequisites

- Python 3.11+
- Required packages: numpy, pandas, scipy, scikit-learn, matplotlib, pyranges, umap-learn
- **Basic Sans** font family (bundled in the `fonts/` directory — registered automatically)

## Installation

```bash
# Clone the repository
git clone https://github.com/barthel-lab/KaryoScope-analysis.git
cd KaryoScope-analysis

# Install dependencies
pip3 install -r requirements.txt
```

## Example Data

The `data/` directory contains feature-annotated telomeric reads from three cell lines:

| Sample | Group | Description |
|--------|-------|-------------|
| HeLa | Telomerase | Telomerase-positive cancer cell line |
| U2OS | ALT | Alternative Lengthening of Telomeres (ALT) cancer cell line |
| IMR90 | Primary | Primary fibroblast cell line |

**Raw BED files** (`data/raw_bed/`):
Each sample has 6 featureset BED files produced by the KaryoScope Snakemake pipeline:

- `{sample}.telogator.1.KS_human_CHM13.subtelomeric.smoothed.features.bed.gz` — Subtelomeric features (canonical/non-canonical telomere, TAR1, ITS)
- `{sample}.telogator.1.KS_human_CHM13.region.smoothed.features.bed.gz` — Region features (satellites, chromosome arms)
- `{sample}.telogator.1.KS_human_CHM13.acrocentric.smoothed.features.bed.gz` — Acrocentric short arm features
- `{sample}.telogator.1.KS_human_CHM13.chromosome.smoothed.features.bed.gz` — Chromosome assignments
- `{sample}.telogator.1.KS_human_CHM13.gene.smoothed.features.bed.gz` — Gene annotations
- `{sample}.telogator.1.KS_human_CHM13.repeat.smoothed.features.bed.gz` — Repeat element annotations

**Color files** (`data/colors/KS_human_CHM13/`):
Feature-to-color mapping files used for visualization.

## Pipeline Overview

The analysis consists of six steps:

1. **Merge Feature Sets**: Combine subtelomeric and region features into a single BED file per sample
2. **Create Sample Metadata**: Define sample groups and colors
3. **Run Clustering**: Hierarchical clustering with enrichment analysis
4. **Sequence Annotation**: Compute per-read feature metrics across all featuresets
5. **Cluster Annotation**: Aggregate per-read metrics to cluster-level labels
6. **Cluster Plot**: Generate publication-quality cluster visualizations

---

## Step 1: Merge Feature Sets

Merge subtelomeric and region features for each sample using `--telomere-satellite-merge`. This creates a combined featureset where telomeric features (canonical/non-canonical telomere, TAR1, ITS) take priority from the subtelomeric file, and remaining positions are filled with satellite/region features:

```bash
mkdir -p results/merged_bed

for sample in HeLa U2OS IMR90; do
  python3 scripts/KaryoScope_merge_beds.py \
    --bed data/raw_bed/${sample}.telogator.1.KS_human_CHM13.subtelomeric.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.region.smoothed.features.bed.gz \
    --output results/merged_bed/${sample}.region_subtelomere_flat.merged.bed.gz \
    --telomere-satellite-merge
done
```

**What this does:**

- Takes telomeric features (canonical_telomere, noncanonical_telomere, TAR1, ITS) from the subtelomeric file as priority
- Fills remaining positions with features from the region file (satellites, chromosome arms)
- Outputs a single merged BED file per sample

**Expected output:**
```
=== Merging HeLa ===
  Subtelomeric intervals: 217,989
  Satellite intervals: 1,810,375
  Common reads: 3,474
  Priority intervals: 118,788
  Merged intervals: 1,922,572
```

---

## Step 2: Create Sample Metadata

The sample metadata file (`results/samples.tsv`) defines sample groups and colors for visualization. It is included with the example data:

```
sample	group	color
HeLa	Telomerase	#F07167
U2OS	ALT	#FBBF24
IMR90	Primary	#60A5FA
```

To create this file with proper tab delimiters, run:

```bash
printf 'sample\tgroup\tcolor\nHeLa\tTelomerase\t#F07167\nU2OS\tALT\t#FBBF24\nIMR90\tPrimary\t#60A5FA\n' > results/samples.tsv
```

---

## Step 3: Run Clustering Analysis

Run hierarchical clustering with block-weighted normalization and automatic k-selection:

```bash
python3 scripts/KaryoScope_cluster_analysis.py \
  --bed results/merged_bed/HeLa.region_subtelomere_flat.merged.bed.gz \
       results/merged_bed/U2OS.region_subtelomere_flat.merged.bed.gz \
       results/merged_bed/IMR90.region_subtelomere_flat.merged.bed.gz \
  --output-prefix results/Core3 \
  --sample-metadata results/samples.tsv \
  --comparison-mode per-sample \
  --matrix-type count_log1p_zscore_blockweight \
  --k-selection cosine-silhouette \
  --min-k 10 \
  --max-k 80 \
  --reduce-dims 50 \
  --min-sequence-length 3000 \
  --max-sequence-length 300000 \
  --edges symmetric \
  --background both \
  --umap
```

**Key parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--comparison-mode` | `per-sample` | Test each sample independently for enrichment |
| `--matrix-type` | `count_log1p_zscore_blockweight` | Log-normalize, z-score, and block-reweight features |
| `--k-selection` | `cosine-silhouette` | Select optimal k using cosine silhouette score |
| `--min-k` / `--max-k` | `10` / `80` | Range of cluster counts to evaluate |
| `--reduce-dims` | `50` | Reduce feature matrix to 50 SVD components |
| `--min-sequence-length` | `3000` | Minimum read length in bp |
| `--max-sequence-length` | `300000` | Maximum read length in bp |
| `--edges` | `symmetric` | Encode feature transitions symmetrically |
| `--background` | `both` | Generate both white and dark background plots |
| `--umap` | — | Generate UMAP projection |

**Expected runtime:** ~5-10 minutes for ~9,400 sequences

---

## Step 4: Sequence Annotation

Compute per-read feature metrics from all featureset BED files. This step analyzes feature densities, interspersion, and coverage across all featuresets for each read:

```bash
for sample in HeLa U2OS IMR90; do
  python3 scripts/KaryoScope_sequence_annotate.py \
    --bed data/raw_bed/${sample}.telogator.1.KS_human_CHM13.acrocentric.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.chromosome.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.gene.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.region.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.repeat.smoothed.features.bed.gz \
         data/raw_bed/${sample}.telogator.1.KS_human_CHM13.subtelomeric.smoothed.features.bed.gz \
         results/merged_bed/${sample}.region_subtelomere_flat.merged.bed.gz \
    --output results/${sample}.sequence_annotations.tsv.gz
done
```

**What this does:**

- Reads all 7 BED files per sample (6 raw featuresets + 1 merged)
- Computes per-read feature fractions, window densities (1 kb sliding window), and interspersion metrics
- Outputs a wide-format TSV with one row per read

!!! note "Pre-computed annotations available"
    Pre-computed sequence annotations are included in `data/sequence_annotations/` for participants who want to skip this step. To use them, copy the files:
    ```bash
    cp data/sequence_annotations/*.tsv.gz results/
    cp data/sequence_annotations/*.adaptive_thresholds.tsv results/
    ```

---

## Step 5: Cluster Annotation

First, combine per-sample sequence annotations into a single file, then annotate clusters with dominant features and auto-generated labels:

```bash
# Combine per-sample annotations
python3 scripts/KaryoScope_cluster_annotate.py \
  --prefix results/Core3 \
  --sequence-annotations results/Core3.combined_annotations.tsv \
  --output results/Core3.cluster_annotations.tsv \
  --featuresets region_subtelomere_flat,region,subtelomeric,chromosome,acrocentric,repeat,gene \
  --auto-label \
  --select-representatives 3 \
  --alt-samples U2OS
```

**Key parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--prefix` | `results/Core3` | Path prefix matching the clustering output |
| `--featuresets` | 7 featuresets | Annotation levels to analyze |
| `--auto-label` | — | Auto-assign cluster names from enrichment patterns |
| `--select-representatives` | `3` | Select 3 representative reads per cluster |
| `--alt-samples` | `U2OS` | Flag U2OS as the ALT sample for ALT-specific labeling |

---

## Step 6: Cluster Plot

Generate a publication-quality cluster visualization showing representative reads per cluster organized by the hierarchical clustering dendrogram:

```bash
python3 scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core3 \
  --bed results/merged_bed/HeLa.region_subtelomere_flat.merged.bed.gz \
       results/merged_bed/U2OS.region_subtelomere_flat.merged.bed.gz \
       results/merged_bed/IMR90.region_subtelomere_flat.merged.bed.gz \
  --colors data/colors/KS_human_CHM13 \
  --database KS_human_CHM13 \
  --output results/Core3.cluster_plot.svg \
  --featuresets region_subtelomere_flat \
  --background both \
  --vertical \
  --show-matrix \
  --orient-telomere-top \
  --n-per-cluster 1 \
  --png \
  --feature-mode transition \
  --min-feature-width 0.5 \
  --enrichment-grid \
  --label-column cluster_name
```

**Key parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--database` | `KS_human_CHM13` | Database name for color file lookup |
| `--vertical` | — | Vertical layout (reads as columns) |
| `--orient-telomere-top` | — | Place telomere end at top of each read |
| `--n-per-cluster` | `1` | Show 1 representative read per cluster |
| `--feature-mode` | `transition` | Color transitions between features |
| `--enrichment-grid` | — | Show enrichment grid alongside dendrogram |
| `--label-column` | `cluster_name` | Use auto-generated cluster names |
| `--show-matrix` | — | Show feature composition matrix |

---

## Output Files

### Clustering (Step 3)

| File | Description |
|------|-------------|
| `Core3.cluster_analysis.tsv` | Cluster statistics, enrichment p-values, sample composition |
| `Core3.sequence_assignments.tsv` | Per-sequence cluster assignments |
| `Core3.feature_matrix.npz` | Feature matrix with SVD components and linkage data |
| `Core3.cluster_analysis.pdf` | Dendrogram and composition bar charts |
| `Core3.k_selection.pdf` | k-optimization diagnostics |
| `Core3.svd_scree.pdf` | SVD variance explained plot |
| `Core3.umap.pdf` | UMAP projection colored by sample/cluster |
| `Core3.enrichment_bubble.pdf` | Enrichment significance bubble plot |
| `Core3.sample_percentage.pdf` | Sample composition heatmap |
| `Core3.log` | Full parameter log and runtime details |

### Sequence Annotation (Step 4)

| File | Description |
|------|-------------|
| `{sample}.sequence_annotations.tsv.gz` | Per-read feature metrics (one row per read) |
| `{sample}.sequence_annotations.adaptive_thresholds.tsv` | Per-feature adaptive thresholds |

### Cluster Annotation (Step 5)

| File | Description |
|------|-------------|
| `Core3.cluster_annotations.tsv` | Cluster-level feature summaries with auto-labels |

### Cluster Plot (Step 6)

| File | Description |
|------|-------------|
| `Core3.cluster_plot.svg` | Vector cluster visualization |
| `Core3.cluster_plot.png` | Rasterized cluster visualization (white background) |
| `Core3.cluster_plot_dark.png` | Rasterized cluster visualization (dark background) |

---

## Interpreting Results

### Cluster Summary

The analysis produces ~52 valid clusters (from k=61 selected by cosine silhouette) with enrichment categories:

```
Summary
============================================================
Total sequences: 9,367
Number of clusters: 52
  - HeLa-enriched: 9
  - U2OS-enriched: 21
  - IMR90-enriched: 13
  - mixed: 9
```

### Reading cluster_analysis.tsv

Key columns:

| Column | Description |
|--------|-------------|
| `Cluster` | Cluster ID |
| `Size` | Number of sequences in cluster |
| `{Sample}%` | Percentage of cluster from each sample |
| `P-value` | Fisher's exact test p-value |
| `Q-value` | FDR-corrected q-value |
| `Enrichment` | Enrichment call (e.g., "U2OS-enriched", "mixed") |
| `Centroid` | Sample of the centroid sequence |

### Enrichment Categories

- **Sample-enriched** (q < 0.05): Cluster significantly over-represented by one sample
- **mixed**: No significant enrichment after FDR correction

### Understanding the Cluster Plot

The cluster plot shows:

- **Dendrogram** (left): Hierarchical relationship between clusters
- **Feature bars** (center): Representative reads colored by genomic feature, oriented with telomere at top
- **Enrichment grid** (right): Sample enrichment significance per cluster
- **Cluster labels**: Auto-generated names based on dominant features (e.g., chromosome arm, satellite type)

---

## Next Steps

1. **Explore the output figures**: Review PDFs for detailed cluster structure and UMAP projections
2. **Compare features across samples**: Use `KaryoScope_enrichment_bubbles.py` for detailed enrichment visualization
3. **Investigate specific clusters**: Use `KaryoScope_plot_reads.py` to visualize all reads from clusters of interest
4. **Compare clustering runs**: Use `KaryoScope_compare_clusterings.py` to compare different parameter choices
