# KaryoScope Cluster Plot: Visualization Guide

## Overview

`KaryoScope_cluster_plot.py` is the visualization engine of the KaryoScope analysis pipeline. It generates publication-quality SVG images showing clustered telomeric reads as colored feature bars organized by hierarchical clustering. Each read is rendered as one or more horizontal bars where colors represent genomic features (telomeric repeats, satellite regions, chromosome arms, etc.), and reads are grouped by their cluster assignment with optional dendrogram visualization, enrichment annotations, and sample metadata.

**Key capabilities:**
- Horizontal and vertical layout modes
- Cluster-level and read-level dendrograms
- Per-sample enrichment visualization (bubbles or grid)
- Multi-featureset display (side-by-side or stacked)
- Density and combined tracks for fiberseq data
- Custom labels, filtering, and publication styling options
- SVG and PNG output with white/black/both backgrounds

---

## Table of Contents

**1. [Prerequisites and Workflow Context](#1-prerequisites-and-workflow-context)**

   - 1.1 [Pipeline Position](#11-pipeline-position)
   - 1.2 [Required Input Files](#12-required-input-files)
   - 1.3 [Required Resources](#13-required-resources)
   - 1.4 [BED Data Sources](#14-bed-data-sources)

**2. [Quick Start Examples (Core-4 Dataset)](#2-quick-start-examples-core-4-dataset)**

   - 2.1 [Basic Horizontal Plot](#21-basic-horizontal-plot)
   - 2.2 [Vertical Plot with Column Tracks](#22-vertical-plot-with-column-tracks)
   - 2.3 [Compact Single-Track Plot](#23-compact-single-track-plot)
   - 2.4 [With Pre-Selected Representatives](#24-with-pre-selected-representatives)

**3. [Layout Modes](#3-layout-modes)**

   - 3.1 [Horizontal Mode (Default)](#31-horizontal-mode-default)
   - 3.2 [Vertical Mode](#32-vertical-mode)
   - 3.3 [Column Tracks Mode](#33-column-tracks-mode)
   - 3.4 [Sample Matrix Mode](#34-sample-matrix-mode)

**4. [Complete Parameter Reference](#4-complete-parameter-reference)**

   - 4.1 [Required Parameters](#41-required-parameters)
   - 4.2 [Data Source Parameters](#42-data-source-parameters)
   - 4.3 [Display and Layout Parameters](#43-display-and-layout-parameters)
   - 4.4 [Dendrogram Parameters](#44-dendrogram-parameters)
   - 4.5 [Read Selection and Filtering Parameters](#45-read-selection-and-filtering-parameters)
   - 4.6 [Labels and Annotation Parameters](#46-labels-and-annotation-parameters)
   - 4.7 [Feature Rendering Parameters](#47-feature-rendering-parameters)
   - 4.8 [Enrichment Visualization Parameters](#48-enrichment-visualization-parameters)
   - 4.9 [Output Parameters](#49-output-parameters)

**5. [Feature Rendering](#5-feature-rendering)**

   - 5.1 [Transition Mode (Default)](#51-transition-mode-default)
   - 5.2 [Smooth Mode](#52-smooth-mode)
   - 5.3 [Color System](#53-color-system)

**6. [Dendrogram System](#6-dendrogram-system)**

   - 6.1 [Cluster-Level Dendrogram](#61-cluster-level-dendrogram)
   - 6.2 [Full Read-Level Dendrogram](#62-full-read-level-dendrogram)
   - 6.3 [Cutting and Grouping](#63-cutting-and-grouping)
   - 6.4 [Reordering](#64-reordering)

**7. [Enrichment Visualization](#7-enrichment-visualization)**

   - 7.1 [Single Bubble Mode (Default)](#71-single-bubble-mode-default)
   - 7.2 [Enrichment Grid Mode](#72-enrichment-grid-mode)

**8. [Advanced Examples (Core-4 Dataset)](#8-advanced-examples-core-4-dataset)**

   - 8.1 [Publication-Ready Vertical Plot](#81-publication-ready-vertical-plot)
   - 8.2 [Significant Clusters Only](#82-significant-clusters-only)
   - 8.3 [Sample-Specific Enrichment View](#83-sample-specific-enrichment-view)
   - 8.4 [Dendrogram with Cut Groups](#84-dendrogram-with-cut-groups)
   - 8.5 [Vertical Mode with Sample Matrix](#85-vertical-mode-with-sample-matrix)
   - 8.6 [Full Read-Level Dendrogram](#86-full-read-level-dendrogram)

**9. [Density and Combined Tracks](#9-density-and-combined-tracks)**

   - 9.1 [Density Featuresets](#91-density-featuresets)
   - 9.2 [Density Line Plots](#92-density-line-plots)
   - 9.3 [Rect Plots](#93-rect-plots)
   - 9.4 [Fiberseq Auto-Discovery](#94-fiberseq-auto-discovery)

**10. [Structural Mode](#10-structural-mode)**

**11. [Output Files](#11-output-files)**

**12. [Troubleshooting](#12-troubleshooting)**

**13. [Reference Tables](#13-reference-tables)**

   - 13.1 [Featureset Display Names](#131-featureset-display-names)
   - 13.2 [Complete Default Parameter Values](#132-complete-default-parameter-values)

---

## 1. Prerequisites and Workflow Context

### 1.1 Pipeline Position

`KaryoScope_cluster_plot.py` is the third step of the analysis pipeline. It consumes outputs from `KaryoScope_cluster_analysis.py` and optionally from `KaryoScope_select_representatives.py`:

```
Step 1: KaryoScope_cluster_analysis.py
  └── Produces: .cluster_analysis.tsv, .sequence_assignments.tsv,
                .feature_matrix.npz, .sample_metadata.tsv

Step 2 (optional): KaryoScope_select_representatives.py
  └── Produces: .representative_reads.tsv, .representative_reads.reads.txt

Step 3: KaryoScope_cluster_plot.py
  └── Produces: .svg (and optionally .png, .log)
```

### 1.2 Required Input Files

The script auto-discovers input files from the `--cluster-analysis-prefix`. Given a prefix like `results/Core4_telomere_region`, it looks for:

| File | Extension | Description |
|------|-----------|-------------|
| Sequence assignments | `.sequence_assignments.tsv` | Per-read cluster assignments with sample, rank, centroid distance (required) |
| Cluster analysis | `.cluster_analysis.tsv` | Cluster statistics, enrichment p-values, sample composition |
| Feature matrix | `.feature_matrix.npz` | Adjacency matrix, cluster centroids, cluster-level linkage for dendrogram computation |
| Sample metadata | `.sample_metadata.tsv` | Sample names, groups, and colors |

All files are auto-discovered by appending the extension to the prefix. The sequence assignments file is required; others are used when present.

### 1.3 Required Resources

The `--colors` directory must contain color definition files for each featureset being plotted. These files map feature names to hex colors:

```
# Format: {database}.{featureset}.colors.txt
# Tab-separated: feature_name\thex_color
canonical_telomere	#FF6B6B
noncanonical_telomere	#FF9999
TAR1	#4ECDC4
ITS	#C7F464
```

For the `KS_human_CHM13` database, available color files include:

| Featureset | # Colors | Contents |
|------------|----------|----------|
| `chromosome` | 31 | Chromosome-specific features (chr1-22, X, Y) |
| `subtelomeric` | 6 | Telomeric features (canonical/noncanonical telomere, TAR1, ITS, etc.) |
| `region` | 21 | Satellite regions, chromosome arms (p_arm, q_arm, satellites) |
| `repeat` | ~20 | Interspersed repeat families (Alu, LINE, SINE, etc.) |
| `acrocentric` | varies | Acrocentric chromosome features |

**Note:** Feature names with a `_specific` suffix in the BED data are automatically stripped to their base name for color lookup (e.g., `chr1_specific` matches the `chr1` color entry).

### 1.4 BED Data Sources

There are two ways to provide the underlying BED feature data:

**Option A: Explicit paths (`--bed`)**

```bash
--bed results/BJ.telogator.1.KS_human_CHM13.region.smoothed.features.bed.gz \
     results/HeLa.telogator.1.KS_human_CHM13.region.smoothed.features.bed.gz
```

**Option B: Auto-discovery (`--input-bed-prefix` + `--database`)**

```bash
--input-bed-prefix results --database KS_human_CHM13
```

Auto-discovery expects the directory structure:
```
{prefix}/{sample}/telogator/1/KaryoScope/{database}/
  └── {sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed.gz
```

Sample names are read from the `.sample_metadata.tsv` file.

---

## 2. Quick Start Examples (Core-4 Dataset)

All examples use the Core-4 dataset with four cell lines:

| Sample | Group | Description | Color |
|--------|-------|-------------|-------|
| BJ | Primary | Primary fibroblast cell line | `#40D392` |
| IMR90 | Primary | Primary fibroblast cell line | `#60A5FA` |
| HeLa | Telomerase | Telomerase-positive cancer cell line | `#F07167` |
| U2OS | ALT | Alternative Lengthening of Telomeres cancer cell line | `#FBBF24` |

### 2.1 Basic Horizontal Plot

The simplest invocation showing all clusters with up to 5 reads each, dendrogram visible, on a white background:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets chromosome,subtelomeric,region \
  --background white \
  --show-dendrogram \
  --n-per-cluster 5 \
  --output results/Core4_telomere_region.cluster_plot.svg
```

**Expected output:** An SVG image approximately 10,874 x 1,168 pixels showing ~51 clusters with 254 total reads. The dendrogram is drawn at the top, reads appear as vertical columns with stacked feature bars (chromosome, subtelomeric, region), and cluster brackets with enrichment labels appear below.

### 2.2 Vertical Plot with Column Tracks

Vertical orientation with each featureset in its own column, 3 reads per cluster:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --vertical \
  --column-tracks \
  --show-dendrogram \
  --n-per-cluster 3 \
  --output results/Core4.vertical.cluster_plot.svg
```

In vertical mode, reads are stacked top-to-bottom as rows with feature bars extending to the right. The dendrogram appears on the left.

### 2.3 Compact Single-Track Plot

One representative per cluster with telomere orientation and enrichment grid:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --vertical \
  --column-tracks \
  --n-per-cluster 1 \
  --enrichment-grid \
  --orient-telomere-top \
  --show-dendrogram \
  --output results/Core4.compact.cluster_plot.svg
```

This produces a compact overview where each cluster is represented by a single read, all oriented with telomeric features at the top. The enrichment grid shows per-sample bubbles next to each cluster.

### 2.4 With Pre-Selected Representatives

For curated figures, use `KaryoScope_select_representatives.py` to choose the best reads, then plot:

```bash
# Step 1: Select representative reads
python scripts/KaryoScope_select_representatives.py \
  --cluster-analysis results/Core4_telomere_region.cluster_analysis.tsv \
  --read-assignments results/Core4_telomere_region.sequence_assignments.tsv \
  --bed-prefix results \
  --n-per-cluster 5 \
  --preferred-min-length 20000 \
  --preferred-max-length 30000 \
  --output results/Core4_telomere_region.reps.tsv

# Step 2: Plot with pre-selected reads
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets chromosome,subtelomeric,region \
  --reads-file results/Core4_telomere_region.reps.reads.txt \
  --show-dendrogram \
  --output results/Core4.with_reps.cluster_plot.svg
```

The representative selection script scores reads based on feature matching (0.5 weight), read length (0.4 weight), and centroid proximity (0.1 weight), preferring reads in the 20-30 kb range.

---

## 3. Layout Modes

### 3.1 Horizontal Mode (Default)

In horizontal mode, reads are arranged left-to-right as vertical columns. This is the default orientation.

```
┌─────────────────────────────────────────────────┐
│                    Legends                       │
│                                                  │
│   ┌──── Dendrogram (optional) ─────┐            │
│   │        ┌─┐  ┌─┐                │            │
│   │     ┌──┤ ├──┤ ├──┐             │            │
│   │  ┌──┤  └─┘  └─┘  ├──┐         │            │
│   └──┘                   └──       │            │
│                                                  │
│  [Enrichment bubbles/grid]                       │
│  [Cluster brackets + labels]                     │
│  [Sample color bars]                             │
│                                                  │
│  ┌──┬──┬──┬──┬──┐  ┌──┬──┬──┐  ┌──┬──┐         │
│  │  │  │  │  │  │  │  │  │  │  │  │  │  ...     │
│  │FeatureBars    │  │FeatureBars│  │  │          │
│  │(stacked)      │  │(stacked)  │  │  │          │
│  └──┴──┴──┴──┴──┘  └──┴──┴──┘  └──┴──┘         │
│    Cluster 1         Cluster 2    Cluster 3      │
│                                                  │
│  [Featureset labels]                             │
│  [Scale bar]                                     │
│  [Color legend]                                  │
└─────────────────────────────────────────────────┘
```

**Best for:** Few reads per cluster, wide displays, traditional dendrogram views.

### 3.2 Vertical Mode

Enabled with `--vertical`. Reads are stacked top-to-bottom as rows. The dendrogram appears on the left.

```
┌────────────────────────────────────────────────────┐
│ Dendrogram │ Enrichment │ Feature Bars      │Labels│
│  (left)    │  bubbles   │ (rows = reads)    │      │
│            │            │                    │      │
│    ┐       │    ●       │ ▓▓▓▓░░▓▓▓▓▓░░▓▓  │ C1   │
│    ├──     │    ●       │ ▓▓░░░▓▓▓▓░░▓▓▓▓  │      │
│    │       │            │                    │      │
│    ├──     │    ●       │ ░░▓▓▓▓▓▓░░▓▓▓▓▓  │ C2   │
│    │  ┐    │    ●       │ ░░▓▓▓▓░░░▓▓▓▓▓▓  │      │
│    │  ├──  │    ●       │ ░▓▓▓░░▓▓▓▓▓░░▓▓  │      │
│    │  │    │            │                    │      │
│    ├──┘    │    ●       │ ▓▓▓▓▓▓▓░░░░░▓▓▓  │ C3   │
│    │       │            │                    │      │
│   ...      │   ...      │ ...                │ ...  │
│            │            │                    │      │
│            │            │ [Featureset labels] │      │
│            │            │ [Scale bar]         │      │
└────────────────────────────────────────────────────┘
```

**Best for:** Many clusters, publication figures, compact display with enrichment annotations.

### 3.3 Column Tracks Mode

Enabled with `--column-tracks`. Works with both horizontal and vertical orientations. Instead of stacking featuresets within each read, each featureset gets its own spatial region. In vertical mode, each featureset becomes a separate column. In horizontal mode, each featureset becomes a separate row.

```bash
--vertical --column-tracks --featuresets subtelomeric,region
```

Produces a layout where subtelomeric and region features are plotted side-by-side, making it easy to compare annotations across featuresets for the same read.

### 3.4 Sample Matrix Mode

Enabled with `--show-matrix` (requires `--vertical`). Adds a sample-by-cluster read count matrix between the dendrogram/bubbles and the feature bars. Includes a sample dendrogram, column bar plot (sample totals), and row bar plot (cluster totals).

```bash
--vertical --show-matrix --show-dendrogram
```

The matrix is drawn as a heatmap where cell intensity represents the number of reads from each sample in each cluster.

---

## 4. Complete Parameter Reference

### 4.1 Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--cluster-analysis-prefix` | string | Prefix from `cluster_analysis.py` outputs. Auto-discovers `.sequence_assignments.tsv`, `.feature_matrix.npz`, `.sample_metadata.tsv`, `.cluster_analysis.tsv` |
| `--colors` | path | Full path to colors database directory containing `{database}.{featureset}.colors.txt` files |
| `--output` | path | Output SVG file path |

### 4.2 Data Source Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--bed` | paths (multiple) | None | Full paths to BED files. Alternative to auto-discovery |
| `--input-bed-prefix` | path | None | Base directory for auto-discovery of BED files. Expects structure: `{prefix}/{sample}/telogator/1/KaryoScope/{database}/` |
| `--database` | string | auto-detected | Database name (e.g., `KS_human_CHM13`). Auto-detected from `--bed` paths if not provided; required with `--input-bed-prefix` |
| `--featuresets` | string | `chromosome,subtelomeric,region` | Comma-separated list of feature sets to plot |
| `--custom-beds` | NAME:PATH (multiple) | None | Custom BED files as additional feature tracks. Format: `featureset_name:/path/to/file.bed`. Requires matching color file in `--colors` directory |
| `--fiberseq` | path | None | Directory containing fiberseq BED files. Auto-discovers `*.FIRE.bed`, `*.LINKER.bed`, `*.m6A.bed`, `*.5mC.bed` and configures combined tracks |
| `--density-featuresets` | string | None | Comma-separated featuresets to render as density tracks with binning (e.g., `fiberseq_m6A,fiberseq_5mC`) |
| `--density-bin-size` | int | 300 | Bin size in bp for density computation |
| `--density-line-plot` | string | None | Combine featuresets as overlaid density line plots. Format: `fiberseq_m6A:fiberseq_5mC`. Uses `--density-bin-size` for binning |
| `--rect-plot` | string | None | Combine featuresets as stacked rectangles showing exact feature regions. Format: `fiberseq_FIRE:fiberseq_LINKER` |
| `--smoothness` | string | `smoothed` | Smoothness level for BED file auto-discovery (matches filename pattern) |

### 4.3 Display and Layout Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--background` | choice | `black` | Background color: `white`, `black`, or `both` (generates two SVGs) |
| `--png` | flag | False | Also export PNG alongside SVG (requires `rsvg-convert` at 4x resolution) |
| `--bar-width` | int | 8 | Width of each feature bar in pixels |
| `--bar-spacing` | int | 0 | Spacing in pixels between bars within a read group |
| `--read-spacing` | int | 12 | Spacing in pixels between read groups |
| `--cluster-spacing` | int | 30 | Spacing in pixels between clusters |
| `--ratio` | float | 1/300 (0.00333) | Scaling factor: pixels per base pair. A 30 kb read at default ratio = 100 pixels |
| `--oversample` | int | 1 | Internal oversampling factor for feature rasterization in `smooth` mode. Higher values resolve smaller features without changing image size |
| `--target-width` | int | None | Target image width in pixels. Auto-calculates `--ratio` to fit. Overrides `--ratio` if set |
| `--target-height` | int | None | Target image height in pixels. Auto-calculates `--ratio` to fit. Overrides `--ratio` if set |
| `--vertical` | flag | False | Rotate plot 90 degrees (dendrogram on left, reads as horizontal rows) |
| `--column-tracks` | flag | False | Display featuresets as separate columns (vertical) or rows (horizontal) instead of stacked within each read |
| `--show-matrix` | flag | False | Show sample-by-cluster read count matrix. Only works with `--vertical` |
| `--font-family` | string | `sans-serif` | Font family for text labels |

### 4.4 Dendrogram Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--show-dendrogram` | flag | False | Show hierarchical clustering dendrogram |
| `--hide-dendrogram` | flag | False | Completely hide the dendrogram (sets height to 0). Cluster ordering from the dendrogram is still applied |
| `--hide-brackets` | flag | False | Hide cluster brackets and labels for a cleaner dendrogram view |
| `--no-reorder` | flag | False | Disable dendrogram reordering. Keeps reads grouped by original cluster order (enrichment tier, then p-value) |
| `--full-dendrogram` | flag | False | Show complete hierarchical tree down to individual reads instead of cluster-level dendrogram. Computes Ward linkage from the adjacency matrix |
| `--dendro-cut` | string | None | Cut dendrogram into groups with extra spacing. Use `n:K` for K groups (e.g., `n:5`) or a number for distance threshold (e.g., `32`) |
| `--dendro-cluster-gap` | int | 0 | Extra gap in pixels between cluster groups in full dendrogram mode |
| `--show-threshold` | flag | False | Visualize the structural distance threshold on the dendrogram |
| `--structural-threshold` | float | 0.25 | Threshold for structural outlier clustering |

### 4.5 Read Selection and Filtering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--reads-file` | path | None | File with read names to include (one per line). Only these reads are plotted. Use `KaryoScope_select_representatives.py` to generate this file |
| `--n-per-cluster` | int | None | Maximum number of sequences per cluster. When no `--reads-file` is provided, selects this many reads per cluster using centroid-based proportional sampling |
| `--max-qvalue` | float | None | Filter to clusters with q_value <= this threshold (e.g., `0.05` for significant clusters only) |
| `--filter-enrichment` | strings (multiple) | None | Show only clusters matching these enrichment labels (e.g., `U2OS-enriched mixed`) |
| `--curated-reps` | path | None | TSV file with curated representative selection. Must have `cluster_id` and `curated_rep_i` columns indicating which rank (1-based) to plot per cluster |
| `--priority-samples` | string | None | Comma-separated sample names to prioritize when selecting representatives (structural mode) |

### 4.6 Labels and Annotation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cluster-labels` | path | None | TSV or Excel file with custom cluster labels. Must have `cluster_id` column and a label column |
| `--label-column` | string | `curated_annotation` | Column name for custom labels in `--cluster-labels` file |
| `--show-read-indices` | flag | False | Show read index labels (1, 2, 3, ...) next to each read |
| `--hide-read-labels` | flag | False | Hide read name labels above each read bar |
| `--show-cluster-numbers` | flag | False | Show cluster numbers below each read instead of featureset names |
| `--show-clade-id` | flag | False | Show clade ID in structural plot labels (e.g., `[C20]`) |
| `--show-clade-count` | flag | False | Show count of reads in each clade (e.g., `[n=15]`) |
| `--orient-telomere-top` | flag | False | Reorient reads so telomeric features (canonical/noncanonical telomere) are always at the top of the visualization |

### 4.7 Feature Rendering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--feature-mode` | choice | `transition` | `transition`: direct bp-to-pixel scaling preserving all feature boundaries. `smooth`: windowed majority vote for overview visualization |
| `--min-feature-width` | float | 0.5 | Minimum pixel width per feature in transition mode. Ensures small features remain visible |
| `--min-width-exclude` | patterns | `novel *arm* ct*` | Glob patterns for feature names to exclude from `--min-feature-width` enforcement. These features (e.g., chromosome arms) are allowed to render at their natural width even if below the minimum. Use `--min-width-exclude` with no arguments to clear all exclusions |

### 4.8 Enrichment Visualization Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--enrichment-grid` | flag | False | Show enrichment as a grid of bubbles (one per sample per cluster) instead of a single bubble per cluster. Requires `per-sample` comparison mode in cluster analysis. See [Section 7.2](#72-enrichment-grid-mode) for details |

### 4.9 Output Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--log-file` | flag | True | Save console output to `{output}.log`. Disable with `--no-log-file` |
| `--png` | flag | False | Also export PNG alongside SVG. Requires `rsvg-convert` (part of librsvg). Renders at 4x resolution |

---

## 5. Feature Rendering

The cluster plot converts BED feature intervals (base pair coordinates) into colored pixel rectangles. Two rendering modes are available, selected via `--feature-mode`.

### 5.1 Transition Mode (Default)

`--feature-mode transition`

Direct bp-to-pixel scaling that preserves every feature boundary. Each feature interval is scaled by the `--ratio` parameter and placed contiguously. The algorithm:

1. **Scale:** Each feature's pixel width = `(bp_end - bp_start) * ratio`
2. **Enforce minimum width:** Features below `--min-feature-width` pixels are expanded to the minimum (except those matching `--min-width-exclude` patterns)
3. **Redistribute:** If the total width exceeds the bar length after minimum enforcement, large features shrink proportionally to absorb the excess
4. **Place contiguously:** Features are placed end-to-end with no gaps or overlaps
5. **Run-length encode:** Adjacent same-color features are merged into single rectangles for efficient SVG output

**Example:** A 30 kb read at `--ratio 1/300` produces a 100-pixel bar. An interspersed 150 bp feature would naturally be 0.5 pixels but is guaranteed at least `--min-feature-width` (0.5 px by default).

**When to use:** Default mode. Best for detailed visualization where feature boundaries matter.

### 5.2 Smooth Mode

`--feature-mode smooth`

Windowed majority vote that summarizes features at pixel resolution. The algorithm:

1. For each pixel, define a bp window of size `1/ratio`
2. Tally bp coverage per color within the window
3. The color with the most coverage wins (plurality vote)
4. Run-length encode contiguous same-color pixels

When `--oversample` is greater than 1, rasterization happens at the oversampled resolution internally, then coordinates are divided back. This lets features smaller than `1/ratio` bp survive the majority vote and appear as fractional-pixel SVG rectangles.

**When to use:** Overview visualization where individual feature boundaries are less important than overall composition.

### 5.3 Color System

Colors are loaded from `{database}.{featureset}.colors.txt` files in the `--colors` directory. Each file is tab-separated with `feature_name\thex_color` pairs.

**Key behaviors:**

- **`_specific` suffix stripping:** Features with a `_specific` suffix in BED data (e.g., `chr1_specific`) are stripped to their base name (`chr1`) for color lookup
- **Unmapped features:** Features not found in the color file render as gray (`#444444`). A warning is printed listing uncolored features
- **Feature opacity:** Each feature rectangle has a default opacity of 1.0, which can be reduced for density tracks

---

## 6. Dendrogram System

### 6.1 Cluster-Level Dendrogram

Enabled with `--show-dendrogram`. Draws a dendrogram connecting clusters based on their structural similarity.

**How it works:**

1. Cluster centroids (mean feature vectors) are loaded from the `.feature_matrix.npz` file, pre-computed by `cluster_analysis.py`
2. The cluster-level linkage matrix (Ward's method on centroid distances) is also pre-computed and stored in the NPZ file
3. If some clusters are filtered out (via `--max-qvalue` or `--filter-enrichment`), the full tree is pruned using Bio.Phylo to extract only the displayed clusters while preserving tree topology
4. Leaf order from the dendrogram determines the display order of clusters (left-to-right in horizontal mode, top-to-bottom in vertical mode)

The dendrogram includes:
- Branch lines (horizontal and vertical connectors)
- A scale axis showing "Cluster Distance" with tick marks
- Optional cluster brackets and enrichment labels below

### 6.2 Full Read-Level Dendrogram

Enabled with `--full-dendrogram`. Shows a complete hierarchical tree down to individual reads instead of the cluster-level summary.

**How it works:**

1. The adjacency matrix is loaded from `.feature_matrix.npz`
2. The matrix is subset to only the displayed reads
3. Pairwise distances are computed with `scipy.spatial.distance.pdist`
4. Ward linkage is computed with `scipy.cluster.hierarchy.linkage`
5. Leaf ordering is optimized with `optimal_leaf_ordering` for cleaner visualization

This mode produces much taller/wider plots since every read is a leaf node.

### 6.3 Cutting and Grouping

The `--dendro-cut` parameter adds visual separation between major dendrogram groups:

- **`n:K`** - Cut the dendrogram into exactly K groups using `fcluster(method='maxclust')`
  ```bash
  --dendro-cut n:5  # Cut into 5 groups
  ```

- **distance threshold** - Cut at a specific distance value using `fcluster(method='distance')`
  ```bash
  --dendro-cut 32  # Cut at distance 32
  ```

Groups are separated by `2x --cluster-spacing` pixels. The `--dendro-cluster-gap` parameter adds additional spacing between cluster groups in full dendrogram mode.

### 6.4 Reordering

By default, cluster order is determined by the dendrogram's leaf ordering, which uses `optimal_leaf_ordering` to minimize crossings.

`--no-reorder` disables this and keeps the original order from `cluster_analysis.py`, which sorts clusters by:
1. Enrichment tier (Tier 0: 100% enriched, Tier 1: 80%+, Tier 2: all others)
2. Then by p-value (ascending)

---

## 7. Enrichment Visualization

### 7.1 Single Bubble Mode (Default)

In vertical mode, a single bubble is drawn next to each cluster. The bubble encodes three statistical dimensions:

| Visual Property | Data Mapping |
|----------------|-------------|
| **Size** (radius) | Cluster size (number of reads). Linear scaling from `min_radius=2` to `max_radius=8` pixels |
| **Color** | `log2(odds_ratio)` mapped to a diverging colormap: blue (depleted, OR < 1) to white (neutral, OR = 1) to red (enriched, OR > 1). Capped at ±4 |
| **Opacity** | `-log10(q_value)` mapped to 0.0 (q=1, not significant) to 1.0 (q <= 1e-10). More significant clusters are more opaque |

### 7.2 Enrichment Grid Mode

Enabled with `--enrichment-grid`. Draws a grid of bubbles with one column per sample and one row per cluster. Requires `per-sample` comparison mode in the upstream cluster analysis.

| Visual Property | Data Mapping |
|----------------|-------------|
| **Size** (radius) | Sample percentage in the cluster. `pct / 100 * bubble_radius`, minimum 30% of max radius |
| **Color** | `log2(odds_ratio)` on a diverging colormap: blue (depleted) to white to red (enriched) |
| **Opacity** | `-log10(p_value)` mapped to 0.0 (not significant) to 1.0. Samples with < 5% contribution are hidden |

**Core-4 example:** With 4 samples, the grid shows 4 columns of bubbles next to each cluster, making it easy to see which samples drive each cluster.

```bash
--enrichment-grid
```

---

## 8. Advanced Examples (Core-4 Dataset)

### 8.1 Publication-Ready Vertical Plot

Custom labels, enrichment grid, both backgrounds, PNG export, telomere orientation:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --vertical \
  --column-tracks \
  --n-per-cluster 1 \
  --orient-telomere-top \
  --enrichment-grid \
  --hide-read-labels \
  --show-cluster-numbers \
  --hide-brackets \
  --show-dendrogram \
  --cluster-labels results/Core4_telomere_region.cluster_annotations.tsv \
  --background both \
  --png \
  --output results/Core4.publication.cluster_plot.svg
```

This generates four files:
- `Core4.publication_white.svg` and `Core4.publication_black.svg`
- `Core4.publication_white.png` and `Core4.publication_black.png`

### 8.2 Significant Clusters Only

Filter to clusters with FDR-corrected q-value <= 0.05:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets chromosome,subtelomeric,region \
  --max-qvalue 0.05 \
  --show-dendrogram \
  --n-per-cluster 3 \
  --output results/Core4.significant_only.cluster_plot.svg
```

Only clusters where the enrichment test reached statistical significance are displayed. The dendrogram is automatically re-computed for the subset of displayed clusters.

### 8.3 Sample-Specific Enrichment View

Show only clusters enriched for specific cell lines:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --filter-enrichment U2OS-enriched HeLa-enriched \
  --show-dendrogram \
  --n-per-cluster 5 \
  --output results/Core4.cancer_enriched.cluster_plot.svg
```

This filters to only U2OS-enriched and HeLa-enriched clusters, useful for focusing on cancer-specific telomeric structural patterns.

### 8.4 Dendrogram with Cut Groups

Cut the dendrogram into 5 major groups with visual separation:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets chromosome,subtelomeric,region \
  --show-dendrogram \
  --dendro-cut n:5 \
  --cluster-spacing 0 \
  --hide-brackets \
  --n-per-cluster 1 \
  --output results/Core4.dendro_cut.cluster_plot.svg
```

Groups are separated by `2x cluster-spacing` pixels. Combined with `--cluster-spacing 0` and `--hide-brackets`, this produces a clean visualization where only the major dendrogram groups are visually separated.

### 8.5 Vertical Mode with Sample Matrix

Add a sample-by-cluster read count matrix:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --vertical \
  --show-matrix \
  --show-dendrogram \
  --n-per-cluster 5 \
  --output results/Core4.matrix.cluster_plot.svg
```

The matrix shows a heatmap of how many reads from each sample (BJ, IMR90, HeLa, U2OS) fall into each cluster, with sample dendrogram, column bar plot, and row bar plot.

### 8.6 Full Read-Level Dendrogram

Show the complete hierarchical tree down to individual reads:

```bash
python scripts/KaryoScope_cluster_plot.py \
  --cluster-analysis-prefix results/Core4_telomere_region \
  --input-bed-prefix results \
  --database KS_human_CHM13 \
  --colors resources/KS_human_CHM13 \
  --featuresets subtelomeric,region \
  --vertical \
  --full-dendrogram \
  --n-per-cluster 5 \
  --output results/Core4.full_dendro.cluster_plot.svg
```

This computes Ward linkage from the adjacency matrix for all displayed reads (not just cluster centroids). The tree shows individual read-level branching, which is useful for inspecting within-cluster heterogeneity.

---

## 9. Density and Combined Tracks

These features allow visualization of fiberseq epigenetic data and other dense feature tracks alongside standard feature annotations.

### 9.1 Density Featuresets

Render featuresets as density heatmaps instead of individual feature rectangles:

```bash
--density-featuresets fiberseq_m6A,fiberseq_5mC \
--density-bin-size 300
```

Small features are binned into windows of `--density-bin-size` base pairs and colored by density level rather than individual feature identity.

### 9.2 Density Line Plots

Overlay multiple featuresets as line plots in a single track:

```bash
--density-line-plot fiberseq_m6A:fiberseq_5mC
```

Each featureset's density is computed at `--density-bin-size` resolution and drawn as a line. Lines are colored by the first color in each featureset's color file. The individual featuresets are removed from the main featureset list and replaced with a combined `density_line` track.

### 9.3 Rect Plots

Combine multiple featuresets as stacked rectangles showing exact feature positions:

```bash
--rect-plot fiberseq_FIRE:fiberseq_LINKER
```

Unlike density plots, rect plots show the exact feature regions as colored rectangles. Useful for features with well-defined boundaries like FIRE elements and linker regions.

### 9.4 Fiberseq Auto-Discovery

The `--fiberseq` shortcut auto-discovers and configures all fiberseq tracks:

```bash
--fiberseq /path/to/fiberseq_dir/
```

**What it does:**
1. Scans the directory for `*.FIRE.bed`, `*.LINKER.bed`, `*.m6A.bed`, `*.5mC.bed`
2. Creates a combined `FIRE_LINKER` file and adds it as a rect plot track
3. Adds m6A and 5mC as custom BED files
4. Auto-configures `--density-line-plot fiberseq_m6A:fiberseq_5mC` if both are found

This is equivalent to manually setting up:
```bash
--custom-beds fiberseq_FIRE_LINKER:/path/combined.bed \
              fiberseq_m6A:/path/m6A.bed \
              fiberseq_5mC:/path/5mC.bed \
--rect-plot fiberseq_FIRE:fiberseq_LINKER \
--density-line-plot fiberseq_m6A:fiberseq_5mC
```

---

## 10. Structural Mode

Structural mode is automatically activated when the `.feature_matrix.npz` file contains `mode=structure`. This mode is used for per-chromosome structural analysis rather than the standard cluster visualization.

**Key differences from standard mode:**

- Generates one SVG per chromosome plus an `all_chromosomes.svg` grid layout
- Reads are classified as "Major" (consensus) or "Outlier" (divergent) types
- Outliers are sorted by `raw_divergence` (descending)
- `--priority-samples` controls which sample's reads are preferred as representatives
- Major reads are always shown first, followed by outliers

**Output structure:**
```
{output_base}.{chromosome}.svg      # Per-chromosome SVG
{output_base}.all_chromosomes.svg   # Combined grid of all chromosomes
```

---

## 11. Output Files

| File | Condition | Description |
|------|-----------|-------------|
| `{output}.svg` | Always | Primary output. Scalable vector graphics |
| `{output}.png` | `--png` | Rasterized at 4x resolution via `rsvg-convert` |
| `{output}.log` | `--log-file` (default: True) | Console output including auto-discovery paths, loading statistics, image dimensions, full parameter table, and original command |
| `{base}_white.svg` | `--background both` | White background variant |
| `{base}_black.svg` | `--background both` | Black background variant |

**Image dimensions** depend on:
- Number of clusters displayed
- Reads per cluster (`--n-per-cluster`)
- Bar width (`--bar-width`), spacing (`--bar-spacing`, `--read-spacing`, `--cluster-spacing`)
- Scaling ratio (`--ratio`) and read lengths
- Whether dendrogram, enrichment grid, or matrix is shown

Typical sizes range from 1,000 to 15,000 pixels wide and 800 to 6,000 pixels tall.

---

## 12. Troubleshooting

### "Error: Read assignments file not found"

The `.sequence_assignments.tsv` file is required. Verify the `--cluster-analysis-prefix` is correct:
```bash
ls results/Core4_telomere_region.sequence_assignments.tsv
```

### "Error: Could not determine database from BED file paths"

When using `--bed`, the database name is parsed from the BED file path. Ensure paths follow the expected naming convention. Alternatively, provide `--database` explicitly.

### "Warning: features not in colors file"

Features present in the BED data but missing from the color file render as gray. Common expected warnings include `novel` and `unknown`. If unexpected features appear, check that the `--colors` directory contains the correct color file for each featureset.

### Image is too wide

Reduce width by:
- Lowering `--n-per-cluster` (fewer reads per cluster)
- Increasing `--ratio` (fewer pixels per bp)
- Using `--target-width` to auto-calculate ratio
- Using `--max-qvalue 0.05` to show only significant clusters
- Switching to `--vertical` mode (trades width for height)

### Image is too tall

Reduce height by:
- Using fewer featuresets (e.g., `--featuresets region` instead of three)
- Reducing `--bar-width` from 8 to 4 or 6
- Using `--column-tracks` to arrange featuresets horizontally

### PNG export fails

`rsvg-convert` is required for PNG export. Install via:
```bash
# macOS
brew install librsvg

# Ubuntu/Debian
sudo apt-get install librsvg2-bin
```

### "Not enough reads for full dendrogram"

The `--full-dendrogram` option requires at least 2 reads with adjacency matrix data. Ensure the displayed reads exist in the `.feature_matrix.npz` file.

### Empty or missing reads

Check the log file for `Loaded data for N reads`. If N is lower than expected:
- Verify BED files exist at auto-discovered paths (check log for `{sample} -> {path}`)
- Ensure `--smoothness` matches the BED filename pattern
- Verify read names in `--reads-file` match those in the sequence assignments

---

## 13. Reference Tables

### 13.1 Featureset Display Names

| Internal Name | Display Name |
|--------------|-------------|
| `chromosome` | Chromosome |
| `subtelomeric` | Subtelomere |
| `region` | Satellite |
| `acrocentric` | Acrocentric |
| `repeat` | Interspersed repeat |
| `fiberseq` | Fiberseq |
| `fiberseq_FIRE` | FIRE |
| `fiberseq_m6A` | m6A |
| `fiberseq_5mC` | 5mC |
| `fiberseq_LINKER` | Linker |
| `fiberseq_FIRE_LINKER` | FIRE/Linker |
| `density_line` | m6A/5mC |
| `rect_plot` | FIRE/Linker |
| `telomere_region` | Telomere Region |

### 13.2 Complete Default Parameter Values

| Parameter | Default |
|-----------|---------|
| `--featuresets` | `chromosome,subtelomeric,region` |
| `--background` | `black` |
| `--bar-width` | `8` |
| `--bar-spacing` | `0` |
| `--read-spacing` | `12` |
| `--cluster-spacing` | `30` |
| `--ratio` | `0.00333` (1/300) |
| `--oversample` | `1` |
| `--smoothness` | `smoothed` |
| `--feature-mode` | `transition` |
| `--min-feature-width` | `0.5` |
| `--min-width-exclude` | `novel *arm* ct*` |
| `--density-bin-size` | `300` |
| `--dendro-cluster-gap` | `0` |
| `--structural-threshold` | `0.25` |
| `--label-column` | `curated_annotation` |
| `--font-family` | `sans-serif` |
| `--log-file` | `True` |
| `--png` | `False` |
| `--vertical` | `False` |
| `--column-tracks` | `False` |
| `--show-dendrogram` | `False` |
| `--hide-brackets` | `False` |
| `--hide-dendrogram` | `False` |
| `--no-reorder` | `False` |
| `--full-dendrogram` | `False` |
| `--show-matrix` | `False` |
| `--show-read-indices` | `False` |
| `--hide-read-labels` | `False` |
| `--show-cluster-numbers` | `False` |
| `--show-clade-id` | `False` |
| `--show-clade-count` | `False` |
| `--enrichment-grid` | `False` |
| `--orient-telomere-top` | `False` |
| `--show-threshold` | `False` |
