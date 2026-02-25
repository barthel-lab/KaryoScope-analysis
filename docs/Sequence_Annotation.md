# Sequence Annotation

## Overview

`KaryoScope_sequence_annotate.py` computes **per-read feature annotations** from BED files and outputs a wide-format TSV with one row per read. This is the first step in the two-step annotation pipeline:

1. **`KaryoScope_sequence_annotate.py`** (this script): per-read BED → sequence annotations TSV
2. **`KaryoScope_cluster_annotate.py`**: sequence annotations → cluster annotations TSV

By separating the resource-intensive per-read computations (window densities, interspersion, feature fractions) from the lightweight cluster-level aggregations, the per-read step only needs to run once and its output can be reused across multiple cluster annotation runs.

### Inputs

| Input | Source | Description |
|-------|--------|-------------|
| `--prefix` | cluster_analysis output | Path prefix for `.read_assignments.tsv` |
| `--bed-dir` | KaryoScope BED output | Comma-separated directories containing per-sample BED files |
| `--output` | User-specified | Output TSV path |

### Outputs

| Output | Description |
|--------|-------------|
| `{output}` | Main per-read annotations TSV (one row per read) |
| `{output_prefix}.adaptive_thresholds.tsv` | Per-feature adaptive thresholds |
| `{output_prefix}.log` | Log file (when `--log-file` enabled) |

---

## 1. Per-Read Metrics

### 1.1 Feature Fractions and BP Counts

For each featureset and each feature within it:

| Column | Type | Description |
|--------|------|-------------|
| `{fs}_frac__{feat}` | float (0–1) | Feature fraction: feature bp / total annotated bp |
| `{fs}_bp__{feat}` | int | Feature base pairs |
| `{fs}_total_bp` | int | Total annotated bp for this read in this featureset |

### 1.2 Window Density Metrics

Computed by `compute_per_read_window_densities_bulk()`. For each read, a binary coverage array is constructed per feature, and a sliding window of 1 kb is applied. Values are stored as raw fractions (0–1).

| Column | Type | Description |
|--------|------|-------------|
| `{fs}_dmax__{feat}` | float (0–1) | Max window density |
| `{fs}_dmin__{feat}` | float (0–1) | Min window density |
| `{fs}_dmedian__{feat}` | float (0–1) | Median window density |
| `{fs}_dfirst__{feat}` | float (0–1) | First 1 kb density |
| `{fs}_dlast__{feat}` | float (0–1) | Last 1 kb density |
| `{fs}_dterminal__{feat}` | float (0–1) | max(first, last) |
| `{fs}_dterminal_min__{feat}` | float (0–1) | min(first, last) |
| `{fs}_max_block_bp__{feat}` | int (bp) | Longest contiguous feature block (gaps ≤ 100 bp merged) |

**Short reads (< 1 kb):** All window stats default to the overall coverage fraction of the read.

### 1.3 Interspersion Metrics

Computed from the `telomere_region` featureset (fallback: `region`). Features are classified into categories (canonical, noncanonical, satellite, arm, ITS_TAR1, ct, other) and transitions between adjacent categories are counted per kilobase.

| Column | Type | Description |
|--------|------|-------------|
| `interspersion_total` | float | All category transitions per kb |
| `interspersion_can_ncan` | float | canonical ↔ noncanonical transitions per kb |
| `interspersion_tel_sat` | float | telomere ↔ satellite transitions per kb |
| `interspersion_arm_tel` | float | arm ↔ telomere transitions per kb |

### 1.4 Optional Alignment Statistics

When `--readnames-dir` is provided, alignment statistics from readnames.txt and stats.tsv files are joined onto each read:

| Column | Type | Description |
|--------|------|-------------|
| `sequencing_approach` | str | From readnames.txt (e.g., "hifi", "ont") |
| `n_alignments` | int | Total number of alignments |
| `n_secondary` | int | Number of secondary alignments |
| `n_supplementary` | int | Number of supplementary alignments |
| `primary_mapq` | int | Mapping quality of primary alignment |
| `primary_de` | float | Divergence of primary alignment |
| `primary_align_len` | int | Aligned bases in primary alignment |
| `primary_align_fraction` | float | Fraction of read aligned in primary alignment |
| `total_align_len` | int | Sum of aligned bases across all alignments |
| `total_align_fraction` | float | Total aligned / read length |

---

## 2. Adaptive Thresholds

`compute_adaptive_thresholds()` sets a per-feature significance threshold:

```
threshold = clamp(median_nonzero / 3, min=0.001, max=0.05)
```

The thresholds TSV is saved alongside the main output and used by `cluster_annotate` to compute `readpct` scores.

| Column | Description |
|--------|-------------|
| `featureset` | Featureset name |
| `feature` | Feature name |
| `threshold` | Computed threshold value |
| `median_nonzero` | Median of nonzero fractions |
| `n_nonzero` | Number of reads with nonzero fraction |
| `n_total` | Total number of reads |

---

## 3. Relationship to Cluster Annotation

The sequence annotations TSV is consumed by `KaryoScope_cluster_annotate.py` via the `--sequence-annotations` argument. The cluster annotator:

1. Reads the pre-computed per-read columns
2. Aggregates them to cluster-level metrics (medians, sums, percentages)
3. Scales density values from 0–1 to 0–100 for the cluster output

This separation means you can re-run cluster annotation (e.g., with different `--auto-label` settings or `--min-size` filters) without recomputing the expensive per-read metrics.

---

## 4. CLI Reference

```
python KaryoScope_sequence_annotate.py \
    --prefix ANALYSIS_PREFIX \
    --bed-dir BED_DIRS \
    --output OUTPUT_TSV \
    [--featuresets FEATURESETS] \
    [--database DATABASE] \
    [--smoothness SMOOTHNESS] \
    [--window-size WINDOW_SIZE] \
    [--readnames-dir READNAMES_DIR] \
    [--reference REFERENCE] \
    [--log-file LOG_FILE]
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--prefix` | Yes | — | Analysis prefix (finds `{prefix}.read_assignments.tsv`) |
| `--bed-dir` | Yes | — | Comma-separated BED directories |
| `--output` / `-o` | Yes | — | Output TSV path |
| `--featuresets` | No | `region,subtelomeric,chromosome,acrocentric,repeat,gene` | Comma-separated featuresets |
| `--database` | No | `KS_human_CHM13` | Database name |
| `--smoothness` | No | `smoothed` | BED smoothness |
| `--window-size` | No | `1000` | Window size in bp |
| `--readnames-dir` | No | `None` | Enables alignment stats from readnames.txt + stats.tsv |
| `--reference` | No | `CHM13` | Reference genome name for stats.tsv files |
| `--log-file` | No | `True` | Save .log file |

---

## 5. Example Commands

### Basic usage

```bash
python KaryoScope_sequence_annotate.py \
    --prefix results/analysis \
    --bed-dir results \
    --output results/analysis.sequence_annotations.tsv
```

### With alignment statistics

```bash
python KaryoScope_sequence_annotate.py \
    --prefix results/analysis \
    --bed-dir results \
    --readnames-dir /data/samples \
    --output results/analysis.sequence_annotations.tsv
```

### Specific featuresets

```bash
python KaryoScope_sequence_annotate.py \
    --prefix results/analysis \
    --bed-dir results \
    --featuresets region,subtelomeric,chromosome \
    --output results/analysis.sequence_annotations.tsv
```

### Full two-step pipeline

```bash
# Step 1: Per-read annotations (run once)
python KaryoScope_sequence_annotate.py \
    --prefix results/analysis \
    --bed-dir results \
    --output results/analysis.sequence_annotations.tsv

# Step 2: Cluster annotations (can re-run with different settings)
python KaryoScope_cluster_annotate.py \
    --prefix results/analysis \
    --sequence-annotations results/analysis.sequence_annotations.tsv \
    --auto-label \
    --output results/analysis.cluster_annotations.tsv
```
