# KaryoScope Preprocessing: Telomeric Read Filtering

## Overview

This document describes the preprocessing steps that prepare long-read telomeric sequencing reads for KaryoScope clustering analysis. The primary preprocessing step is **telomeric read filtering** using Telogator2, which enriches for reads containing telomeric and subtelomeric sequences. While not strictly required, this filtering dramatically improves computational efficiency and signal-to-noise ratio for telomeric structural variant detection.

---

## Table of Contents

**1. [Telomeric Read Filtering (Telogator2)](#1-telomeric-read-filtering-telogator2)**
   - 1.1 [Purpose](#11-purpose)
   - 1.2 [Algorithm](#12-algorithm)
   - 1.3 [Output Files](#13-output-files)
   - 1.4 [Read Type Annotations](#14-read-type-annotations)
   - 1.5 [Integration with KaryoScope Workflow](#15-integration-with-karyoscope-workflow)
   - 1.6 [Filtering Impact on Dataset Size](#16-filtering-impact-on-dataset-size)

**2. [Alternative: Direct Input Without Telogator](#2-alternative-direct-input-without-telogator)**

**3. [Quality Control Recommendations](#3-quality-control-recommendations)**

**4. [Citation](#4-citation)**

---

## 1. Telomeric Read Filtering (Telogator2)

### 1.1 Purpose

Telogator2 identifies and extracts reads containing telomeric sequences from whole-genome long-read sequencing data. This enrichment step:

- Reduces dataset size by 10-100x depending on sequencing depth
- Focuses downstream analysis on structurally informative telomeric reads
- Removes reads lacking telomeric content
- Enables efficient telomere-focused structural variant detection

### 1.2 Algorithm

Telogator2 uses k-mer-based classification to identify telomeric reads:

**Primary filtering criterion:**

Reads are classified as telomeric if they contain ≥ `MINIMUM_CANON_HITS` occurrences of tandem canonical telomeric 12-mers.

**Default parameters:**
- `MINIMUM_CANON_HITS`: 8 (minimum 12-mer hits)
- `MINIMUM_READ_LEN`: 4000 bp (minimum read length)
- Base k-mer size: 6 bp (doubled to 12 bp for tandem search)

**Canonical telomeric 12-mers (tandem 6-mers):**

- **Primary (all platforms):** `CCCTAACCCTAA` (forward), `TTAGGGTTAGGG` (reverse complement)
  - Constructed by repeating the canonical 6-mer `CCCTAA` / `TTAGGG` twice
- **Additional (Nanopore only):** Tandem repeats of `CCTGG` and `TTTTAA` (improve sensitivity for Nanopore error profiles)

**Search algorithm:**

The algorithm searches for exact matches to these 12-mer patterns (6-mer repeated twice) using simple string counting. Reads with ≥8 occurrences of the canonical 12-mer pattern are retained.

**Example:**
```
Read sequence:  ...TTAGGGTTAGGGTTAGGGTTAGGG...CCCTAACCCTAA...
12-mer matches:   [TTAGGGTTAGGG]
                        [TTAGGGTTAGGG]
                              [TTAGGGTTAGGG]
                                    [TTAGGGTTAGGG] [CCCTAACCCTAA]

Total 12-mer hits: 5 (below threshold of 8, read excluded)

Read sequence:  ...TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG...
12-mer matches:   [TTAGGGTTAGGG]
                        [TTAGGGTTAGGG]
                              [TTAGGGTTAGGG]
                                    [TTAGGGTTAGGG]
                                          [TTAGGGTTAGGG]
                                                [TTAGGGTTAGGG]
                                                      [TTAGGGTTAGGG]
                                                            [TTAGGGTTAGGG]

Total 12-mer hits: 8 (meets threshold, read retained)
```

**Note:** The algorithm uses Python's `.count()` method to search for the 12-mer string `CCCTAACCCTAA` or its reverse complement `TTAGGGTTAGGG`, counting overlapping occurrences.

### 1.3 Output Files

Telogator2 generates multiple output files. The key file for KaryoScope input is:

**`tel_reads.fa.gz`** - Filtered telomeric reads in FASTA format
- Contains all reads meeting the tandem canonical k-mer threshold
- Ready for input to KaryoScope feature annotation pipeline
- Typical size reduction: 10-100x smaller than input

**Additional outputs (optional, used for telomere length analysis):**
- `tlens_by_allele.tsv` - Allele-specific telomere lengths
- `all_final_alleles.png` - Plots of all telomeric alleles
- `violin_atl.png` - Violin plots of allelic telomere lengths

### 1.4 Read Type Annotations

KaryoScope workflows support organizing reads by sequencing technology and library preparation using user-defined type annotations:

**Common PacBio type annotations:**
- `hifi_fiber`: High-fidelity reads from Fiber-seq libraries (chromatin accessibility profiling)
- `hifi_notfiber`: Standard PacBio HiFi reads

**Common Oxford Nanopore type annotations:**
- `ONT_UL`: Ultra-long reads (typically >50 kb)
- `ONT_nonfrag`: Non-fragmented standard-length reads
- `ONT_frag`: Fragmented or sheared reads

**Implementation:**

Type annotations are specified in two locations:

1. **Workflow configuration (`config.yaml`):**
   - Users manually assign type labels when defining input files
   - The `{type}` wildcard in workflow files (e.g., `{sample}.{type}.{i}`) refers to these user-defined labels
   - Example:
   ```yaml
   data:
     sample1:
       ONT_UL:
         1: /path/to/ultralong_reads.fastq.gz
       hifi_fiber:
         1: /path/to/fiberseq_reads.bam
   ```

2. **Sequencing approach annotation (`readnames.txt` files):**
   - For post-clustering annotation, users create `{sample}.readnames.txt` files
   - Tab-separated format: `read_name\tsequencing_approach`
   - Values: `hifi`, `ont`, or other user-defined approaches
   - Location: `{readnames_dir}/{sample}/telogator/{sample}.readnames.txt`
   - Used by `KaryoScope_annotate_sequences.py` to add sequencing approach metadata to clustering results

**Note:** Type classifications are independent of Telogator2 filtering and based entirely on user knowledge of library preparation and sequencing strategy.

### 1.5 Integration with KaryoScope Workflow

The telogator-filtered `tel_reads.fa.gz` file serves as input to the KaryoScope feature annotation pipeline:

**Workflow structure:**
```
tel_reads.fa.gz (from Telogator2)
    ↓
KaryoScope feature annotation (get_features_reads script)
    ↓
results/{sample}/telogator/{i}/features/{database}/
    {sample}.telogator.{i}.{database}.combined.featureIDs.clustered.bed.gz
    ↓
Feature smoothing (smooth_features.py)
    ↓
results/{sample}/telogator/{i}/KaryoScope/{database}/
    {sample}.telogator.{i}.{database}.{featureset}.smoothed.features.bed
    ↓
Feature merging (KaryoScope_merge_beds.py)
    ↓
Clustering analysis (KaryoScope_cluster_analysis.py)
```

**Directory structure:**
```
results/{sample}/telogator/{i}/
├── features/{database}/           # Raw feature annotations
│   └── {sample}.telogator.{i}.{database}.combined.featureIDs.clustered.bed.gz
└── KaryoScope/{database}/         # Processed BED files for clustering
    ├── {sample}.telogator.{i}.{database}.chromosome.smoothed.features.bed
    ├── {sample}.telogator.{i}.{database}.region.smoothed.features.bed
    ├── {sample}.telogator.{i}.{database}.repeat.smoothed.features.bed
    ├── {sample}.telogator.{i}.{database}.subtelomeric.smoothed.features.bed
    └── {sample}.telogator.{i}.{database}.{merged_name}.smoothed.merged.bed
```

**Iteration number `{i}`:**

The `{i}` parameter is an index for multiple sequencing runs on the same sample. For example:
- `telogator/1/`: First telogator run or first input file
- `telogator/2/`: Second telogator run or second input file (e.g., different sequencing batch)

This allows processing multiple datasets or parameter sets in parallel within the same workflow.

### 1.6 Filtering Impact on Dataset Size

**Example: Core4 Dataset**

Post-telogator filtering results:
- **BJ:** 2,681 telomeric reads
- **HeLa:** 3,474 telomeric reads
- **IMR90:** 2,197 telomeric reads
- **U2OS:** 4,005 telomeric reads
- **Total:** 12,357 telomeric reads

**Technology distribution:**
- ~50% PacBio HiFi, ~50% ONT ultra-long
- Classification breakdown: hifi_fiber (25%), hifi_notfiber (28%), ONT_UL (47%)

Read counts vary by sample based on:
- Biological telomere content (e.g., ALT-positive cells may have more telomeric reads)
- Sequencing depth
- Telomere length distribution

---

## 2. Alternative: Direct Input Without Telogator

KaryoScope clustering can be performed on any set of annotated long reads without telogator preprocessing. However, this approach requires significantly more computational resources (10-100x more reads to process).

To use non-filtered reads, provide the reads directly to the KaryoScope feature annotation pipeline (`get_features_reads` script) as described in the KaryoScope Comprehensive Guide.

---

## 3. Quality Control Recommendations

**Before telogator filtering:**
- Ensure reads are basecalled with recent algorithms (Guppy 5+, Dorado SUP, or SMRTLink 13+)
- Older Nanopore basecallers (Guppy <5) have high error rates in telomeric regions and may fail

**After telogator filtering:**
- Verify expected telomeric read counts match sample type and sequencing depth
- Check for platform biases in read count distributions
- Inspect `tlens_by_allele.tsv` for reasonable telomere length distributions (if generated)

**Warning signs:**
- Very low telomeric read counts (<100 reads) may indicate:
  - Insufficient sequencing depth
  - Sample quality issues
  - Incorrect telogator parameters
- Very high telomeric read counts (>50,000 reads) may indicate:
  - Telomere-enriched sequencing
  - ALT-positive sample

---

## 4. Citation

If using Telogator2 for preprocessing, please cite:

Stephens, Z., & Kocher, J. P. (2024). Characterization of telomere variant repeats using long reads enables allele-specific telomere length estimation. *BMC Bioinformatics*, 25(1), 194.
