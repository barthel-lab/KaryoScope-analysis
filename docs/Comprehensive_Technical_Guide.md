# KaryoScope-Analyze: Comprehensive Technical Guide

## Overview

KaryoScope analysis is a graph-based method for analyzing sequencing data that identifies structural patterns based on genomic feature abundance and adjacencies as annotated by KaryoScope. This guide provides a detailed technical explanation of the clustering algorithm, from input data to final enrichment calls.

---

## Table of Contents

**1. [Input Data Structure](#1-input-data-structure)**

   - 1.1 [K-mer Coordinate System](#11-k-mer-classification-overview)
   - 1.2 [Feature Sets](#12-feature-sets)
   - 1.3 [Graph Representation](#14-graph-representation)
   - 1.4 [Sequence Filtering and Quality Control](#15-sequence-filtering-and-quality-control)

**2. [Feature Set Merging Strategies](#2-feature-set-merging-strategies)**

   - 2.1 [Overview and Rationale](#21-overview-and-rationale)
   - 2.2 [Concatenation Merging (Default Mode)](#22-concatenation-merging-default-mode)
   - 2.3 [Priority-Based Override Merging](#23-priority-based-override-merging)
   - 2.4 [Choosing a Merging Strategy](#24-choosing-a-merging-strategy)

**3. [Graph-to-Matrix Conversion](#3-graph-to-matrix-conversion)**

   - 3.1 [Feature Matrix Overview](#31-feature-matrix-overview)
     - 3.1.1 [Matrix Building Modes for Concatenation Merging](#311-matrix-building-modes-for-concatenation-merging)
   - 3.2 [Edge Variables](#32-edge-variables)
   - 3.3 [Abundance Variables](#33-abundance-variables)

**4. [Dimensionality Reduction (SVD)](#4-dimensionality-reduction-svd)**

   - 4.1 [When SVD is Applied](#41-when-svd-is-applied)
   - 4.2 [Truncated SVD Algorithm](#42-truncated-svd-algorithm)
   - 4.3 [Why Truncated SVD?](#43-why-truncated-svd)
   - 4.4 [Variance Explained Analysis](#44-variance-explained-analysis)

**5. [Hierarchical Clustering](#5-hierarchical-clustering)**

   - 5.1 [Distance Metric](#51-distance-metric)
   - 5.2 [Linkage Method](#52-linkage-method)
   - 5.3 [Dendrogram Structure](#53-dendrogram-structure)

**6. [Optimal k Selection](#6-optimal-k-selection)**

   - 6.1 [Grid Search Strategy](#61-grid-search-strategy)
   - 6.2 [Metrics Computed for Each k](#62-metrics-computed-for-each-k)
   - 6.3 [Composite Score](#63-composite-score)
   - 6.4 [K-Selection Methods](#64-k-selection-methods)
   - 6.5 [Example: Core-4 Dataset](#65-example-core-4-dataset)

**7. [Enrichment Testing](#7-enrichment-testing)**

   - 7.1 [Comparison Modes](#71-comparison-modes)
   - 7.2 [Enrichment Strength Categories](#72-enrichment-strength-categories)
   - 7.3 [Centroid Assignment](#73-centroid-assignment)

**8. [FDR Correction](#8-fdr-correction)**

   - 8.1 [Multiple Testing Problem](#81-multiple-testing-problem)
   - 8.2 [Benjamini-Hochberg Procedure](#82-benjamini-hochberg-procedure)
   - 8.3 [Impact Example (Core-4 Dataset)](#83-impact-example-core-4-dataset)
   - 8.4 [Alternative: Benjamini-Yekutieli](#84-alternative-benjamini-yekutieli)

**9. [Parameter Selection Guidelines](#9-parameter-selection-guidelines)**

   - 9.1 [Feature Matrix Parameters](#91-feature-matrix-parameters)
   - 9.2 [Dimensionality Reduction](#92-dimensionality-reduction)
   - 9.3 [Clustering Parameters](#93-clustering-parameters)
   - 9.4 [Enrichment Parameters](#94-enrichment-parameters)
   - 9.5 [Sequence Filtering](#95-sequence-filtering)

**Note:** For information on filtering and preprocessing the input sequences (e.g. for selecting telomere sequences) see `KaryoScope_Preprocessing.md`.

---

## 1. Input Data Structure {#1-input-data-structure}

### 1.1 K-mer Coordinate System {#11-k-mer-classification-overview}

A typical bed file is defined on the individual nucleotides of an input sequence. KaryoScope on the other hand is defined on the k-mers of an input sequence (by default, k=31). A position in this k-mer coordinate system corresponds to the k-mer which **starts** at that genomic position in the input. Note that a bed file in the k-mer coordinate system will include all but the final k-1 positions of every sequence in the input.

After k-mer annotation, consecutive k-mers with the same feature label are merged into intervals and represented in BED format:

```
sequence_id    start    end    feature
uuid           0        19     novel
uuid           19       182    q_arm_specific
uuid           182      183    ct_specific
uuid           183      314    q_arm_specific
```

In the above, the first 19 k-mers (starting at positions 1-19) are `novel` k-mers not found in the database, the next 163 k-mers (starting at positions 20-182) are `q_arm_specific`, the next k-mer (starting at position 183) is `ct_specific`, and the final 131 k-mers (starting at positions 184-314) are `q_arm_specific`. Note that in this representation, a feature "length" (i.e. `end` - `start`) actually corresponds to the number of k-mers.

### 1.2 Feature Sets {#12-feature-sets}

Features are organized into multiple independent feature sets, each providing a different annotation layer. KaryoScope is built on the 2,512,377,669 distinct k-mers of the CHM13 genome (excluding chrM). Each feature set is a database that assigns each of those 2,512,377,669 k-mers to a single feature. Any k-mer not in the database is assigned as the `novel` feature.

A feature set is built from a bed file with every position annotated as one of several initial features. However, k-mers may occur in multiple of these initial features (for example a `hsat1A` k-mer may also be a `hsat1B` k-mer). To ensure that k-mers are assigned to a single feature, KaryoScope constructs corresponding "feature_specific" and "feature_multigroup" sets of k-mers which are disjoint from each other.

For every initial feature, KaryoScope will construct a corresponding main "feature_specific" feature which consists of k-mers in that initial feature and no other initial feature. For k-mers belonging to multiple features, KaryoScope defines "feature_multigroup" features using a feature phylogeny. Each "feature_multigroup" feature consists of k-mers in a subset of initial features and not in any of the other initial features. For example, the `hsat1_multigroup1` feature consists of k-mers in `hsat1A` and `hsat1B` but no other initial features. 

KaryoScope consists of 6 feature sets:

**1. Chromosome feature set:** Chromosomal assignment

- Main Features: `chr1_specific`, `chr2_specific`, ..., `chr22_specific`, `chrX_specific`, `chrY_specific`
- Phylogeny-informed Features: `chromosome_multigroup1`, `autosome_multigroup1`, `acrocentric_multigroup1`, `metacentric_multigroup1`, `submetacentric_multigroup1`, `sex_multigroup1`
- File suffix: `.chromosome.bed`

**2. Region feature set:** Chromosomal regions and centromeric satellites

This feature set consist of two main groups: arm features and centromere features.

- Arm Features: `p_arm_specific`, `q_arm_specific`
- Centromere Features: `active_specific` (active aSat HOR), `inactive_specific` (inactive aSat HOR), `divergent_specific` (divergent aSat HOR), `monomeric_specific` (monomeric aSat), `hsat1A_specific`, `hsat1B_specific`, `hsat2_specific`, `hsat3_specific`, `bsat_specific` (beta satellite), `gsat_specific` (gamma satellite), `censat_specific` (other centromeric satellites), `ct_specific` (centromeric transition regions)
- rDNA Feature: `rDNA_specific`
- Phylogeny-informed Features: `arm_multigroup1`, `centromere_multigroup1`, `asat_multigroup1`, `hor_multigroup1`, `hsat_multigroup1`, `hsat1_multigroup1`
- File suffix: `.region.bed`

**3. Subtelomeric feature set:** Telomeric and subtelomeric sequences

- Main Features: `canonical_telomere_specific`, `noncanonical_telomere_specific`, `ITS_specific` (interstitial telomeric sequence), `TAR1_specific`
- Non-features: `nonsubtelomeric_specific`
- Phylogeny-informed Features: `telomere_like_multigroup1`
- File suffix: `.subtelomeric.bed`

**4. Repeat feature set:** Transposable elements and repetitive sequences

- Main Features: `rRNA_specific`, `scRNA_specific`, `snRNA_specific`, `tRNA_specific`, `LINE_specific`, `Retroposon_specific`, `SINE_specific`, `LTR_specific`, `DNA_specific`, `RC_specific`, `Unknown_specific`, `D20S16_specific`
- Non-features: `nonrepeat_specific`
- Phylogeny-informed Features: `LINE-dependent_Retroposon_multigroup1`, `Class_I_Retrotransposition_multigroup1`, `Class_II_DNA_Transposition_multigroup1`, `Transposable_Element_multigroup1`, `Interspersed_Repeat_multigroup1`, `repeat_multigroup1`
- File suffix: `.repeat.bed`

**5. Gene feature set:** Coding region annotations

- Main Features: `exon_specific`, `intron_specific`
- Non-features: `intergenic_specific`
- File suffix: `.gene.bed`

**6. Acrocentric feature set:** Acrocentric sequences

- Main Features: `DJ_specific` (distal junction of rDNA array), `PJ_specific` (proximal junction of rDNA array), `rDNA_specific`, `SST1_specific`, `PHR_specific` (pseudo-homologous regions)
- Non-features: `nonacrocentric_specific`
- Phylogeny-informed Features: `array_multigroup1`
- File suffix: `.acrocentric.bed`

**Note:** KaryoScope annotations are output to a separate BED file for each feature set. One may wish to merge multiple feature sets such that every genomic interval is annotated with the corresponding features from all the feature sets being merged (see [Feature Set Merging Strategies](#2-feature-set-merging-strategies)).

### 1.3 Graph Representation {#14-graph-representation}

We construct a **path graph** for each annotated sequence:

- **Nodes** = features
- **Edges** = transitions between adjacent features along the sequence
- **Edge weights** = average length (i.e. number of k-mers) of the two adjacent features forming the transition: (length_feature_i + length_feature_j) / 2

**Example path graph:**
```
Sequence: uuid
Path: novel → q_arm_specific → ct_specific → q_arm_specific

Nodes: {novel, q_arm_specific, ct_specific}
Edges: {novel→q_arm_specific, q_arm_specific→ct_specific, ct_specific→q_arm_specific}
Edge weights: {19, 82.5, 0.5}
```

This graph representation enables clustering sequences by their structural patterns (which features appear and how they connect).

### 1.4 Sequence Filtering and Quality Control {#15-sequence-filtering-and-quality-control}

Sequences undergo filtering to ensure optimized input for clustering (default parameters shown below, modifiable via script arguments):

**Length filtering (default):**

- Minimum: 10,000 bp (10 kb) - removes short, potentially fragmented sequences
- Maximum: 100,000 bp (100 kb) - removes rare ultra-long outliers
- Rationale: Focuses on high-quality long sequences with sufficient feature content

**Feature exclusion:**
Sequences containing only the following features are excluded:

- `novel` - sequences not present in the KaryoScope database
- `canonical_telomere*` - canonical telomeric repeats (using wildcard matching)
  - Rationale: Canonical telomere abundance is sample-dependent and can bias clustering results unnecessarily.

**Example filtering impact (Core-4 telomeric sequences dataset):**
```
Total records: 6,850,068
After feature exclusion: 6,476,740 (94.5%)
Sequences before length filter: 12,357
Sequences after length filter: 8,422 (68.1%)
Final dataset: 8,422 sequences for clustering
```

This dual filtering ensures:

1. Structurally informative sequences (not just canonical telomeres)
2. High-quality sequences suitable for statistical analysis
3. Computationally tractable datasets

---

## 2. Feature Set Merging Strategies {#2-feature-set-merging-strategies}

### 2.1 Overview and Rationale {#21-overview-and-rationale}

When combining multiple KaryoScope feature sets (e.g., chromosome + region + repeat), overlapping intervals from different feature sets must be merged to create unified feature annotations for each genomic interval.

**Why merge feature sets?**

Merging enables the measurement of associations between independent genomic annotation layers, revealing biologically meaningful structural patterns:

- **Juxtapositions:** Detecting adjacency between satellite DNA (region feature set) and subtelomeric features (subtelomeric feature set)
- **Context-specific repeats:** Identifying repeat elements (repeat feature set) enriched in specific chromosomal locations (chromosome feature set)
- **Multi-dimensional patterns:** Capturing complex genomic architectures that span multiple annotation types

**Example biological question:**
"Are certain satellite arrays (e.g. hsat, bsat) preferentially found adjacent to specific subtelomeric elements (e.g. noncanonical_telomere, TAR1)?"

This requires merging the region feature set (containing satellite annotations) with the subtelomeric feature set (containing non-canonical telomere/TAR1 annotations) to identify juxtapositions.

### 2.2 Concatenation Merging (Default Mode) {#22-concatenation-merging-default-mode}

The default strategy preserves information from all feature sets by concatenating feature labels.

**Algorithm:**

1. Identify overlapping intervals between feature sets
2. For each overlapping region, concatenate feature labels using a separator (default: `:`)
3. Create a new interval with the combined label

**Comprehensive example showing increased granularity:**

Since each feature set annotates every k-mer (using `novel` as necessary), the merged result will annotate every k-mer with a feature from each input feature set.

```
Position:       0         100       200       300       400       500       600
                |         |         |         |         |         |         |
Chromosome:     [chr7───────────────────────────────────────────────────────]
Region:         [p_arm──────────────][bsat──────────────][q_arm─────────────]
Repeat:         [nonrepeat][LINE────][SINE────][nonrepeat───────────────────]

Merged result (fragmented into non-overlapping intervals):
Interval 1:  0-100    → chr7:p_arm:nonrepeat
Interval 2:  100-200  → chr7:p_arm:LINE
Interval 3:  200-300  → chr7:bsat:SINE
Interval 4:  300-400  → chr7:bsat:nonrepeat
Interval 5:  400-600  → chr7:q_arm:nonrepeat
```

**Key observations:**

- **Every interval contains a feature from each input feature set**
- The merged BED file is **more granular than** any single input BED file
- Intervals are split wherever feature boundaries from different feature sets differ
- Number of merged intervals is determined by union of all boundary positions across all feature sets

**Output characteristics:**

- Preserves complete information from all input feature sets
- Typical format: `chromosome:region:repeat` (e.g. `chr7:p_arm:LINE`)
- Increases feature dimensionality (more unique feature combinations)
- Higher fragmentation = more granular annotation

**Implementation:** `KaryoScope_merge_beds.py` (default mode)

### 2.3 Priority-Based Override Merging {#23-priority-based-override-merging}

An alternative strategy is to use conditional logic to select features from one feature set over another based on biological priority rules. This reduces dimensionality at the cost of information loss.

#### 2.3.1 Telomere-Satellite Priority Merge (`--telomere-satellite-merge`)

This mode prioritizes subtelomeric features over region/centromeric features.

**Priority features from the subtelomeric feature set:**

- `canonical_telomere_specific`
- `noncanonical_telomere_specific`
- `TAR1_specific`
- `ITS_specific`

**Algorithm:**

1. Identify priority features from subtelomeric feature set (BED1)
2. Identify all features from satellite/region feature set (BED2)
3. Keep priority features in their entirety
4. Fill remaining (non-priority) positions with satellite/region features
5. Discard any satellite/region features that overlap with priority features

**Example showing merging decisions at every position:**

Since both feature sets annotate every position, a merging decision must be made at every interval.

```
Position:       0          100        200        300        400        500        600        700
                |          |          |          |          |          |          |          |
Subtelomeric:   [canonical_telomere───][nonsubtelomeric─────][ITS──────][nonsubtelomeric─────]
Region:         [hsat1A─────────────────────────────────────][bsat───────────────────────────]

Decision:       priority               non-priority          priority   non-priority
                feature                (use region)          feature    (use region)

Merged result:  [canonical_telomere───][hsat1A──────────────][ITS─────][bsat─────────────────]
                0-200                  200-400               400-500   500-700
```

Merged result:
```
sequence_id    start    end    feature
uuid           0        200    canonical_telomere_specific
uuid           200      400    hsat1A_specific
uuid           400      500    ITS_specific
uuid           500      700    bsat_specific
```

#### 2.3.2 3-Way Priority Merge (`--priority-merge`)

This mode prioritizes subtelomeric features over region features and region features over repeat features with a few additional conditional rules.

**Conditional rules:**

- `nonsubtelomeric` + `ct` + repeat (not `nonrepeat`) → use repeat feature
- `nonsubtelomeric` + `rDNA` + `rRNA` → `rRNA`
- `nonsubtelomeric` + `p_arm` + repeat → use repeat feature
- `nonsubtelomeric` + `q_arm` + repeat → use repeat feature

**Algorithm:**

1. Extract priority subtelomeric features (`canonical_telomere`, `noncanonical_telomere`, `ITS`, `TAR1`, `telomere_like_multigroup1`)
2. Merge region + repeat using conditional rules
3. Subtract priority subtelomeric regions from region+repeat merged intervals
4. Combine priority subtelomeric + remaining region/repeat intervals

**Output characteristics:**

- Reduces feature complexity by selecting single labels
- Loses information from lower-priority feature sets
- Biologically-informed feature selection
- Lower feature dimensionality
- Simpler interpretation

**Implementation:** `KaryoScope_merge_beds.py` (with `--priority-merge` or `--telomere-satellite-merge` flags)

### 2.4 Choosing a Merging Strategy {#24-choosing-a-merging-strategy}

**Use concatenation merging (default) when:**

- You want to preserve all information
- Computational resources allow higher dimensionality
- Discovering novel feature associations is the goal
- You need to distinguish `chr7:hsat:LINE` from `chr7:hsat:SINE`

**Use priority-based merging when:**

- You want simpler, lower-dimensional features
- Specific features (e.g. subtelomeric) are biologically more important
- Computational efficiency is critical
- You're willing to sacrifice information for interpretability

---

## 3. Graph-to-Matrix Conversion {#3-graph-to-matrix-conversion}

Each sequence is converted into a numerical feature matrix that captures **composition** (what features are present) alongside **structure** (how they connect). This dual representation enables clustering sequences by their telomeric architecture patterns.

### 3.1 Feature Matrix Overview {#31-feature-matrix-overview}

The feature matrix is constructed by concatenating two types of variables:

```
Feature Matrix = [Edge Variables | Abundance Variables]
Shape: (n_sequences, n_edge_variables + n_abundance_variables)
```

**Edge variables** represent transitions between consecutive genomic features along the sequence, capturing structural patterns.

**Abundance variables** represent the proportion of each genomic feature within the sequence, capturing compositional patterns.

**Example dimensions:**

```
- 8,422 sequences
- 276 edge variables (unique transitions)
- 24 abundance variables (unique genomic features)
- Final shape: (8422, 300)
```

**Sparsity:** Typically 90-98% sparse (most edges don't occur in most sequences).

### 3.1.1 Matrix Building Modes for Concatenation Merging {#311-matrix-building-modes-for-concatenation-merging}

When using **concatenation merging** (which produces colon-separated multi-layer features), KaryoScope supports two matrix building strategies:

**Note:** These modes are **not applicable** to priority-based merging, which produces single-layer features without colons (e.g., `canonical_telomere`, `hsat`, `ITS`).

**Combined Mode (default):**

- Treats merged features as atomic units
- Example: `chr7:p_arm:LINE` is a single feature
- Creates edges like: `chr7:p_arm:LINE` → `chr7:q_arm:SINE`
- **Advantage:** Preserves exact feature combinations from merging
- **Trade-off:** Higher dimensionality (more unique feature combinations)

**Layered Mode:**

- Splits colon-separated features into independent layers
- Example: `chr7:p_arm:LINE` splits into 3 layers: chromosome=`chr7`, region=`p_arm`, repeat=`LINE`
- Builds separate edge + abundance matrices per layer
- Concatenates layer matrices: `[chr_matrix | region_matrix | repeat_matrix]`
- **Advantage:** Lower dimensionality (features split across layers)
- **Trade-off:** Loses information about specific feature co-occurrences

**Example comparison:**

```
Input: Two reads with merged features
Read 1: chr7:p_arm:LINE → chr7:q_arm:SINE
Read 2: chr7:p_arm:SINE → chr7:q_arm:LINE

Combined mode creates edges:
  - chr7:p_arm:LINE → chr7:q_arm:SINE
  - chr7:p_arm:SINE → chr7:q_arm:LINE
  (2 distinct edges, captures full feature combinations)

Layered mode creates edges (3 layers):
  Layer 1 (chromosome): chr7 → chr7, chr7 → chr7
  Layer 2 (region): p_arm → q_arm, p_arm → q_arm
  Layer 3 (repeat): LINE → SINE, SINE → LINE
  (Same edges within each layer, loses cross-layer associations)
```

**When to use each:**

- **Combined (default)**: When feature co-occurrences are important (e.g., distinguishing `chr7:hsat:LINE` from `chr7:hsat:SINE`)
- **Layered**: When dimensionality is too high or you want layer-independent analysis

**Parameter:** `--matrix-mode` (default: `combined`)

### 3.2 Edge Variables {#32-edge-variables}

Edge variables capture the structural organization of features along each sequence by representing transitions between consecutive feature pairs.

#### Edge Extraction

For each sequence's path graph, consecutive feature pairs are extracted to create **transition edges**.

#### Edge Modes

**Directional** (preserves transition order):
```
Sequence path: q_arm → ct → q_arm

Edges counted:
  q_arm → ct (1st transition)
  ct → q_arm (2nd transition)

Each transition creates a distinct, ordered edge.
```

**Bidirectional** (counts both forward and reverse):
```
Sequence path: q_arm → ct

Edges counted:
  q_arm → ct (forward edge)
  ct → q_arm (reverse edge, also added)

Each physical transition generates two directional edges in the matrix.
```

**Symmetric** (default, order-independent):
```
Sequence path: q_arm → ct → q_arm

Edges counted:
  ct->q_arm (alphabetically sorted: ct < q_arm)
  ct->q_arm (same edge, merged)

Result: Single edge with weight = 2
```

Symmetric mode reduces dimensionality by treating A→B and B→A as equivalent, while preserving biological transition patterns.

#### Edge Weighting

Three matrix types are available:

**Binary Matrix:**
```
edge_value = 1 if transition occurs at least once in the sequence, else 0
```

Indicates presence/absence of the transition regardless of frequency. Useful for comparing sequence structural patterns without length or repeat count bias.

**Count Matrix:**
```
edge_value = number of times transition occurs in the sequence
```

Counts the total number of times each transition is observed within the sequence. Captures transition frequency information but may be biased by sequence length or repetitive structures.

**Length-Weighted Matrix (default):**
```python
for each transition (feature_i → feature_j):
    edge_weight = length_i / total_sequence_length
```

For **bidirectional mode**, edges use asymmetric weights:
```python
forward_edge (feature_i → feature_j):
    edge_weight = length_i / total_sequence_length

reverse_edge (feature_j → feature_i):
    edge_weight = (length_i + length_j) / (2 × total_sequence_length)
```

**Normalization:** All weights are divided by sequence length to make sequences of different lengths comparable.

**Edge accumulation:** When the same transition occurs multiple times in a sequence:

- **Binary**: Set to 1 (no accumulation, presence/absence only)
- **Count**: Increment by 1 each occurrence (accumulation: `edge_value += 1`)
- **Length-weighted**: Add normalized weight each occurrence (accumulation: `edge_value += weight`)

### 3.3 Abundance Variables {#33-abundance-variables}

Abundance variables capture the compositional content of each sequence by measuring the proportion of each genomic feature:

```python
for each feature in sequence:
    abundance[feature] = total_feature_length / sequence_length
```

This provides complementary information to edge variables: while edges capture feature connectivity patterns, abundances capture overall feature content.

---

## 4. Dimensionality Reduction (SVD) {#4-dimensionality-reduction-svd}

### 4.1 When SVD is Applied {#41-when-svd-is-applied}

SVD is applied **only when** the feature dimensionality exceeds a threshold (default: 500 dimensions).

```python
# Determine effective number of components (constrained by data)
n_components = min(reduce_dims, n_sequences - 1, n_features)

# Apply SVD only if reduction would actually occur
if n_features > n_components:
    apply_truncated_svd(n_components=n_components)
```

**Constraints:**

- Cannot exceed `n_sequences - 1` (mathematical limit of SVD)
- Cannot exceed `n_features` (no point in "reducing" to more dimensions)
- Target dimensions specified by `--reduce-dims` (default: 500)

**Examples:**

- 8,000 sequences, 1,200 features, target=500 → reduce to 500
- 8,000 sequences, 300 features, target=500 → no reduction (300 < 500)
- 400 sequences, 1,200 features, target=500 → reduce to 399 (n_sequences - 1)

### 4.2 Truncated SVD Algorithm {#42-truncated-svd-algorithm}

**Input:** Feature matrix X[n × d], target dimensions k=500

**Decomposition:**
```
X ≈ U × Σ × V^T
```

Where:

- **U[n × k]**: Sequence coordinates in SVD space (used for clustering)
- **Σ[k × k]**: Singular values (diagonal, ordered by importance)
- **V^T[k × d]**: Feature loadings (how original features combine)

**Output:** Transformed matrix U × Σ[n × k]

### 4.3 Why Truncated SVD? {#43-why-truncated-svd}

**Computational advantages:**

1. **Efficiency:** Computes only top k singular vectors (not full decomposition)
2. **No centering:** Works directly on sparse matrices (unlike PCA)
3. **Dimensionality reduction:** 80%+ reduction (e.g., 2415 → 500 dims)
4. **Noise filtering:** Minor components discarded

**Benefits for clustering:**

1. **Better distance metrics:** Euclidean distance more meaningful in lower dimensions
2. **Computational efficiency:** Faster distance matrix calculation
3. **Curse of dimensionality mitigation:** Reduces distance concentration
4. **Redundancy removal:** Correlated features combined

### 4.4 Variance Explained Analysis {#44-variance-explained-analysis}

Example from Core-6 dataset:
```
Original dimensions: 2415
Target dimensions: 500
Total explained variance: 100.0%
Components for 50% variance: 2
Components for 90% variance: 8
Components for 95% variance: 15
Top 5 singular values: 65.7, 53.9, 26.1, 20.3, 14.2
Top 5 variance %: 23.0%, 40.2%, 11.4%, 7.0%, 3.4%
```

**Interpretation:** Only 2 components capture 50% of variance, indicating strong low-dimensional structure in telomere sequence patterns.

---

## 5. Hierarchical Clustering {#5-hierarchical-clustering}

### 5.1 Distance Metric {#51-distance-metric}

**Euclidean distance** in the (optionally SVD-reduced) feature space:

```python
dist_matrix = pdist(feature_matrix, metric='euclidean')
```

For two sequences i and j:

```
d(i,j) = sqrt(Σ(x_i,k - x_j,k)^2) for k=1 to n_features
```

### 5.2 Linkage Method {#52-linkage-method}

**Ward's method** (default) minimizes within-cluster variance:

```python
linkage_matrix = linkage(dist_matrix, method='ward')
```

**Ward's criterion:** At each merge, choose the pair of clusters that results in minimum increase in total within-cluster sum of squares.

**Distance measure:** Ward uses **squared Euclidean distance** internally.

### 5.3 Dendrogram Structure {#53-dendrogram-structure}

The linkage matrix encodes the hierarchical tree:

```
[cluster_1, cluster_2, distance, n_members]
```

This can be cut at different heights to produce k clusters.

---

## 6. Optimal k Selection {#6-optimal-k-selection}

### 6.1 Grid Search Strategy {#61-grid-search-strategy}

**Range:** min_k to max_k (default: 40 to 300)

**Constraint:** max_k ≤ n_sequences / 10 (ensures ≥10 sequences per cluster on average)

**Early stopping:** Stops if no improvement in composite score for N iterations (default: 150)

Example:
```
Testing cluster counts from 40 to 300 (261 values)
Early stopping at k=200 (no improvement for 150 iterations)
```

### 6.2 Metrics Computed for Each k {#62-metrics-computed-for-each-k}

For each candidate k, the dendrogram is cut and the following metrics are calculated:

#### A. Clustering Quality Metrics

**Silhouette Score** (range: -1 to 1, higher is better):
```python
silhouette = silhouette_score(feature_matrix, cluster_labels)
```
Measures cluster cohesion and separation.

**Calinski-Harabasz Index** (higher is better):
```python
CH = calinski_harabasz_score(feature_matrix, cluster_labels)
```
Ratio of between-cluster to within-cluster variance.

**Davies-Bouldin Index** (optional, lower is better):
```python
DB = davies_bouldin_score(feature_matrix, cluster_labels)
```
Average similarity between clusters and their most similar neighbor.

#### B. Biological Enrichment Metrics

For each cluster, perform Fisher's exact test (or Chi-square for multi-group):

```python
# Fast enrichment check during k-optimization
for cluster in clusters:
    odds_ratio = calculate_odds_ratio(cluster_samples, all_samples)
    if odds_ratio > 1.5:  # Meaningful enrichment threshold
        enrichment_type = "sample-enriched"
    else:
        enrichment_type = "mixed"

    # Classify by purity
    purity = max_sample_percentage(cluster)
    if purity >= 0.95:
        category = "perfect"
    elif purity >= 0.80:
        category = "strong"
    else:
        category = "any"
```

**Counts:**

- `any_enriched`: Clusters with p < 0.05 (any enrichment)
- `strong_enriched`: Clusters with ≥80% from one sample
- `perfect_enriched`: Clusters with ≥95% from one sample

**Ratios:**
```python
enriched_ratio = any_enriched / valid_clusters
strong_ratio = strong_enriched / valid_clusters
perfect_ratio = perfect_enriched / valid_clusters
```

### 6.3 Composite Score {#63-composite-score}

Weighted combination balancing statistical quality and biological interpretability:

```python
silhouette_norm = (silhouette + 1) / 2  # Normalize from [-1,1] to [0,1]

composite_score = (0.5 × silhouette_norm +
                   0.1 × enriched_ratio +
                   0.4 × perfect_ratio)
```

**Weights:**

- 50%: Cluster quality (silhouette)
- 10%: Any enrichment ratio
- 40%: Perfect enrichment ratio (emphasizes high-purity clusters)

### 6.4 K-Selection Methods {#64-k-selection-methods}

Four methods available for selecting optimal k:

#### Method 1: Maximum Silhouette
```python
k_optimal = argmax(silhouette_scores)
```
Favors fewer, tighter clusters. Often selects lower k values.

#### Method 2: Maximum Calinski-Harabasz
```python
k_optimal = argmax(CH_scores)
```
Favors compact, well-separated clusters.

#### Method 3: Maximum Composite
```python
k_optimal = argmax(composite_scores)
```
Balances quality and enrichment. May favor higher k if it improves enrichment.

#### Method 4: Composite-Knee (default)

**Kneedle Algorithm** finds the "elbow" point where increasing k yields diminishing returns:

```python
# 1. Normalize both axes to [0,1]
k_norm = (k_values - k_min) / (k_max - k_min)
score_norm = (scores - score_min) / (score_max - score_min)

# 2. Calculate distance from diagonal
knee_distance = score_norm - k_norm

# 3. Apply smoothing (20% rolling window)
window_size = max(3, int((k_max - k_min) * 0.2))
knee_smooth = rolling_mean(knee_distance, window=window_size)

# 4. Find maximum distance
k_optimal = k_values[argmax(knee_smooth)]
```

**Geometric interpretation:** The knee point is where the curve is farthest from a straight diagonal line.

**Why this works:**

- Before the knee: rapid improvement in quality per additional cluster
- At the knee: optimal trade-off point

**Important note on normalization:** The knee point detection normalizes based on the **observed range** of k values and scores during the search. This means the detected knee position may vary slightly if `max_k` is changed, because normalization bounds depend on the actual scores computed. This is expected behavior - the knee represents the diminishing returns point within the explored range.

- After the knee: diminishing returns (overfitting)

### 6.5 Example: Core-4 Dataset {#65-example-core-4-dataset}

```
Testing k from 40 to 300, stopped at k=200

Optimal k by metric:
  Silhouette:        k=42 (score=0.1947)
  Calinski-Harabasz: k=40 (score=2263.6)
  Composite:         k=50 (score=0.4426)
  Composite-knee:    k=51 (smoothed, window=32)

At k=51:
  Valid clusters:    51
  Any enriched:      43 (84.3%)
  Strong enriched:   9 (17.6%)
  Perfect enriched:  8 (15.7%)

Selected: k=51 (composite-knee method)
```

**Rationale:** k=51 achieves 84.3% enriched clusters while balancing cluster quality. Beyond this point, quality improvements are marginal.

---

## 7. Enrichment Testing {#7-enrichment-testing}

### 7.1 Comparison Modes {#71-comparison-modes}

Three modes available depending on experimental design:

#### Mode 1: Two-Group Comparison

**Use case:** Control vs. Treatment (e.g., normal vs. tumor)

**Statistical test:** Fisher's exact test

**Contingency table:**
```
                In Cluster    Out of Cluster
Control         a             b
Treatment       c             d
```

**Test:**
```python
odds_ratio, p_value = fisher_exact([[a, c], [b, d]])
```

**Enrichment assignment:**
```python
if p_value < 0.05:
    if odds_ratio > 1:
        enrichment = "control-enriched"
    elif odds_ratio < 1:
        enrichment = "treatment-enriched"
else:
    enrichment = "mixed"
```

**Key advantage:** Odds ratio correctly handles unbalanced group sizes (e.g., 90% tumor, 10% normal).

#### Mode 2: Multi-Group Comparison

**Use case:** 3+ groups without a designated control

**Statistical test:** Chi-square test of independence

**Contingency table:**
```
               Group1  Group2  Group3  ...
In cluster     n1      n2      n3    ...
Out cluster    m1      m2      m3    ...
```

**Test:**
```python
chi2, p_value, dof, expected = chi2_contingency(contingency)
```

**Enrichment assignment:**
```python
if p_value < 0.05:
    dominant_group = max(groups, key=lambda g: percentage_in_cluster[g])
    enrichment = f"{dominant_group}-enriched"
else:
    enrichment = "mixed"
```

#### Mode 3: Per-Sample Comparison

**Use case:** Each sample tested individually vs. all others (most granular)

**Statistical test:** Multiple Fisher's exact tests (one per sample)

**For each sample:**
```python
# 2×2 contingency table
in_cluster = sample_count_in_cluster
out_cluster = sample_total - in_cluster
other_in = cluster_size - in_cluster
other_out = total_sequences - sample_total - other_in

odds, p_val = fisher_exact([[in_cluster, other_in],
                             [out_cluster, other_out]])
```

**Enrichment assignment:**
```python
# Find most significant ENRICHED sample (p < 0.05 AND odds > 1)
enriched_samples = {s: p for s, p in p_values.items()
                    if p < 0.05 and odds_ratios[s] > 1}

if enriched_samples:
    best_sample = min(enriched_samples, key=lambda s: p_values[s])
    enrichment = f"{best_sample}-enriched"
else:
    enrichment = "mixed"
```

**Important:** Only considers **over-representation** (odds > 1), not depletion, for enrichment calls.

### 7.2 Enrichment Strength Categories {#72-enrichment-strength-categories}

Clusters are classified by purity:

```python
purity = max_sample_percentage / 100

if purity >= 0.95:
    category = "perfect"    # ≥95% from one sample
elif purity >= 0.80:
    category = "strong"     # ≥80% from one sample
elif p_value < 0.05:
    category = "any"        # Significant but mixed
else:
    category = "mixed"      # Not significant
```

**Example clusters:**

**Perfect enrichment (100% U2OS):**
```
Cluster 26: 134 sequences
  U2OS: 134 (100.0%)
  p-value: 4.91e-68
  Category: perfect
```

**Strong enrichment (84% BJ):**
```
Cluster 4: 306 sequences
  BJ: 113 (36.9%)
  IMR90: 46 (15.0%)
  HeLa: 78 (25.5%)
  U2OS: 69 (22.5%)
  p-value: 4.92e-09
  Category: strong
```

**Mixed (no enrichment):**
```
Cluster 21: 507 sequences
  BJ: 128 (25.2%)
  IMR90: 94 (18.5%)
  HeLa: 164 (32.3%)
  U2OS: 121 (23.9%)
  p-value: 4.38e-05 (but no dominant sample)
  Category: mixed
```

### 7.3 Centroid Assignment {#73-centroid-assignment}

For each cluster, the **centroid sequence** is identified:

```python
# Calculate cluster centroid (mean feature vector)
cluster_matrix = feature_matrix[cluster_members]
centroid = cluster_matrix.mean(axis=0)

# Find sequence closest to centroid
distances = euclidean_distance(cluster_matrix, centroid)
centroid_sequence = cluster_members[argmin(distances)]
```

**Purpose:** Representative sequence for visualization and validation.

**Warning flag:** If centroid sequence's sample differs from enrichment label, flagged with ⚠️ (indicates heterogeneous cluster).

---

## 8. FDR Correction {#8-fdr-correction}

### 8.1 Multiple Testing Problem {#81-multiple-testing-problem}

With 51 clusters, we perform 51 statistical tests. At α=0.05, we expect ~2.5 false positives by chance.

**Solution:** False Discovery Rate (FDR) correction using Benjamini-Hochberg method.

### 8.2 Benjamini-Hochberg Procedure {#82-benjamini-hochberg-procedure}

```python
# 1. Collect all raw p-values
raw_pvals = [cluster_1_pval, cluster_2_pval, ..., cluster_51_pval]

# 2. Apply FDR correction
q_values = false_discovery_control(raw_pvals, method='bh')

# 3. Reclassify enrichment based on q-values
for cluster in clusters:
    if q_value < 0.05:
        enrichment = original_enrichment  # Keep direction
    else:
        enrichment = "mixed"  # Lost significance
```

**Algorithm:**
1. Sort p-values: p₁ ≤ p₂ ≤ ... ≤ pₘ
2. Find largest i such that: pᵢ ≤ (i/m) × α
3. Reject hypotheses 1 through i
4. Calculate q-values: qᵢ = pᵢ × (m/i)

### 8.3 Impact Example (Core-4 Dataset) {#83-impact-example-core-4-dataset}

```
Raw p < 0.05: 43 enriched clusters
FDR q < 0.05: 41 enriched clusters
2 clusters lost significance after FDR correction
```

**Interpretation:** 95% confidence that ≤5% of the 41 enriched clusters are false discoveries.

### 8.4 Alternative: Benjamini-Yekutieli {#84-alternative-benjamini-yekutieli}

For dependent tests (more conservative):

```python
q_values = false_discovery_control(raw_pvals, method='by')
```

**When to use:** If clusters are not independent (e.g., nested hierarchical structure).

---

## 9. Parameter Selection Guidelines {#9-parameter-selection-guidelines}

### 9.1 Feature Matrix Parameters {#91-feature-matrix-parameters}

**Matrix Type:**

- `length_weighted` (default): Edge weight = source feature length normalized by sequence length. Best for biological interpretation as it reflects the relative genomic span of transitions.
- `count`: Edge weight = number of times transition occurs in the sequence. Use if feature lengths are unreliable or when transition frequency is more informative than length.
- `binary`: Edge weight = 1 if transition occurs at least once in the sequence, else 0. Use for presence/absence patterns only, ignoring frequency and length.

**Matrix Mode (for concatenation-merged features only):**

- `combined` (default): Treats colon-separated features as atomic (e.g., `chr7:p_arm:LINE` is one feature)
- `layered`: Splits into independent layers and concatenates layer matrices
- **Note:** Only applies to concatenation merging; not relevant for priority-merged features

**Edge Mode:**

- `symmetric` (default): Reduces dimensions, order-independent
- `directional`: Use if transition direction matters (e.g., strand-specific)
- `bidirectional`: Use to emphasize transitions (creates asymmetric weights)

**Abundance:**

- `True` (default): Recommended for comprehensive feature representation
- `False`: Use only if interested purely in transitions

### 9.2 Dimensionality Reduction {#92-dimensionality-reduction}

**SVD threshold:**

- `500` (default): Good balance
- `300`: Use for smaller datasets (faster, less smoothing)
- `1000`: Use for very large feature spaces (more robust)
- `0`: Disable (only if d < 300)

### 9.3 Clustering Parameters {#93-clustering-parameters}

**Linkage method:**

- `ward` (default): Minimizes variance, produces balanced clusters
- `average`: More robust to outliers
- `complete`: Creates compact clusters (high within-cluster similarity)
- `single`: Creates elongated clusters (use with caution)

**K-selection:**

- `composite-knee` (default): Best for discovery (finds diminishing returns)
- `composite`: Use if you want maximum enrichment
- `silhouette`: Use if quality more important than enrichment
- `calinski`: Similar to silhouette but emphasizes separation

**K-range:**

- `min_k=40` (default): Ensures granular clustering
- `max_k=300` (default): Broad search range
- Adjust based on dataset size: aim for 10-100 sequences per cluster on average

**Early stopping:**

- `150` (default): Good for most cases
- `50`: Faster, less thorough
- `0`: Disable (test full range)

### 9.4 Enrichment Parameters {#94-enrichment-parameters}

**Comparison mode:**

- `per-sample`: Most granular, recommended for initial exploration
- `two-group`: Use if clear control/treatment design
- `multi-group`: Use for 3+ groups without designated control

**FDR threshold:**

- `0.05` (default): Standard stringency
- `0.10`: More permissive (discovery mode)
- `0.01`: More stringent (validation mode)

**Enrichment thresholds:**

- `perfect_threshold=0.95`: Very stringent
- `strong_threshold=0.80`: Moderate stringency

### 9.5 Sequence Filtering {#95-sequence-filtering}

**Length filters:**

- `min_sequence_length=10000`: Removes noisy short sequences
- `max_sequence_length=100000`: Removes rare long outliers
- Adjust based on expected biological sequence length distribution

**Feature exclusion:**

- Exclude low-quality annotations: `novel,unknown`
- Exclude canonical telomeres: `canonical_telomere*`
- Use wildcards for pattern matching

---

## Summary Pipeline

```
1. Input: BED files (path graphs)
   ↓
2. Extract edge and abundance variables
   → Edge matrix (transitions, length-weighted)
   → Abundance matrix (feature proportions)
   ↓
3. Combine: Feature matrix [n_sequences × n_variables]
   → Typical: 8000 sequences × 300 variables (7% dense)
   ↓
4. [Optional] SVD dimensionality reduction
   → If n_features > 500: reduce to 500 dims
   → Retain 100% variance (noise filtered)
   ↓
5. Hierarchical clustering (Ward + Euclidean)
   → Build linkage matrix (dendrogram)
   ↓
6. K-selection (test k=40 to 300)
   → For each k: compute quality + enrichment metrics
   → Composite score: 50% silhouette + 10% enrichment + 40% purity
   → Find knee point (diminishing returns)
   ↓
7. Cut dendrogram at optimal k
   → Typical: k=40-60 clusters
   ↓
8. Enrichment testing
   → Fisher's exact or Chi-square per cluster
   → Classify: perfect/strong/any/mixed
   ↓
9. FDR correction (Benjamini-Hochberg)
   → q-value threshold: 0.05
   → Typical: 80-90% clusters remain significant
   ↓
10. Output: Cluster assignments + enrichment calls
```

---

## References

- Ward, J. H. (1963). Hierarchical grouping to optimize an objective function. *Journal of the American Statistical Association*, 58(301), 236-244.
- Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics*, 20, 53-65.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society*, Series B, 57(1), 289-300.
- Satopaa, V., et al. (2011). Finding a "Kneedle" in a Haystack: Detecting Knee Points in System Behavior. *IEEE ICDCSW*, 166-171.
