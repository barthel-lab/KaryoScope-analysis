# Cluster Annotation

## Overview

`KaryoScope_cluster_annotate.py` takes the output of `KaryoScope_cluster_analysis.py` and annotates each cluster with structural metrics derived from the raw BED feature data. It computes density profiles, coverage statistics, interspersion rates, and optionally assigns automatic cluster names via a rule-based decision tree.

### Inputs

| Input | Source | Description |
|-------|--------|-------------|
| `--prefix` | cluster_analysis output | Path prefix for `.sequence_assignments.tsv` and `.cluster_analysis.tsv` |
| `--bed-dir` | KaryoScope BED output | Directories containing per-sample BED files |
| `--output` | User-specified | Output TSV path |

### Output

A TSV file with one row per cluster containing: basic cluster info, enrichment statistics, density metrics, coverage metrics, interspersion metrics, and optionally auto-assigned cluster names and representative reads.

---

## 1. Density Metrics

!!! note "Smoothed BED files"
    All density, coverage, and interspersion metrics are computed from the **smoothed** BED files (controlled by `--smoothness`, default `"smoothed"`). The same smoothed BED files feed all window-based, coverage-based, and interspersion computations.

### 1.1 Window-Based Density Metrics

Computed by `compute_cluster_window_densities()`. For each feature, a **binary coverage array** is constructed per read (1 where the feature is present, 0 elsewhere), then a sliding window of 1 kb is applied. Per-read statistics are collected, and the **cluster-level value is the median across all reads** (zero-padded for reads missing the feature). All values are scaled 0--100 (percent).

| Metric | Column Pattern | Formula | Meaning |
|--------|---------------|---------|---------|
| **dmax** | `{fs}_dmax__{feat}` | max(window_sums) / window_size x 100 | Peak density in any 1 kb window |
| **dmin** | `{fs}_dmin__{feat}` | min(window_sums) / window_size x 100 | Lowest density in any 1 kb window |
| **dmedian** | `{fs}_dmedian__{feat}` | median(window_sums) / window_size x 100 | Central tendency across all windows |
| **dfirst** | `{fs}_dfirst__{feat}` | coverage[:1000].sum() / 1000 x 100 | Density in first 1 kb of the read |
| **dlast** | `{fs}_dlast__{feat}` | coverage[-1000:].sum() / 1000 x 100 | Density in last 1 kb of the read |
| **dterminal** | `{fs}_dterminal__{feat}` | max(dfirst, dlast) x 100 | Stronger terminal end |
| **dterminal_min** | `{fs}_dterminal_min__{feat}` | min(dfirst, dlast) x 100 | Weaker terminal end |
| **max_block_bp** | `{fs}_max_block_bp__{feat}` | Longest contiguous covered stretch (gaps <= 100 bp merged) | Longest unbroken feature block (in bp, not %) |

**Short reads (< 1 kb):** All window stats default to the overall coverage fraction of the read.

**Aggregation:** Per-read values are collected and zero-padded for reads missing a feature. The **median** is then taken across all reads in the cluster.

#### Worked Example: Window-Based Metrics

Consider a single 5 kb read with a feature covering positions 500--2500 (2 kb block) and 4000--4800 (0.8 kb block):

```
0         500                          2500       4000         4800  5000
|---------|============FEATURE==========|---gap----|===FEAT=====|----|
0000000000111111111111111111111111111111110000000000011111111111110000
```

**Step 1: Binary coverage array.** 5000 positions, with 1s at positions 500--2499 and 4000--4799.

**Step 2: Sliding 1 kb window sums.** There are 4001 windows (positions 0--999, 1--1000, ..., 4000--4999). Each window sum is the count of covered positions in that 1 kb span.

Example windows:

| Window start | Positions covered | Sum |
|-------------|-------------------|-----|
| 0 | 500--999 | 500 |
| 500 | 500--1499 | 1000 |
| 1500 | 1500--2499 | 1000 |
| 2000 | 2000--2499 | 500 |
| 2500 | (none) | 0 |
| 3500 | 4000--4499 | 500 |
| 4000 | 4000--4799 | 800 |

**Step 3: Compute metrics.**

| Metric | Calculation | Value |
|--------|-------------|-------|
| **dmax** | max window sum / 1000 x 100 = 1000/1000 x 100 | **100.0%** |
| **dmin** | min window sum / 1000 x 100 = 0/1000 x 100 | **0.0%** |
| **dmedian** | median of all 4001 window sums / 1000 x 100 | **~50%** |
| **dfirst** | coverage[0:1000].sum() / 1000 x 100 = 500/1000 x 100 | **50.0%** |
| **dlast** | coverage[4000:5000].sum() / 1000 x 100 = 800/1000 x 100 | **80.0%** |
| **dterminal** | max(50.0, 80.0) | **80.0%** |
| **dterminal_min** | min(50.0, 80.0) | **50.0%** |
| **max_block_bp** | Longest contiguous block (gap of 1500 bp > 100 bp breaks it) | **2000 bp** |

!!! note "Cluster aggregation"
    For a cluster of N reads, each read produces these per-read values. The cluster metric is the **median** across reads, with zeros padded for reads that lack the feature entirely.

---

### 1.2 Coverage-Based Metrics

#### readpct

`score_cluster_features()` computes the percentage of reads in a cluster where a feature exceeds an adaptive significance threshold.

**Column:** `{fs}_readpct__{feat}`

**Formula:** `100 x (reads with feature fraction > adaptive_threshold) / cluster_size`

#### bppct

`compute_cluster_bp_scores()` computes the base-pair-weighted percentage of a feature across all reads in the cluster.

**Column:** `{fs}_bppct__{feat}`

**Formula:** `100 x feature_bp / total_bp`

#### Adaptive Thresholds

`compute_adaptive_thresholds()` sets a per-feature significance threshold to distinguish real signal from noise:

```
threshold = clamp(median_nonzero / 3, min=0.001, max=0.05)
```

Where `median_nonzero` is the median of all nonzero feature fractions across reads.

#### Worked Example: Coverage-Based Metrics

Consider a 3-read cluster where Feature X covers varying fractions:

```
Read 1 (4 kb): |==X (2kb)==|----------|                    fraction = 0.50
Read 2 (6 kb): |X(1k)|-------------------------|          fraction = 0.167
Read 3 (5 kb): |-------------------------|                 fraction = 0.00
```

Each character between `|` delimiters represents 200 bp (Read 1 = 20 chars, Read 3 = 25 chars, Read 2 = 30 chars). `=` = feature X, `-` = non-feature.

**Adaptive threshold:**

1. Nonzero fractions: 0.50, 0.167
2. Median of nonzero = median(0.50, 0.167) = 0.333
3. Threshold = 0.333 / 3 = 0.111
4. Clamped to [0.001, 0.05] &rarr; **0.05**

**readpct:**

- Read 1: 0.50 > 0.05 &check;
- Read 2: 0.167 > 0.05 &check;
- Read 3: 0.00 > 0.05 &cross;
- readpct = 2/3 x 100 = **66.7%**

**bppct:**

- Total feature bp: 2000 + 1000 + 0 = 3000
- Total bp: 4000 + 6000 + 5000 = 15000
- bppct = 3000 / 15000 x 100 = **20.0%**

---

### 1.3 Interspersion Metrics

Computed by `compute_cluster_interspersion()`. Counts typed transitions between adjacent BED features per kilobase. Reports the cluster median across reads.

Features are first classified into categories by `classify_bed_feature()`:

| Category | Examples |
|----------|---------|
| canonical | canonical_telomere |
| noncanonical | noncanonical_telomere |
| satellite | active, monomeric, bsat, hsat1A, hsat2, hsat3, gsat, censat |
| arm | p_arm, q_arm, arm_multigroup |
| ITS_TAR1 | ITS, TAR1 |
| ct | ct (centric transition) |
| other | Everything else (excluded from transition counts) |

| Metric | Column | Transition types counted |
|--------|--------|------------------------|
| **interspersion_total** | `interspersion_total` | Any category change (excluding 'other') |
| **interspersion_can_ncan** | `interspersion_can_ncan` | canonical &harr; noncanonical |
| **interspersion_tel_sat** | `interspersion_tel_sat` | (canonical or noncanonical) &harr; satellite |
| **interspersion_arm_tel** | `interspersion_arm_tel` | arm &harr; (canonical or noncanonical or ITS_TAR1) |

#### Worked Example: Interspersion Metrics

Consider a single 5 kb read with this feature layout:

```
|--CAN--|--ncan--|--SAT--|--CAN--|--ARM--|
0      1000    2000    3000    4000    5000
```

Categories: canonical, noncanonical, satellite, canonical, arm

Transitions (excluding 'other'):

| # | From | To | total | can_ncan | tel_sat | arm_tel |
|---|------|----|-------|----------|---------|---------|
| 1 | canonical | noncanonical | x | x | | |
| 2 | noncanonical | satellite | x | | x | |
| 3 | satellite | canonical | x | | x | |
| 4 | canonical | arm | x | | | x |

span_kb = 5.0

| Metric | Calculation | Value |
|--------|-------------|-------|
| **interspersion_total** | 4 / 5.0 | **0.80 transitions/kb** |
| **interspersion_can_ncan** | 1 / 5.0 | **0.20 transitions/kb** |
| **interspersion_tel_sat** | 2 / 5.0 | **0.40 transitions/kb** |
| **interspersion_arm_tel** | 1 / 5.0 | **0.20 transitions/kb** |

!!! note "Cluster aggregation"
    For a cluster, each read produces these per-read rates. The cluster metric is the **median** across reads.

---

## 2. Output Columns

### Fixed columns

| Column | Source | Type |
|--------|--------|------|
| `cluster_id` | assignments file | int |
| `size` | Number of reads in cluster | int |
| `cluster_name` | auto_label or empty for user curation | str |
| `curated_rep_i` | User curation or 1 if reps selected | int or '' |

### Enrichment statistics (from cluster_analysis.tsv)

| Column | Type | Notes |
|--------|------|-------|
| `enrichment` | str | e.g. 'U2OS-enriched', 'mixed' |
| `q_value` | str (scientific notation) | Multiple-testing corrected |
| `log2_fc` | float | log2(odds_ratio), only if odds_ratio column exists |
| `{sample}_count` | int | Per-sample read count |
| `{sample}_pct` | float | Per-sample percentage (1 decimal) |
| `{sample}_pval` | str (scientific notation) | Per-sample p-value |
| `{sample}_odds` | float | Per-sample odds ratio (2 decimals) |

### Interspersion (telomere_region featureset only)

`interspersion_total`, `interspersion_can_ncan`, `interspersion_tel_sat`, `interspersion_arm_tel`

### Per-featureset summary

| Column | Type | Example |
|--------|------|---------|
| `{fs}_top` | str | `"canonical_telomere(45.2%); active(23.1%); ..."` |

### Per-feature metrics (for every feature in every featureset)

| Column | Description |
|--------|-------------|
| `{fs}_readpct__{feat}` | % of reads with feature above threshold |
| `{fs}_bppct__{feat}` | % of total bp covered by feature |
| `{fs}_dmax__{feat}` | Peak 1 kb window density |
| `{fs}_dmin__{feat}` | Minimum 1 kb window density |
| `{fs}_dmedian__{feat}` | Median 1 kb window density |
| `{fs}_dfirst__{feat}` | First 1 kb density |
| `{fs}_dlast__{feat}` | Last 1 kb density |
| `{fs}_dterminal__{feat}` | max(dfirst, dlast) |
| `{fs}_dterminal_min__{feat}` | min(dfirst, dlast) |
| `{fs}_max_block_bp__{feat}` | Longest contiguous block in bp |

### Representative reads (optional, with `--select-representatives N`)

`representative_read_1` through `representative_read_N`

---

## 3. Cluster Name Assignment

When `--auto-label` is enabled, `auto_label_cluster()` assigns each cluster a structural label using a rule-based decision tree.

### 3.1 Thresholds

| Constant | Value | Purpose |
|----------|-------|---------|
| CAN_ECTR | 70 | Canonical terminal density for ECTR (each end) |
| NCAN_ECTR | 10 | Noncanonical terminal density for ECTR (each end) |
| CAN_SUB | 15 | Canonical terminal density for Subtelomere (one end) |
| NCAN_SUB | 5 | Noncanonical terminal density for Subtelomere (one end) |
| DMAX_HIGH | 35 | dmax to count feature as "present" |
| ENRICH_DMAX | 35 | dmax threshold for enrichment qualifiers |
| NCAN_VARIANT | 25 | Noncanonical dmax for variant-enriched qualifier |
| SAT_DOMINANT | 80 | readpct for satellite-dominant (rule 5) |
| CT_ENRICH | 20 | {pfx}_bppct__ct for SegDup enrichment |
| ALT_BLOCK_BP | 6000 | max_block_bp for Type II ALT subtelomere |
| ARM_PRESENT | 30 | Arm dmax to count as having arm sequence |

### 3.2 Decision Tree

Rules are evaluated in priority order. The first matching rule determines the label.

!!! info "Chromosomal anchor"
    **Chromosomal anchor** means the read extends beyond the subtelomeric region into a chromosome arm (p_arm or q_arm). Detected when `max(p_arm_dmax, q_arm_dmax) >= 30` or `{pfx}_bppct__ct >= 20`. Reads with a chromosomal anchor have subtelomeric structure on one end and chromosomal sequence on the other, distinguishing Subtelomeres from ECTRs.

An end is considered **telomeric** if its canonical density meets the canonical threshold OR its noncanonical density meets the noncanonical threshold. ECTR uses stricter thresholds (CAN_ECTR=70, NCAN_ECTR=10) checked independently at each end; Subtelomere uses looser thresholds (CAN_SUB=15, NCAN_SUB=5) on the stronger terminal end.

#### Rule 1: ECTR

**Rule 1 — ECTR:** Telomere detected at **both** ends of the read (canonical >= 70% or noncanonical >= 10% at each end). Also includes reads with telomere at one end but **no chromosomal anchor** — these are free-floating telomeric DNA rather than chromosome-attached subtelomeres.

#### Rule 2: Subtelomere

**Rule 2 — Subtelomere:** At least one end telomeric (canonical >= 15% or noncanonical >= 5%), **with** a chromosomal anchor. If the longest contiguous canonical telomere block is >= 6 kb, the cluster is labeled "Type II ALT subtelomere" instead of plain "Subtelomere". Enrichment qualifiers are appended.

#### Rule 3: Interstitial telomere

**Rule 3 — Interstitial telomere:** High telomere dmax (canonical or noncanonical >= 35) but no telomere at either terminal end. These reads contain telomeric sequence buried in the interior. Enrichment qualifiers are appended.

#### Rule 4: Interstitial ITS/TAR1

**Rule 4 — Interstitial ITS/TAR1:** No telomere structure, but ITS and/or TAR1 enriched (dmax >= 35). The label is determined by which features are present:

- **Both** ITS and TAR1 dmax >= 35 → "Interstitial ITS/TAR1"
- **Only** TAR1 dmax >= 35 → "Interstitial TAR1"
- **Only** ITS dmax >= 35 → "Interstitial ITS"

Satellite, rDNA, and SegDup qualifiers are appended. The variant qualifier is **not** applied for interstitial ITS/TAR1 clusters.

#### Rule 5: Satellite dominant

**Condition:** A single satellite feature is present in >= 80% of reads.

```
IF max_satellite_readpct >= 80:
    label = satellite_name    (e.g. "HSat1A", "active aSat")
```

Note: No enrichment qualifiers are appended for satellite-dominant clusters.

#### Rule 6: Unlabeled

If no rule matches, the cluster receives an empty label (`""`).

### 3.3 Enrichment Qualifiers

`_enrichment_qualifiers()` builds a list of enrichment tags appended to the base label. Checked in order:

1. **variant** -- noncanonical dmax > 25
2. **{satellite name}** -- max satellite dmax >= 35 (picks highest by dmax)
3. **TAR1/ITS** -- both >= 35 together, or individually
4. **rDNA** -- acrocentric_dmax__rDNA >= 35
5. **SegDup** -- {pfx}_bppct__ct >= 20

**Format:** `" (variant-enriched)"` for one qualifier, or `" (variant-, HSat3-enriched)"` for multiple.

### 3.4 Post-Processing: Type I ALT Relabeling

When `--alt-samples` is provided, clusters where the ALT sample percentage exceeds the threshold (default 80%) have "Type I ALT " prepended to their label. This applies only to clusters that are **not** already labeled ECTR or Type II ALT.

Example: "Subtelomere (variant-enriched)" &rarr; "Type I ALT subtelomere (variant-enriched)"

---

## 4. Naming Examples

### Example A: ECTR

A cluster of reads with canonical telomere at both ends:

```
|====CAN====|---ncan---|====CAN====|   (5 kb read)
0           1500       3000        5000
```

| Metric | Value |
|--------|-------|
| can_dfirst | 100.0 |
| can_dlast | 100.0 |
| ncan_dfirst | 0.0 |
| ncan_dlast | 0.0 |
| can_dterminal | 100.0 |
| has_arm | False |

**Decision trace:**

1. end_is_tel(first): can_dfirst=100 >= CAN_ECTR=70 &rarr; True
2. end_is_tel(last): can_dlast=100 >= CAN_ECTR=70 &rarr; True
3. ectr = True AND True = **True**
4. **Label: "ECTR"**

### Example B: Subtelomere

A cluster with telomere at one end and arm sequence at the other:

```
|===ARM===|--ct--|=======CAN=======|   (5 kb read)
0         1500  2500               5000
```

| Metric | Value |
|--------|-------|
| can_dfirst | 0.0 |
| can_dlast | 100.0 |
| can_dterminal | 100.0 |
| p_arm_dmax | 100.0 |
| {pfx}_bppct__ct | 20.0 |
| has_arm | True (p_arm_dmax=100 >= ARM_PRESENT=30) |

**Decision trace:**

1. ectr: end_is_tel(first) = False (can_dfirst=0 < 70) &rarr; ectr = False
2. sub: end_is_tel(terminal) = True (can_dterminal=100 >= 15)
3. sub AND NOT has_arm = True AND False = False &rarr; Rule 1 fails
4. Rule 2: sub=True, has_arm=True &rarr; check max_block_bp
5. can_max_block < 6000 &rarr; **Label: "Subtelomere (SegDup-enriched)"**

### Example C: Interstitial Telomere

A cluster with telomere in the middle, flanked by chromosome arm sequence:

```
|---arm---|=====CAN=====|---arm---|   (5 kb read)
0         1500          3500      5000
```

| Metric | Value |
|--------|-------|
| can_dfirst | 0.0 |
| can_dlast | 0.0 |
| can_dmax | 100.0 |
| p_arm_dmax | 100.0 |

**Decision trace:**

1. ectr: both ends False &rarr; ectr = False
2. sub: dterminal = max(0, 0) = 0 &rarr; sub = False
3. Rule 1 and 2 fail (no terminal telomere)
4. Rule 3: tel_dmax = max(can_dmax=100, ncan_dmax=0) = 100 >= 35 &rarr; True
5. **Label: "Interstitial telomere"**

### Example D: Satellite Dominant

A cluster of reads that are almost entirely HSat3:

```
|===================HSat3===================|   (5 kb read)
0                                           5000
```

| Metric | Value |
|--------|-------|
| can_dfirst | 0.0 |
| can_dlast | 0.0 |
| can_dmax | 0.0 |
| hsat3_readpct | 100.0 |

**Decision trace:**

1. Rules 1--4 all fail (no telomere, no ITS/TAR1)
2. Rule 5: max_sat_readpct = 100 >= SAT_DOMINANT=80 &rarr; True
3. **Label: "HSat3"**

### Example E: Type II ALT Subtelomere

A cluster with arm sequence and a long contiguous canonical telomere block:

```
|==ARM==|=============CAN (7 kb block)==============|   (9 kb read)
0       2000                                        9000
```

| Metric | Value |
|--------|-------|
| can_dlast | 100.0 |
| can_dterminal | 100.0 |
| can_max_block_bp | 7000 |
| has_arm | True |

**Decision trace:**

1. ectr: end_is_tel(first) = False &rarr; ectr = False
2. sub = True (can_dterminal=100 >= 15)
3. Rule 2: sub=True, has_arm=True, can_max_block=7000 >= ALT_BLOCK_BP=6000
4. **Label: "Type II ALT subtelomere"**
