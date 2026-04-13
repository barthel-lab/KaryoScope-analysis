# KaryoScope — HPRC centromere clustering analysis

## Project description
K-mer based annotation of centromere repeat features across the HPRC dataset.
KaryoScope (KS) annotates genomic sequences by pattern-matching against a
k-mer database. Input: FASTQ/FASTA. Output: BED file (4 columns below).

BED output columns:
  1. sequence name (read name, chromosome, or any sequence identifier)
  2. start position
  3. end position
  4. annotation (feature label, e.g. Hsat1, aSat, ct)

## Repositories
- KaryoScope tool:     /Users/ychen/Documents/GitHub/KaryoScope
- Analysis repo:       /Users/ychen/Documents/GitHub/KaryoScope-analysis
  (working directory — all new scripts and outputs go here)

## Current task
1. Build an all-chromosome dendrogram that identifies structural outliers
across the HPRC pangenome centromere dataset. The goal is to find rare
centromere structural variants (deletions, inversions, rearrangements)
validated against FISH data, while avoiding false outliers from normal
population variation.

2. A stacked bar chart showing total number of haplotypes per chromosome.
   Each bar shows Major (blue) and Outlier (red) counts.
   Status: DONE — agent_results/allchr_barplot.svg
   Script: scripts/KS_allchr_barplot.py

3. Update filter.pdf (in KaryoScope repo) with human-readable NucFlag
   QC filter names and heatmap panel.
   Status: DONE — flowchart order: Erroneous -> Collapsed -> Collapsed
   (with Variants). Panel 2 replaced with % retained heatmap.
   Script: KaryoScope/results/figureA/generate_filter_flowchart.R

4. Allele-specific outlier analysis for chr3, chr8, chr11, chr12.
   Pairs h1/h2 by sample to check if outliers are monoallelic or biallelic.
   Status: DONE — agent_results/allchr_allele_heatmap.svg
   Script: scripts/KS_allchr_allele_heatmap.py
   TSVs: allchr_allele_summary.tsv, allchr_allele_pairs.tsv,
          allchr_allele_cooccurrence.tsv
# Data

## Input
- BED: /Users/ychen/Documents/GitHub/KaryoScope/local_data/centromere_region_beds/
        pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed
- Colors: /Users/ychen/Documents/GitHub/KaryoScope/resources/databases/KS_human_CHM13

## Output
- Directory: /Users/ychen/Documents/GitHub/KaryoScope-analysis/agent_results/
- Always export: SVG (primary) + PNG (for self-validation screenshot)

## Excluded chromosomes
Never include these — fewer than 100 haplotypes in HPRC:
chr13, chr14, chr15, chr21, chr22, chrY (all acrocentrics + chrY)

# Pipeline

## Pipeline logic (how the steps connect)

Step 1 and Step 2 each build their OWN feature matrix from the raw BED
file independently. Step 1's matrix is discarded after producing labels.

  Step 1: BED --> [per-chr matrix in RAM] --> cluster labels (TSV)
  Step 2: BED + labels (TSV) --> [new global matrix in RAM] --> dendrogram (SVG)

Step 1 loops chromosome-by-chromosome. Each chromosome gets its own
matrix, its own clustering, its own k selection. The output TSV has
every haplotype labelled Major or Outlier with a divergence score.

Step 2 reads that TSV for labels only (no matrix values carried over).
It picks one representative per cluster (centroid-proximal for Major,
most divergent for Outlier), pools all 18 chromosomes into one global
matrix, and computes a single unified dendrogram.

Why different matrix types between steps:
  - Step 1 uses blockweight (equalises edge vs abundance within a
    single chromosome's feature space)
  - Step 2 uses plain zscore (blockweight across chromosomes would
    distort per-chromosome structural signal; plain zscore standardises
    each transition column without chromosome-specific reweighting)

## Step 1: Clustering (KaryoScope_cluster_analysis.py)

Per-chromosome clustering to determine Major/Outlier labels.
Loops through each chromosome independently: builds a feature matrix
encoding directional transition counts between adjacent repeat blocks,
applies Ward's linkage, selects optimal k via silhouette (k=2-10).
Largest cluster = Major; rest = Outlier.

Required flags — never omit these:

  --analysis-mode structure
  --edges directional
  --no-abundance
  --max-sequence-length 50000000
  --exclude-features "novel"
  --matrix-type count_log1p_zscore_blockweight
  --k-selection silhouette

## Step 2: Dendrogram (KS_allchr_dendrogram.py)

Takes Step 1 labels, selects one representative per cluster, builds a
NEW global matrix from raw BED across all chromosomes, computes a
unified dendrogram. No per-chromosome clustering loop here — the main
dendrogram is one global computation.

Three-stage outlier refinement (applied before dendrogram construction):
  1. Silhouette-optimal k per chromosome (from Step 1 labels)
  2. Silhouette threshold filter — collapse weak splits to k=1
  3. Centroid scan — flag Major members > N SD from centroid (stage 2)

Representative selection logic:
  - Major: haplotype with lowest raw_divergence (most typical)
  - Outlier: haplotype with highest raw_divergence (most extreme)
  - Override with --all-haplotypes to show every haplotype

Recommended parameters:
  --matrix-type count_log1p_zscore
  --sil-threshold 0.5
  --centroid-sd 5

Why this combination:
  - zscore_blockweight clustering (Step 1) produces clean Major/Outlier splits
  - sil-threshold 0.5 suppresses forced splits on chromosomes without
    genuine structure (chr1, chr2, chr4, chr5, chr6, chr7, chr9, chr10,
    chr17, chr20 all get collapsed to k=1)
  - centroid-sd 5 rescues rare individual outliers missed by clustering
    (e.g., chr5 HG00558#1 HSat3 deletion at 10.6 SD)

Why NOT count_log1p alone:
  - Passes all validations but massively over-splits (63% of haplotypes
    become "outliers" — most look visually identical to Major)
  - Many outlier clusters of n=50-100 are normal population subtypes,
    not structural variants

## Step 3: Annotation (KS_allchr_annotate.py)

Compares each outlier representative's BED features against its
chromosome's Major. Reports:
  - Block order changes (rearrangements, inversions, swaps)
  - Feature gains/losses (gained/lost satellite blocks)
  - Major abundance shifts (>15% change only — user preference)
  - Does NOT report subtle abundance changes; these are noise

## Step 4: Generate 3 plot variants

1. Clean (no annotations)
2. Annotated (all outliers labeled)
3. Filtered (subtle "edge pattern difference" outliers removed entirely)

See full commands in logs/session.md.

## Methods text (for manuscript)

Centromeric sequences from all HPRC haplotypes were annotated using
KaryoScope against the CHM13 reference k-mer database. Analysis was
restricted to 18 non-acrocentric autosomes plus chrX (excluding chr13,
chr14, chr15, chr21, chr22, and chrY due to insufficient haplotype
representation). Annotations labelled "novel" were excluded.

Structural outlier detection proceeded in two stages. First, for each
chromosome independently, we constructed a feature matrix encoding
directional transition counts between adjacent repeat blocks (e.g.,
aSat->HSat3). Counts were log-transformed, z-score normalised, and
reweighted by block number to balance edge and abundance contributions.
Ward's method hierarchical clustering was applied and the optimal k was
selected by maximising the silhouette coefficient (k = 2-10). The
largest cluster was designated Major; remaining clusters were designated
Outlier. Second, a silhouette threshold filter was applied: chromosomes
with silhouette score below 0.5 at k = 2 were collapsed to k = 1,
suppressing forced splits on structurally homogeneous chromosomes. A
centroid-distance scan then flagged individual haplotypes within Major
clusters exceeding 5 standard deviations from the centroid as stage-2
outliers, rescuing rare singletons (prevalence <0.3%) that clustering
alone cannot isolate.

To visualise all chromosomes jointly, one representative per cluster was
selected (centroid-proximal for Major, most divergent for Outlier). A
new global feature matrix was constructed from the raw BED annotations
using log-transformed, z-score-normalised directional transition counts
across all representatives, and a unified dendrogram was computed via
Ward's linkage. Each outlier was structurally annotated by comparing its
repeat block order and composition against the Major representative,
reporting rearrangements, gains/losses, and abundance shifts exceeding
15%.

# Scripts

## Existing scripts (READ ONLY — do not modify)
- scripts/KaryoScope_cluster_analysis.py   (per-chromosome clustering)
- scripts/KaryoScope_cluster_plot.py        (per-chromosome plotting)

## New scripts
- scripts/KS_allchr_dendrogram.py   all-chromosome dendrogram with
                                     sil filter + two-stage outlier detection
- scripts/KS_allchr_annotate.py     structural annotation of outliers
- scripts/KS_allchr_barplot.py      stacked bar chart (Major vs Outlier counts)
- scripts/KS_allchr_allele_heatmap.py  allele-specific outlier co-occurrence
                                        (imports filtering from barplot)

# Validation criteria

The output is only correct if ALL THREE of these appear as expected.
Check after every plot is generated.

FISH-validated structural variants (from screenshor4claude/slides.pptx):

1. chr3  — NA21144#1#CM094092.1
   Deletion of HSat1A and second aSat.
   Must appear as an outlier cluster (stage 1).
   Status: PASS → chr3_Outlier_6

2. chr5  — HG00558#1#CM088494.1
   Deletion of HSat3 (1 in 359 haplotypes — 0.3%).
   Must appear as an outlier (stage 2 centroid scan rescues this).
   Status: PASS → chr5_Outlier_S2_4 (10.6 SD from centroid)
   Note: silhouette picks k=2 for chr5 which cannot isolate a singleton.
   The two-stage approach is essential for this case.

3. chr9  — HG02630#1#CM091811.1
   Inversion (paternal, confirmed by trio FISH on parents).
   Active aSat and HSat3 blocks are swapped.
   Status: PASS → chr9_Outlier_S2_1

Additional FISH samples (slides 5-7, not in primary validation):

4. chr12 — HG03816: Second aSat array (ambiguous in FISH)
   Our result: #1 Outlier, #2 Major — matches FISH observation

5. chr17 — NA20799: Second HSat2/3 (ambiguous, too small to measure)
   Our result: both haplotypes in different outlier clusters

6. chr19 — NA18960: Large 3x interrupted aSat (FISH confirmed)
   Our result: stage 2 outlier at 3.0 SD — caught by centroid scan

## Parameter sensitivity (tested grid)

Silhouette threshold x centroid SD:

  sil=0.3 sd=3:  all PASS, 134 reps
  sil=0.3 sd=5:  all PASS, 100 reps
  sil=0.5 sd=3:  all PASS, 103 reps
  sil=0.5 sd=5:  all PASS, 66 reps   ← recommended (cleanest)
  sil=0.5 sd=10: chr5+chr9 FAIL
  no stage 2:    chr5 always FAIL

Without stage 2, chr5 fails at every configuration.
Without sil filter, chromosomes with sil<0.3 produce forced splits.

## Known data discrepancy: chr5 haplotype count (358 vs 359)

The reference PDF (filter.pdf) shows 359 chr5 haplotypes post-QC, but
the clustering pipeline reports 358. The missing sequence is
HG03050#1#CM098762.1, which spans 123.2 Mb (nearly the full chromosome).
The --max-sequence-length 50000000 flag in the clustering script correctly
filters it out. This is expected behavior, not a bug.

# Known cosmetic issues

See screenshot: screenshor4claude/screenshot1.png

1. Bar height too small — FIXED (--row-height 12, --bar-height 10)
2. Gap between dendrogram and bars — FIXED (flush)
3. Tangled dendrogram — FIXED (global distance recomputation)
4. Labels at wrong position — FIXED (between dendrogram and bars)
5. Chromosome ordering — FIXED (chr1→chrX, --dendro-order for raw)
6. Cluster size in labels — FIXED ("chrN Major n=X" / "chrN [sample] n=X")
7. Non-white background under bars — FIXED (pure white, no fills)
8. Annotation text readability — FIXED (short labels, color-coded,
   positioned at bar end, 9px sans-serif)
9. Add a colored block between the labels and bars. see `/Users/ychen/Documents/GitHub/KaryoScope-analysis/screenshor4claude/chromosome_block.png`. each chromosome should have their unique color. avoid similar next to each other. 

# Behaviour rules

## Before starting any task
1. Read existing scripts and any relevant data files first
2. Work autonomously to completion — do not pause for approval mid-task
3. If genuinely blocked (missing file, ambiguous parameter), make the
   most conservative choice, log the decision in logs/session.md,
   and continue. Do not wait for a response.
4. Never ask clarifying questions — infer intent from context and proceed

## After every script run
1. Print sanity check: dimensions, sample counts, any warnings
2. Confirm output files were written and are non-empty
3. Export a PNG screenshot and READ it to validate the output visually
4. Flag anything biologically unexpected before continuing
5. Append a timestamped entry to logs/session.md:
   - What ran
   - What the output was
   - Any unexpected findings or decisions made

## User preferences
- Structural changes (block rearrangements, deletions, gains) matter
  more than subtle abundance shifts
- Only report abundance changes if >15% difference
- The user prefers fewer false positives over catching everything
- The user cares about biological interpretability, not just statistics
- Do not describe "subtle edge differences" as meaningful outliers

## Do not
- Modify KaryoScope_cluster_analysis.py or KaryoScope_cluster_plot.py
- Include chr13, chr14, chr15, chr21, chr22, or chrY in any analysis
- Omit any of the required flags listed above
- git push (unless explicitly asked)
- Create output files outside of agent_results/
- Proceed past a validation failure without flagging it
- Stop mid-task to ask for confirmation or approval
- Wait for human input once a task has started
- Use Unicode arrows or special characters in TSV output — use ASCII only

## On task completion
Write a final entry to logs/session.md with:
- DONE: [timestamp]
- Files produced (paths + sizes)
- Validation results for all 3 FISH-validated samples
- Any decisions made autonomously and why
- Recommended next steps
