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
Build a NEW all-chromosome dendrogram that shows all non-excluded chromosomes
on a single plot. This is distinct from the per-chromosome plots already done.
The new plot should look like a concatenation of per-chromosome clusters,
because each chromosome has its own centromere composition signature.

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

# Required flags — never omit these

These options must be present in every call to KaryoScope_cluster_analysis.py:

  --analysis-mode structure
  --edges directional
  --no-abundance
  --max-sequence-length 50000000
  --exclude-features "novel"

Matrix type (choose one — all three share the log1p base transform):
  count_log1p                        — log(count+1), no further normalization
  count_log1p_zscore                 — log(count+1) + per-column z-score
  count_log1p_zscore_blockweight     — log(count+1) + z-score + block reweighting

The two z-score variants (count_log1p_zscore, count_log1p_zscore_blockweight)
normalize columns to zero mean / unit variance. This can flatten magnitude
differences needed to detect subtle deletions (e.g. chr5 HSat3 deletion).

K-selection method (choose one):
  --k-selection silhouette | --k-selection calinski

Tested all 6 combinations (3 matrix × 2 k-selection). Results:

  count_log1p + silhouette              → ALL PASS (total K=139) ← USE THIS
  count_log1p + calinski                → chr5 FAIL (total K=89)
  count_log1p_zscore + silhouette       → chr5 FAIL (total K=122)
  count_log1p_zscore + calinski         → chr5 FAIL (total K=137)
  count_log1p_zscore_blockweight + sil  → chr5 FAIL (total K=122)
  count_log1p_zscore_blockweight + cal  → chr5 FAIL (total K=137)

Only count_log1p + silhouette passes all 3 validation criteria.
Z-score normalization flattens magnitude differences needed to detect
the chr5 HSat3 deletion.

# Scripts

## Existing scripts (READ ONLY — do not modify)
- scripts/KaryoScope_cluster_analysis.py   (per-chromosome clustering)
- scripts/KaryoScope_cluster_plot.py        (per-chromosome plotting)

## New scripts to write
Write ALL new code to new files. Never edit the existing scripts above.
Suggested names (not fixed):
- scripts/KS_allchr_dendrogram.py     all-chromosome dendrogram builder
- scripts/KS_allchr_plot.py           all-chromosome plot renderer

Before writing: read the existing scripts to understand data structures,
matrix formats, and plotting conventions. Mirror their style.

# Validation criteria

The output is only correct if ALL THREE of these appear as expected.
Check after every plot is generated:

1. chr3  — NA21144#1#CM094092.1
   Deletion of Hsat1 and second aSat. Must appear as an outlier cluster.
   Reference hap2: NA21144#2#CM094108.1
   Status: PASS with count_log1p + silhouette → chr3_Outlier_8

2. chr5  — HG00558#1#CM088494.1
   Deletion of HSat3. Must appear as an outlier (any outlier with similar
   composition is acceptable — does not need to be a solo cluster).
   Reference hap2: HG00558#2#CM088523.1
   Status: PASS with count_log1p + silhouette → chr5_Outlier_5
   Note: FAILS with all z-score variants and with calinski. Only
   count_log1p + silhouette detects this deletion (see grid above).

3. chr9  — HG02630#1#CM091811.1
   Inversion. Must appear as an outlier cluster.
   Status: PASS with count_log1p + silhouette → chr9_Outlier_4

All 3 pass with count_log1p + silhouette. If a future change breaks
validation, do NOT switch to z-score variants — they cannot detect chr5.

# Known cosmetic issues to fix

See screenshot: screenshor4claude/screenshot1.png

1. Bar height too small — target: at least 8px per haplotype row
   Status: FIXED (--row-height 12, --bar-height 10)
2. Unexplained gap between dendrogram and bar panel — must be zero gap;
   both panels should share the same y-axis with no whitespace between them
   Status: FIXED (dendrogram tips flush with bars)
3. Tangled dendrogram — likely caused by per-chromosome distance calculation
   being reused without recalculation. Recalculate distance globally across
   all chromosomes before drawing the all-chromosome dendrogram.
   Status: FIXED (global distance matrix recomputed)
4. Chromosome labels should appear BETWEEN the dendrogram and the bar panel,
   not at the right end of the bars.
   Status: FIXED (labels in column between dendrogram and bars)
5. Chromosomes should be ordered chr1 → chrX in the dendrogram. Use a
   constrained ordering that groups chromosomes sequentially (chr1, chr2, …).
   It is acceptable if a few outliers slightly break the order.
   Status: FIXED (default: chr1→chrX; use --dendro-order for raw tree order)
6. Labels must show cluster size: "chrN Major n=X" or "chrN [sample] n=X".
   The representative is 1 of n haplotypes in that cluster.
   All n values for a chromosome must sum to the total haplotype count.
   Status: FIXED

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

## Subagent routing
Parallel dispatch (all conditions must be met):
  - Tasks touch different files with no shared state
  - Example: running the clustering script while reading existing code

Sequential dispatch (any condition applies):
  - Next step reads output of current step
  - Shared files or data structures involved

## Do not
- Modify KaryoScope_cluster_analysis.py or KaryoScope_cluster_plot.py
- Include chr13, chr14, chr15, chr21, chr22, or chrY in any analysis
- Omit any of the six required flags listed above
- git push
- Create output files outside of agent_results/
- Proceed past a validation failure without flagging it to me
- Stop mid-task to ask for confirmation or approval
- Wait for human input once a task has started

## On task completion
Write a final entry to logs/session.md with:
- DONE: [timestamp]
- Files produced (paths + sizes)
- Validation results for all 3 FISH-validated samples
- Any decisions made autonomously and why
- Recommended next steps