# Rearrangement detection — model, design, and Engine A spec

> Status: **design agreed; Engine A fully implemented** — measurement
> (`core/colocalization.py`), differential test (`core/rearrangement.py`), and the
> `detect-rearrangements` CLI. The statistics are a **v1 for coauthor review** (the
> **Group A** bucket in `OPEN_QUESTIONS.md`): CMH across length buckets + BH-FDR +
> recurrence/effect/floor gates, with read independence assumed (surfaced as a runtime
> warning, not enforced — open item). **Engine B** (OLC clustering) is **fully implemented**
> — aligner (`core/feature_align.py`), overlap graph + clustering + seed-anchored consensus
> (`core/feature_assembly.py`), and the `cluster` CLI. Rendering of the cluster layouts is
> deferred to the (separately gated) plotting tier.

This document is the mental model and design for detecting **recurrent rearrangements**
from KaryoScope-annotated long/short reads, superseding the old `cluster_analysis`
approach. It is the reference for the `colocalization` / future `detect-rearrangements`
and clustering tools.

## 1. What we are trying to detect

Recurrent **rearrangements**: feature arrangements present in a sample that the normal
(unrearranged) genome does not produce — e.g. a telomere abutting sequence it is never
near normally, or rDNA brought close to α-satellite. "Recurrent" = supported by multiple
independent reads. We want to (a) call these events, (b) compare a sample against a
**matched control** and against a **normal reference** with statistical rigour, and
(c) organize the carrier reads into structural groups with a consensus and a picture.

## 2. Generative mental model

Each latent cluster *c* has a **reference annotated sequence** `R_c`: an ordered list of
`(feature, length)` segments (a string over the feature alphabet, with a length per
symbol). A read is a **noisy, partial observation** of some `R_c`:

1. **Windowing** — a sub-interval `[a, b)` of the reference is observed; `b−a` ranges over
   orders of magnitude (short reads → ultra-long ONT). Two reads from the same cluster may
   overlap only partially, or not at all.
2. **Local edits** — segment length perturbation (same feature, jittered length), segment
   insertions/deletions, and *(rare, deferred)* inversion of a contiguous block.
3. **Annotation noise** — mislabeled or `novel` stretches, boundary jitter, specks.

**Reads are not oriented or anchored.** Telogator emits reads with the telomere at either
end or internal, and many reads have no telomere at all. So there is no common coordinate,
and any method must be **orientation-invariant** — a read and its reverse are the same
object. (Past methods failed exactly here, emitting mirror-image duplicate clusters.)

## 3. The framing this forces

A read is an **ordered sequence of feature-segments**, not a point in feature space.
Composition / fraction vectors (the old approach) are **blind to arrangement** and, worse,
**break on partial overlap**: two reads from different parts of the same reference have
different compositions and look dissimilar. Order matters on a spectrum — k=1 (composition)
< k=2 (adjacency, the old interspersion) < k-grams < full alignment — and only the
order-aware end handles offset, indels, and partial overlap.

Crucially, the annotation has already done a nucleotide-level alignment-free reduction, so
working at the **feature/annotation level** (tens–hundreds of segments) makes sequence
methods both cheap and semantically right: we compare *structure*, and the annotation *is*
the structure.

## 4. Two engines

- **Engine A — colocalization detection (alignment-free).** The primary biological signal.
  Per read, measure feature proximities; flag and test abnormal ones. No alignment;
  orientation-invariant; robust to read length and partial overlap.
- **Engine B — alignment-based clustering + consensus + visualization.** Groups carrier
  reads into structural haplotypes for interpretation. Orientation handled in the metric.

The expensive machinery (alignment) lives only in Engine B, on already-subset reads.

## 5. Colocalization, not just adjacency

The signal is **abnormal colocalization**: two features unusually *close on the same
molecule*, not necessarily adjacent (rDNA and α-sat may have a few segments between them).
Adjacency is the distance-0 special case (`adjacency ⊂ colocalization`).

Why colocalization rather than the junction adjacency itself: a rearrangement is a physical
**breakpoint**, but the feature pair *at* the breakpoint is often non-distinguishing
(`arm|arm`) or cryptic, so the abnormal adjacency is **masked**. The colocalized pair
`{A, B}` *is* distinguishing by construction, so it is detectable even when the exact
junction is not. Colocalization is also **orientation-invariant for free** (proximity is
symmetric), and a same-read colocalization is **direct single-molecule evidence of physical
linkage** — the thing long reads uniquely provide. The proximity window also *is* the
local-vs-global dial: how far apart `A` and `B` may be and still "count."

### Two orthogonal axes of abnormality

1. **Distance/proximity** — `A` and `B` closer than the reference ever shows (novel
   juxtaposition).
2. **Abundance/frequency** — `A` and `B` colocalize at a *normal* distance, but the
   *fraction of molecules* showing it is off (copy-number / dosage).

Both are questions about one quantity — the per-pair **colocalization rate** — at two
window settings (see §7).

## 6. Engine A — measurement layer (implemented)

Input is the **overlay-annotation output** (one coalesced feature-segment BED per read; the
feature resolution / precedence is chosen upstream by `overlay-annotations`).

Per read, for every unordered pair of distinct feature types both present, compute the
**minimum bp gap** between the nearest interval of each (0 if adjacent). A single
left-to-right sweep keeping the most-recent end per feature gives all pair gaps in
`O(segments × distinct-features)`. **bp**, not segment count: bp is physical, invariant to
annotation granularity, and directly comparable to reference genomic distances. A
**`min_occurrence_bp`** denoise treats a feature as absent unless it has an interval at
least that long, so annotation specks don't anchor phantom proximities (the gap is still
measured in true bp, including any filtered specks in between). Pairs are canonical
unordered tuples → orientation-invariant.

This layer is pure and non-statistical; everything below builds on it.

## 7. Engine A — differential test (specified; deferred)

The per-pair, per-sample statistic is the **colocalization rate**

> `r = (reads with min-gap(A,B) ≤ W) / (all reads in the length bucket)`

- **Denominator = all reads in the length bucket.** The all-reads denominator makes `r` an
  *abundance*, so dosage changes show up; conditioning on co-presence would divide that
  signal out. Reads are bucketed by length (short / long, boundaries TBD) and compared
  within bucket; combine buckets with a **stratified test (Cochran–Mantel–Haenszel)** — never
  pool a long-read rate against a short-read denominator.
- **Two window settings, one machinery:** tight `W` (per-pair reference percentile) →
  *distance-novelty*; generous `W` (read-scale, value TBD — report a couple for now) →
  *abundance*. Gaps are measured once, so `r` at several `W` is free.
- **Test:** per candidate pair, a proportion test (Fisher / beta-binomial for
  overdispersion) on exp vs control; **FDR** across the candidate set; gated on a minimum
  **recurrence** and **effect size**, not p-value alone.

### Baselines (they answer different questions)

- **Matched control sample** — the clean differential baseline (reads vs reads, apples to
  apples); the headline `exp vs control` rate test.
- **Annotated CHM13 reads** (real ONT/HiFi/short reads run through KaryoScope) — a
  read-matched *normal* baseline. Three jobs: (1) abundance baseline, (2) an **empirical
  false-positive floor** — a normal genome's reads should show ~0 abnormal colocalizations,
  so whatever rate they do show bounds the artifact rate any real call must clear, (3) with
  the **CHM13 genome** (ground-truth distances), defines the per-pair "abnormally close"
  thresholds.

### Candidate sets

- **reference-abnormal** pairs (rare/never close in the reference) → distance-novelty
  (default).
- **reference-normal** pairs → abundance (where dosage signal lives).
- full set → exploratory option.

### Anti-circularity discipline (Group A)

Thresholds and the candidate set come from the **reference**, fixed **before** looking at
exp/control counts. We never use the tested reads to both nominate and test, never let the
control define its own comparison threshold, and require read independence (dedup) so
recurrence counts molecules, not duplicates.

## 8. Engine B — overlap-layout-consensus clustering

Engine B is **mini-assembly in feature space**: a read is a window of a latent reference
sequence (in unknown orientation), so grouping reads + building a consensus *is* the
overlap-layout-consensus (OLC) problem — on tens–hundreds of feature segments over a
~tens-symbol alphabet, not 10⁶ bases. This resolves partial overlap (transitive linking),
orientation, consensus, and the eventual visualization with one model. Reads are always
pre-subset, so an `O(N²)` pairwise pass is fine.

**Input.** Whatever overlay-annotations BED subset the user hands it (decoupled from Engine
A; one natural input is the carriers of an Engine-A call). Clusters the reads as given, with
an optional `--min-length` filter; no internal length buckets. In practice only **long
reads** are clustered (short reads carry too little structure — they still feed Engine A's
colocalization rates). The OLC graph would correctly absorb a contained short read if mixed
in, so the tool *can* handle short/long/combined; long-only is the default.

### Overlap — local alignment of feature-segment sequences (`core/feature_align.py`)

Each read is an ordered `(feature, length)` sequence (the coalesced overlay). A
**Smith-Waterman local alignment** over segments, with:

* **Substitution score from the hierarchy** (tiered): exact match best; same coarse group
  (`satellite`/`arm`/`ct`/telomere-types — `HSat1A`↔`HSat1B`) partial; `novel` neutral (0);
  otherwise mismatch. Supplied as a `sub_score(a, b)` callable, so the aligner is
  hierarchy-agnostic and trivially testable.
* **Match reward weighted by `min(len_a, len_b)`** — bp-weighted, so large shared blocks
  dominate and specks contribute ~nothing; and **length-change-lenient**, since the
  `|len_a − len_b|` difference is neither rewarded nor penalized (a feature's length may
  drift). The `min_occurrence_bp` denoise also applies.
* **Linear per-bp gap** — skipping a segment costs `gap_factor × its length`, so skipping a
  large structural block (a real insertion/deletion) is penalised in proportion, while a
  feature's length drift is absorbed by the match rule above. (Affine is a later refinement.)
* **Best of forward / reverse** B — a read and its reverse therefore have ~0 distance, which
  kills the mirror-image-cluster failure by construction.

The alignment returns its score, the aligned segment-index columns, and the aligned spans,
which classify the overlap as **dovetail**, **containment**, or **internal**.

### Layout — the overlap graph (`core/feature_assembly.py`)

* **Edges = *proper* overlaps only** (dovetail or containment that clear a minimum overlap
  length and a minimum normalized concordance). Internal-only matches are rejected — they are
  usually shared repeats and are the classic cause of false chaining. This is what makes
  connected-components clustering trustworthy.
* **Clusters = connected components** (v1; community detection is a later knob if chaining
  appears). **Singletons are kept** as size-1 clusters so nothing vanishes silently.
* **Orientation propagation**: each edge carries a relative orientation; BFS from the cluster
  seed assigns each member an orientation (signed/2-coloring), flagging rare conflicts.
* **Bridges = rearrangements.** A read joining two otherwise-separate clusters via proper
  overlaps on different structures is a junction; a **recurrent** bridge is a recurrent
  rearrangement — i.e. Engine A's signal as graph topology, with recurrence separating a real
  rearrangement bridge from a one-off chimera.

### Consensus

**Seed-anchored** (v1): the longest read is the backbone; members are aligned/oriented to it
and the per-position **majority feature** (with per-position support) is the cluster's
consensus structural sequence, extended where members reach past the seed. Members linked
only transitively (not overlapping the seed) don't inform the seed-frame consensus in v1;
progressive layout is the v2 fix.

### Output (data only; rendering deferred)

A clusters table (id, members, sizes), one consensus feature-BED per cluster (with support),
and a member layout TSV (orientation + offset). The plotting tier renders it later (it's
separately gated on the `karyoplot` push-down). **Cluster-level differential testing stays
deferred** — a cluster is a modeled object, so testing it risks circularity.

### Scale

Pre-subset + `min_occurrence_bp` denoise (fewer segments) + a **feature-content prefilter**
(skip pairs whose feature sets barely overlap) prune the `O(N²)` before the DP; pure Python
suffices at subset scale, with minimizer-style seeding as the escape hatch.

## 9. Open knobs / deferred

- generous `W` scale (report a couple for now); length-bucket boundaries; `min_occurrence_bp`.
- explicit copy-number model — **rejected for now** in favour of the rate comparison.
- distance-shift test (compare A–B *distance distributions* exp vs control) — optional refinement.
- self-pairs (`A` close to another `A`, for tandem amplification) — possible extension.
- inversions — treated as rare; revisit if data demands.

## 10. Validation on v2 data (findings)

Both engines were run on the committed `KS_human_CHM13_v2` overlays (priority preset,
`region`/`repeat`/`subtelomeric`).

**Engine A** — `detect-rearrangements` U2OS (ALT) vs IMR90 (primary), length-stratified at
50 kb, ran on the full overlays in ~6 s. Output is well-formed (CMH + BH-FDR); the top
signal is satellite↔telomere / satellite↔repeat colocalization differences (plausible for
ALT). Two caveats confirmed in practice: (a) with no `--reference` the artifact floor is 0,
so `reference_abnormal` is trivially true everywhere and the floor gate is inactive — a real
normal baseline (annotated CHM13 reads) is needed before calls are trustworthy; (b) at this
read count q-values are astronomically small for almost everything, so the **effect-size and
floor gates**, not the p-value, are what actually discriminate. (Dedup warning fired, as
designed.)

**Engine B** — `cluster` on a 250 long-read U2OS subset ran end-to-end (aligner → graph →
clusters → consensus; orientation conflicts flagged; singletons kept) but **over-merged**:
223/250 reads in one cluster, whose consensus is dominated by interspersed repeats
(LINE/SINE/LTR) at high support. The cause is the chaining the design anticipated — long
reads share *ubiquitous* interspersed-repeat content, so they form genuine proper overlaps
through generic repeats. Tightening the edge criteria helps but doesn't resolve it
(overlap ≥ 30 kb, identity ≥ 0.9 → largest cluster 223 → 91; 110 clusters, 20 multi-read).
**v2 fix (needed before Engine B is usable on real data):** down-weight / mask ubiquitous
repeat-class features in the aligner so overlaps must rest on *distinctive* structure, and/or
replace connected components with **community detection** (modularity) to resist chaining.
Runtime was ~50–60 s for 250 reads in pure Python — fine for a subset, but the Jaccard
prefilter pruned little (long reads share many features), so larger subsets will want the
repeat down-weighting (also fewer edges) or minimizer-style seeding.

**v2 fix implemented — feature weighting (`--weight-method`).** Per-feature weights scale
both the match reward and the overlap-length criterion, so an overlap must rest on
*distinctive* shared content. Two builders: `repeat-mask` (default) zeroes the hierarchy's
`Interspersed_Repeat` subtree (LINE/SINE/LTR/DNA/…) plus `nonrepeat`; `idf` down-weights by
read frequency. Findings on the U2OS subset:

- **`idf` (gentle, linear `1 − df/N` + 0.1 floor) did *not* work** — long reads are so
  repeat-dominated that even floored repeat content × huge overlap length clears the
  threshold (largest cluster 223 → **244**, *worse*). Frequency also can't separate the
  drivers cleanly here: `SINE` (0.60) and `canonical_telomere` (0.54) have nearly equal
  prevalence, so any frequency cutoff that removes SINE also removes telomere.
- **`repeat-mask` works for the repeat-chaining** — masking by *biology* (the interspersed-
  repeat class) cleanly keeps telomere (not in the `repeat` featureset) while zeroing
  LINE/SINE/LTR/DNA/nonrepeat: largest cluster 223 → **127**, singletons 11 → **112**.
- The **residual 127-read cluster is structurally coherent**, not repeat chaining: it is
  held together by near-universal *structural* features in this telomeric read set
  (`canonical_telomere`/`TAR1`). Subdividing it into haplotypes is the **community-detection**
  job — the remaining v2 lever (modularity on the within-cluster overlap graph), still open.

**The bigger lever: chromosome identity (`chromosome-telomere-satellite` preset + layer-aware
scorer).** The deeper issue is that *structural features alone can't tell chromosome ends
apart* — reads from different chromosomes look alike (telomere → subtelomere → satellite →
arm). Under `telomere-satellite` (no `repeat`), the U2OS subset chained almost completely
(**249/250 in one cluster**, and `idf` didn't help — the ubiquitous feature is now `arm`).
Adding the chromosome featureset as a **two-layer `chromosome:structural` label** (new preset,
emitting `chr13:canonical_telomere` etc.) is the fix:

- The aligner's substitution scorer is **layer-aware** (`chromosome_aware_substitution`):
  the structural layer is hierarchy-graded as before; the chromosome layer adds a **soft
  per-bp penalty unless the two positions are the same *specific* chromosome** (`chr…`).
  Different chromosomes *and* ambiguous labels (`autosome`/`categorized`/…) are penalized, so
  reads don't chain through unresolved assignments; the penalty is finite, so a translocation
  read still bridges the two chromosomes over its matching halves.
- On the U2OS subset: the 249-blob breaks into **clean per-chromosome clusters** (largest 18;
  e.g. chr4×18, chr18×10, chr7×8, chr20×7), with the **multi-chromosome clusters surfacing as
  candidate translocations** (chr4+chr5, chr14+chr18, …). Runtime ~3 s for 250 reads.
- Treating ambiguous chromosome labels as *neutral* (the first try) re-chained them into an
  `autosome` blob (largest 84); the **strict** rule above resolved it. (Grading ambiguous
  labels through the chromosome hierarchy — `chr5` is consistent with its parent `autosome` —
  is a possible future refinement.)

## 11. Pre-overlay denoise — hierarchical binning (`bin-annotations`)

The raw featureset BEDs are extremely fragmented: telogator's per-base annotation shreds a
single biological region into a chaotic alternation of tiny segments (e.g. an rDNA array
read as hundreds of interleaved 1–93 bp `rDNA`/`ct`/`novel` slivers). That fragmentation
propagates through the overlay into the feature-segment sequences Engine B aligns and the
read plots, making both noisy. `bin-annotations` is a **rolling-window mode filter** that
runs on each featureset BED *before* `overlay-annotations` to collapse this noise while
preserving the C4 partition and the sequence length.

**Why a hierarchical vote, not a flat mode.** A flat mode is fragile when related siblings
split their vote: a window of `aSat` 20 / `bSat` 40 / `arm` 40 is 60 % satellite, yet a flat
count hands it to `arm`. So the per-feature window bp propagate **up the database tree**, and
the call descends from the root into the dominant child while its subtree holds a majority,
stopping at the deepest node that still does (`core/annotation_binning.py`). The same window
descends `categorized → centromeric` (60 of 100) and — with the default node-relative
threshold — on to `bSat` (40 of the **60** at `centromeric`); when subtypes split evenly it
honestly stops at the `centromeric` ancestor instead of guessing a leaf.

- **Knobs.** `--majority-fraction` (τ, default 0.5); `--threshold-scope` = `node` (τ relative
  to the current node — the conditional majority, more specific, the default) or `window` (τ
  relative to the whole window — conservative, climbs to internal nodes more readily). The
  node scope can descend to a leaf holding < τ of the whole window when the top-level split is
  near-even; that is an inherently ambiguous boundary, and τ tunes it.
- **Boundaries stay sharp.** If no top-level group has a majority (e.g. a clean 50/50 boundary
  between two unrelated features), the descent falls back to flat plurality (ties toward the
  deeper label) rather than emitting the generic root — so binning denoises interiors without
  smearing a vague label across every boundary.
- **`novel`** votes as a top-level leaf (a novel-dominated window stays `novel`); every other
  out-of-tree label is the C2 error.

**Validation (U2OS v2, `region`).** Strong fragmentation reduction: at window 101 / node /
τ=0.5, **1,375,812 → 218,558 intervals (6.3×)**, the worst read **6,812 → 896 (7.6×)**, its
shredded rDNA/ct alternation resolved into coherent multi-kb `rDNA`/`ct` blocks; at window
1001 / τ=0, **→ 29,891 intervals (46×)**. Internal-node calls (`centromeric`, `alpha_hor`,
`arm`) appear where subtypes/arms abut at τ>0; **τ=0 always descends to a specific leaf**
(hierarchy-aware plurality), which is what we want for the chromosome layer (always a specific
`chr*`, never an ambiguous `autosome`/`categorized`).

**O(intervals), window-independent.** Between the O(intervals) breakpoints where the window's
entering/leaving base crosses an interval edge, the per-feature counts are linear in the step
offset, so the descent is evaluated O(1) times per segment via `_descent_run`, which returns a
*conservative* lower bound on how many forward steps the call stays constant (recompute-and-
merge keeps the output exact regardless of how tight the bound is). Cost no longer grows with
the window: ~23–28 s for the full `region` BED at *either* window 101 or 1001, vs ~160 s for
the original per-base version. A per-base reference (`bin_intervals_naive`) pins the fast path
in a 600-case property test (τ ∈ {0 … 1}, both scopes), and the fast output is byte-identical
to the original implementation on the real U2OS `region` BED.

## 12. Feature weighting from the reference genome (`genome-weights`, `--weight-method genome-freq`)

The over-merging in §10/the arm-star problem is, at root, clustering on a *ubiquitous* feature:
the chromosome **arm** is shared by every read of a chromosome, so sharing it carries little
information. Read-frequency (`idf`) can't fix this — long telomere reads are uniformly
arm-dominated. The principled signal is **genome-wide coverage**: `core/genome_weights.py`
tallies each feature's bp in the annotated CHM13 reference (one C4 BED per featureset), takes
its fraction `p` of its featureset's partition, and uses information content `-ln(p)` scaled to
`(0, 1]` by the rarest feature across all featuresets. So:

| feature | genome fraction | weight |
| --- | --- | --- |
| `q_arm` | 56.8 % | 0.027 |
| `p_arm` | 26.3 % | 0.064 |
| `ct` | 8.5 % | 0.119 |
| `aSat` | 0.015 % | 0.423 |
| `canonical_telomere` | 0.0028 % | 0.505 |
| `nonsubtelomeric` | 99.8 % | 0.0001 |

This crushes the arm while keeping `canonical_telomere` informative (genome-rare) — telomere is
*not* filler, contrary to a first guess; the genome decides. Because τ=0 binning emits the
**leaves** (`p_arm`/`q_arm`, not the internal `arm`), the down-weighting lands on the labels we
actually produce. A latent bug is fixed alongside: `_weight`/`idf_weights` now key on the
**structural layer** of `chromosome:structural` labels, so weighting applies on composite
overlays at all (previously `repeat-mask`/`idf` silently no-op'd there).

**Validation (U2OS, w1001/τ0 binned, chr-tel-sat overlay, 49-read subset, `--weight-method
genome-freq`).** Clustering ran in ~11 s (vs ~460 s on the w101 binned overlay — w1001 cuts the
per-read segment count that drives the O(N²·L²) aligner). The multi-read clusters are now
chromosome-coherent and surface **recurrent candidate translocations**: a chr4-centric cluster
with 3 reads showing chr4+chr22, and a 2-read chr12+chr9 cluster — where before (uniform
weights) reads were glued by a shared `p_arm` block into a star whose consensus just echoed the
seed. **Distinctive-overlap edge criterion (the `--min-overlap-bp` size gate, measured on
distinctive bp).** Genome-freq weighting alone didn't fully kill arm-chaining: `p_arm` weight
0.064 × a shared 150 kb arm block ≈ 9,600 weighted bp, still above any floor, so pure-arm reads
still formed edges. The fix is to measure the overlap-size gate on matched **distinctive** features (weight ≥
`--distinctive-weight`, default 0.15 — above `arm`/`ct`/`p_arm`/`q_arm`, below the satellites and
telomere); an overlap built only of filler scores 0 here and is rejected. Crucially the bar must
be **small** (~1 kb): a real same-chromosome-arm translocation pair (cluster_2, chr12 q_arm +
chr9 p_arm) shares only ~2 kb of distinctive `mon`/`ITS`, whereas a pure-arm interloper shares 0.

Validated (U2OS w1001/τ0 subset, genome-freq, `--min-overlap-bp 1000` on distinctive bp): the chr4 arm-star
dissolves to singletons, and the all-`p_arm` chr5 read (which had no telomere) drops out of the
chr5 cluster, leaving three structurally-justified pairs — chr12+chr9 (translocation), chr20,
and chr5+`canonical_telomere`. This matches a read-by-read maintainer review of the plots.

Residual limitations: connected-components could still in principle link reads through a *rare*
shared feature — addressed by **community detection** (§13). The layout/consensus weakness (single
seed-relative offset → drift; seed-only consensus) is fixed by the **consensus-coordinate layout**
(§14).

### Recomputing / extending the weights

`data/chm13v2_feature_weights.tsv` is **committed**, so clustering never needs to recompute it. Re-run
`genome-weights` only when the reference changes or to **add a featureset** — the command is generic
(one `--bed FEATURESET=PATH` per featureset; a single O(intervals) streaming pass). The reference is
the annotated CHM13 BEDs `data/raw_bed/chm13v2.0.KS_human_CHM13_v2.<featureset>.smoothed.bed.gz`
(gitignored, ~475 MB — needed *only* here).

```bash
DB=/path/to/KS_human_CHM13_v2                       # hierarchy.tsv (C2 feature validation)
REF=data/raw_bed/chm13v2.0.KS_human_CHM13_v2        # annotated reference BED prefix
karyoscope-analysis genome-weights \
  --bed region=$REF.region.smoothed.bed.gz \
  --bed subtelomeric=$REF.subtelomeric.smoothed.bed.gz \
  --bed chromosome=$REF.chromosome.smoothed.bed.gz \
  --bed repeat=$REF.repeat.smoothed.bed.gz \
  --bed gene=$REF.gene.smoothed.bed.gz \
  --bed acrocentric=$REF.acrocentric.smoothed.bed.gz \
  --hierarchy $DB/hierarchy.tsv \
  -o data/chm13v2_feature_weights.tsv
```

Wrapped as `scripts/compute_genome_weights.sh --db DB [--ref-prefix P] [--featuresets "a b c"]
[-o OUT]`. To add a featureset, give it another `--bed` (or extend `--featuresets`): it is tallied,
normalized on the same `(0, 1]` scale as the rest, and written out. Two requirements — the featureset
must validate against `hierarchy.tsv`, and its reference BED must cover the **whole genome** (the
weight uses the feature's fraction of the featureset total, so a partial reference inflates every `p`).

## 14. Consensus-coordinate layout (`consensus_layout`)

The v1 layout placed each read by a single seed-relative offset, so length-mismatched segments
drifted (a chr20 cluster's shared `HSat3` ended up ~20 kb apart between two reads), and the
consensus only spanned the seed. `feature_assembly.consensus_layout` replaces both: each member's
read→seed **alignment columns** become anchors `(member_bp, consensus_bp)`, and every segment is
mapped through a piecewise-linear function (slope-1 outside the anchors). So aligned features land
on their seed counterpart's coordinate — **they stack vertically** — and a read's overhang
extrapolates beyond the seed, **extending the frame to the union** of all reads. The consensus is
re-derived by majority vote over that union grid. Validated: in the chr20 cluster both reads' `HSat3`
now sits at identical consensus coordinates (415485–424089), their `p_arm` blocks line up, and a
member's telomere overhang extends the consensus past the seed. `layout.tsv` is now per-segment
(consensus coords); `cluster-plot` draws straight from it (no overlay) with length filtering off so
gaps are unambiguous.

## 13. Scaling `cluster` to a whole sample (blocking + memoization + parallelism)

The OLC aligner is all-vs-all: 49 reads (~1.2k pairs) ran in 11 s, but 4005 reads is ~8M pairs,
each a full Smith-Waterman, and was intractable (43 min, unfinished). Three exact (non-lossy)
levers fixed it:

1. **Blocking index (`--block-min-bp`).** Reads are bucketed by the **specific-chromosome leaf**
   of their composite labels (`chr5:p_arm` → `chr5`), summing bp per chromosome; only reads
   sharing a chromosome with ≥ `--block-min-bp` are aligned. A `min_overlap_bp` edge needs
   substantial same-specific-chromosome content (the layer-aware scorer penalizes cross-chromosome),
   so a non-candidate pair can't form an edge. **Translocation reads** are indexed under *every*
   chromosome they span (a chr12+chr9 read is in both buckets), so the translocation signal is
   kept. Bucketing on the chromosome leaf (not the full `chr:struct` composite) is deliberate —
   composite buckets would miss a pair that clusters through a *different* shared element (one read
   `mon`, another `aSat` at the homologous spot). Measured on U2OS: 8M → ~0.9M candidate pairs at
   `--block-min-bp 2000`; ambiguous composites (`categorized:*`, `autosome:*`) are excluded as keys.
2. **Scorer memoization.** The substitution scorer is invoked once per DP cell but has only
   ~hundreds of distinct `(feature, feature)` argument pairs across a run; caching it removes the
   per-cell hierarchy/chromosome recomputation.
3. **Process parallelism (`--workers`/`-j`).** Candidate-pair alignments are independent (edges
   merge afterward in the union-find), so they fan out across processes via `fork` — each worker
   inherits `reads` + the scorer copy-on-write (only pair chunks and edge results cross the
   boundary). Output is identical to serial (edges sorted by `(a, b)`); a unit test pins this.

**Result:** full U2OS (4005 reads, w1001/τ0 overlay) clusters in **~100 s on 8 cores** (≈6.4×;
648 s CPU), into 2509 clusters (40 multi-read) — vs >8 min serial and 43 min unoptimized. More
cores scale further. (`fork` is required, so this is best on Linux/HPC; it also works on macOS for
this pure-Python path.) If a single chromosome bucket is ever too big, a similarity-preserving
sub-bucket (MinHash/LSH on the feature multiset) is the next lever — kept in reserve as it is
*lossy*, unlike the three above.

### Community detection (`--communities`, default on)

Clustering all 4005 reads (connected components) produced a **1058-read mega-cluster spanning
every chromosome** — chaining returned at scale. Diagnosis: it was *transitive bridge-chaining*,
not one runaway chromosome — 708 single-chromosome reads glued by 292 two-chromosome reads and
**52 reads spanning ≥3 chromosomes**. The ≥3-chromosome "hubs" are **noise**, not complex
rearrangements: their chromosome layer is scattered small slivers (biggest block ~18 kb, most
3–7 kb, 8–38 runs/read, ambiguous stretches mixed in), unlike the large clean blocks (~150 kb +
~100 kb) of a real chr4+chr22 translocation read. Blanket-dropping multi-chromosome reads would
risk genuine signal, so the fix is structure-aware: **weighted label propagation**
(`_label_propagation`) subdivides each connected component — a read joins the community it shares
the most overlap weight with, so a sparse bridge attaches to one group instead of merging all.
Orientation parities still come from the component-wide union-find (consistent within each
sub-community). Result on full U2OS: the mega-cluster dissolves into clean per-chromosome groups
(several chromosomes resolve into *multiple* haplotype communities — chr18 → 148 + 83, chr20 →
105 + 48) and 24 recurrent translocation candidates (chr4+chr22 ×43, chr18+chr19 ×34, chr13+chr11
×33, chr1+chr21 ×23, …); 4005 reads → 2544 clusters (75 multi-read), no mega-cluster. (Louvain is
the upgrade if label propagation ever proves too coarse/unstable; it was sufficient here.)

### Runbook — whole-sample clustering (e.g. on HPC)

Inputs needed (the package + `data/chm13v2_feature_weights.tsv` are committed; the large
`data/raw_bed/*.bed.gz` telogator featureset BEDs and the DB `hierarchy.tsv` are **gitignored**,
so copy them to the machine). The genome weights are committed, so the 475 MB CHM13 reference BEDs
are only needed to *recompute* weights, not to cluster.

```bash
DB=/path/to/KS_human_CHM13_v2            # hierarchy.tsv lives here
S=U2OS                                   # sample
RAW=data/raw_bed/$S.telogator.1.KS_human_CHM13_v2

# 1. bin each featureset (window 1001, tau 0 = always descend to a specific leaf)
for fs in region subtelomeric chromosome; do
  karyoscope-analysis bin-annotations --input $RAW.$fs.smoothed.features.bed.gz \
    --hierarchy $DB/hierarchy.tsv --feature-set $fs \
    --window 1001 --majority-fraction 0 -o $S.$fs.binned.bed.gz
done

# 2. overlay with the chromosome-telomere-satellite preset
karyoscope-analysis overlay-annotations \
  --bed chromosome=$S.chromosome.binned.bed.gz --bed region=$S.region.binned.bed.gz \
  --bed subtelomeric=$S.subtelomeric.binned.bed.gz \
  --hierarchy $DB/hierarchy.tsv --preset chromosome-telomere-satellite -o $S.overlay.bed

# 3. cluster ALL reads (scaled): genome-freq weights + distinctive + blocking + parallel
karyoscope-analysis cluster --input $S.overlay.bed --hierarchy $DB/hierarchy.tsv --min-length 0 \
  --weight-method genome-freq --genome-weights data/chm13v2_feature_weights.tsv \
  --min-overlap-bp 1000 --block-min-bp 2000 --workers <NCORES> -o $S.clusters.tsv

# 4. plot the multi-read clusters (cluster-plot draws straight from the layout/consensus;
#    --no-consensus-track omits the union-consensus top row, as in the published figures)
karyoscope-analysis cluster-plot --layout $S.clusters.layout.tsv \
  --consensus $S.clusters.consensus.bed --colors $DB/colors.tsv \
  --min-cluster-size 2 --min-segment-bp 0 --chromosome-track --no-consensus-track \
  -o $S.clusters.svg
```

The whole runbook is wrapped as `scripts/run_cluster_pipeline.sh` (validates inputs, derives
sidecar names, falls back to an auto-palette without `colors.tsv`); pass `--sample`/`--prefix`/
`--db`. Verified to reproduce `plots_preview/U2OS.final.all_clusters.svg` byte-for-byte.
