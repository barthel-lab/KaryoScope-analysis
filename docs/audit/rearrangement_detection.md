# Rearrangement detection — model, design, and Engine A spec

> Status: **design agreed; Engine A fully implemented** — measurement
> (`core/colocalization.py`), differential test (`core/rearrangement.py`), and the
> `detect-rearrangements` CLI. The statistics are a **v1 for coauthor review** (the
> **Group A** bucket in `OPEN_QUESTIONS.md`): CMH across length buckets + BH-FDR +
> recurrence/effect/floor gates, with read independence assumed (surfaced as a runtime
> warning, not enforced — open item). Engine B (alignment-based clustering) is specified
> here but not yet built.

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

## 8. Engine B — clustering / consensus / viz (sketch; later)

Orientation-aware **alignment of feature-segment sequences**: local (Smith-Waterman-style)
on the symbol sequence (= the coalesced overlay), **length as a confidence weight** (+
`min_occurrence_bp` denoise), **substitution cost from hierarchy distance** (`HSat1A↔HSat1B`
cheap, `arm↔telomere` dear; `novel` neutral), affine gaps, and **best of forward/reverse**
so a read and its reverse have ~0 distance (kills mirror-image clusters by construction).
The pairwise distances feed a precomputed-distance clusterer (hierarchical / density-based),
and members are oriented to a seed for a consensus feature sequence and visualization.
Cluster-level differential testing is **deferred** (a cluster is a modeled object → testing
it risks circularity).

## 9. Open knobs / deferred

- generous `W` scale (report a couple for now); length-bucket boundaries; `min_occurrence_bp`.
- explicit copy-number model — **rejected for now** in favour of the rate comparison.
- distance-shift test (compare A–B *distance distributions* exp vs control) — optional refinement.
- self-pairs (`A` close to another `A`, for tandem amplification) — possible extension.
- inversions — treated as rare; revisit if data demands.
