# Changelog

All notable changes to KaryoScope-analysis are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Reorganizing the historical `scripts/` collection into an installable package
(`karyoscope_analysis`) with a unified `karyoscope-analysis` CLI, modeled on the
core KaryoScope engine. See `docs/audit/` for the full audit and decision record.

### Added

- Package skeleton: `src/karyoscope_analysis/` (src layout), hatchling build,
  `karyoscope-analysis` console entry point, `version` subcommand.
- Tooling: `pyproject.toml` (ruff, pytest, coverage), `.pre-commit-config.yaml`,
  GitHub Actions CI (lint + test matrix), issue/PR templates, community files
  (`CONTRIBUTING`, `CODE_OF_CONDUCT`, `CITATION`, `LICENSE`).
- `docs/audit/`: per-script audits, scoping `DECISIONS.md`, `KNOWN_ISSUES.md`,
  `feature_matrix_metrics.md`, `clustering_methods.md`, and `OPEN_QUESTIONS.md`.
- `core/feature_vocab.py`: the hierarchy-derived, **v2-only** feature vocabulary
  (satellite/arm/ct/telomere groups loaded from the database `hierarchy.tsv` —
  satellites = the centromeric subtree minus `ct`; legacy v1 names rejected;
  `novel` the only out-of-taxonomy feature accepted). Single source of truth for
  feature groupings, replacing the old `scripts/_feature_vocab.py` constants.
  Ships with tests and a committed `tests/data/hierarchy.tsv` fixture.
- `core/intervals.py`: pure per-`seq_id` interval algebra (`coalesce`, `refine`,
  `merge_overlapping`, `total_covered`) — the sweep-line workhorses behind
  `overlay-annotations` (no pyranges, decision M5).
- `core/io/bed.py`: annotation-BED reader/writer enforcing the **C4 invariant**
  (rows grouped by `seq_id`; each sequence's intervals form a gapless,
  non-overlapping partition) — malformed input and gaps/overlaps are errors (C2);
  `.gz` (gzip/bgzip) read transparently. Provides both an eager
  `read_annotation_bed` and a streaming `iter_annotation_rows` (one row at a time,
  same C4 validation) for single-pass overlays.
- `core/annotation_resolution.py`: the `overlay-annotations` resolution engine —
  a `precedence` (default winner) + an ordered list of class-based `rules` (M2),
  with `emit` forms passthrough / `{literal}` / `composite`. Specs are validated
  structurally (jsonschema) and semantically against the database hierarchy (every
  featureset / feature / `@class` must exist); `when` keys are lint-checked to be in
  precedence order. Composite labels join in precedence order (e.g. `DJ_TAR1`).
- Built-in overlay presets (`presets/*.yaml`): the four legacy `merge_beds` priority
  modes ported 1:1 — `telomere-satellite`, `priority`, `chromosome-acrocentric`,
  `telomere-acrocentric` — with v1→v2 name translation (`telomere_like_multigroup1`→
  `telomere_like`, `arm_multigroup1`→`arm`, `array_multigroup1`→`array`,
  `acrocentric_multigroup1`→`acrocentric`, `noncentromeric`→`rDNA`). Loaded via
  `load_builtin_preset`; validated against the hierarchy.
- **`bin-annotations` subcommand** + `core/annotation_binning.py`: a **hierarchy-aware
  rolling-window mode filter**, run *before* `overlay-annotations` to denoise a single
  featureset BED. Each base's feature is replaced by the locally dominant feature in a
  centered window (default 101 bp, clipped at sequence ends; output length and the C4
  partition are preserved). The vote is **not flat** — per-feature window bp propagate up
  the database tree and the call is found by descending from the root into the dominant
  child while its subtree holds a majority, stopping at the deepest node that still does;
  so related siblings (e.g. `aSat`/`bSat` under `centromeric`) reinforce each other rather
  than splitting their vote and losing to an unrelated minority, and a split among subtypes
  honestly reports their ancestor. Knobs: `--majority-fraction` (τ, default 0.5; **0 = always
  descend to a specific leaf**, no internal/ambiguous labels) and
  `--threshold-scope` (`node` = conditional majority relative to the current node, more
  specific [default]; `window` = relative to the whole window, conservative); a window with
  no top-level majority falls back to flat plurality (keeps clean boundaries sharp). `novel`
  votes as a top-level leaf; every other label must be in the featureset (C2). Streams one
  sequence at a time (`O(one sequence)` memory), atomic output. Strongly fragmentation-
  reducing on real data (U2OS region: 1.38M → 219k intervals at window 101, 6.3×; → 30k at
  window 1001, 46×). **O(intervals) and window-independent**: between the O(intervals)
  breakpoints where the window's entering/leaving base crosses an interval edge the per-feature
  counts are linear, so the descent is evaluated O(1) times per segment (a conservative
  "constant for the next N steps" bound; recompute-and-merge keeps it exact) — cost no longer
  grows with the window (U2OS region ~23–28 s at *either* window 101 or 1001, vs ~160 s for the
  original per-base version). A per-base reference (`bin_intervals_naive`) pins the fast path in
  a 600-case property test (τ ∈ {0, …, 1}, both scopes).
- **`overlay-annotations` subcommand** (replaces `KaryoScope_merge_beds.py`): a
  **single-pass, streaming k-way overlay**. It reads every per-featureset BED
  concurrently and sweeps a line across the union of track boundaries per `seq_id`,
  resolving each segment via a preset / custom spec / the default basic overlay and
  coalescing — holding only the *current* interval of each track, so it runs in
  `O(featuresets)` memory regardless of input size (peak ~2 MB vs ~420 MB for the
  in-memory approach on a full sample, byte-identical output). Tracks must present
  sequences in the same order (lockstep, validated — no fragile name comparison);
  order/coverage/span disagreements are errors. Validates input features against the
  hierarchy as they are first seen (C2), and writes output atomically (temp file +
  replace, so a mid-stream error leaves no partial file). First fully-migrated tool,
  wired into the `karyoscope-analysis` CLI.
- `core/seq_features.py`: pure per-`seq_id` feature metrics for `build-feature-matrix`
  — coverage (`bp`/`frac`/`total_bp`), 1-bp-step sliding-window density
  (`dmax`/…/`dterminal`), `max_block_bp` (gap-bridged), hierarchy-derived
  interspersion, and the adaptive-threshold computation. Every magic constant
  (window size, block-gap tolerance, threshold factor/bounds) is a parameter (F4).
  The sliding-window density is computed **analytically from the intervals**
  (`O(intervals + window)`, no dense per-feature coverage array, cumsum, or
  full-array median) — byte-identical to the old dense computation, verified by a
  property test over thousands of random partitions.
- **`build-feature-matrix` subcommand** (replaces the matrix-building part of
  `KaryoScope_sequence_annotate.py`): a **single-pass, streaming** build — every
  featureset BED is read concurrently in lockstep by `seq_id` (via
  `core/io/bed.iter_aligned_groups`), so only one sequence's intervals are held at a
  time. Emits the wide per-`seq_id` matrix (`{featureset}__{metric}__{feature}`
  schema, F2) + an adaptive-threshold sidecar (F5). Constants are CLI options (F4);
  alignment-QC columns intentionally omitted (they move to `cluster-diagnostics`, F6).
  On a full sample: peak memory ~1244 MB → ~218 MB (5.7×) and read+compute ~58s → ~21s
  (2.8×), output byte-identical. Second fully-migrated tool; completes the
  data-foundation tier.
- End-to-end tests for `overlay-annotations` and `build-feature-matrix` on real
  `KS_human_CHM13_v2` HeLa data: fast default tests run on tiny committed fixtures
  (`tests/data/v2_subset/`, carved from the full BEDs by a documented, deterministic
  `make_subset.py`), and `@pytest.mark.integration` tests run the same pipelines on
  the full `data/raw_bed/` BEDs (skipped when that large data is absent). They assert
  the real contracts: C4-valid output, `seq_id` conservation, the F2 column schema,
  coverage self-consistency (`bp`/`frac` sum to `total_bp`/1), and telomere
  interspersion over an overlay composite.
- `docs/audit/rearrangement_detection.md`: the agreed mental model and design for
  detecting recurrent rearrangements as **abnormal feature colocalizations** —
  reads as orientation-agnostic ordered feature-segment sequences; the two axes of
  abnormality (proximity vs abundance); the Engine A (alignment-free colocalization
  detection) / Engine B (alignment-based clustering) split; the experiment-vs-control
  + annotated-CHM13-reads baselines; and the anti-circularity discipline (Group A).
- `core/colocalization.py`: Engine A **measurement layer** — per read, the minimum bp
  gap between every co-present feature pair (a single sweep, orientation-invariant,
  `min_occurrence_bp` denoise), plus a streaming reader over an overlay BED. Pure and
  non-statistical.
- `core/rearrangement.py` + **`detect-rearrangements` subcommand**: Engine A differential
  test. Aggregates each sample's per-read gaps into per-(length-bucket, pair) support
  counts, then tests experiment-vs-control colocalization rates per `(pair, window)` with a
  **Cochran-Mantel-Haenszel** stratified test (length buckets as strata), **BH-FDR**, and
  recurrence / effect-size / **reference artifact-floor** gates (the floor = annotated
  CHM13 reads). Reports all tested pairs with a `passes` flag and a `reference_abnormal`
  annotation. A **v1 for coauthor review** (Group A statistics); read independence is
  assumed and surfaced as a runtime warning, not enforced (de-duplication is an open item).
- `docs/audit/rearrangement_detection.md` §8: the full **Engine B** design — clustering as
  overlap-layout-consensus (OLC) assembly over feature-segment sequences (overlap graph on
  *proper* overlaps only, connected-components clustering, seed-anchored consensus,
  orientation propagation; recurrent inter-cluster bridges = rearrangements).
- `core/feature_align.py`: Engine B's **feature-sequence local aligner** — Smith-Waterman
  over `(feature, length)` segments with hierarchy-tiered substitution (`hierarchy_substitution`:
  exact > sibling > unrelated, `novel` neutral), `min(len)`-weighted (length-change-lenient)
  matches, linear per-bp gaps, and best-of-forward/reverse. Classifies overlaps as
  dovetail / containment / internal, with a feature-Jaccard prefilter. Pure; the overlap
  graph + consensus + CLI build on it.
- `core/feature_assembly.py`: Engine B's **overlap graph + clustering** — keeps only
  *proper* overlaps (dovetail/containment) clearing a minimum overlap length and normalized
  identity (internal-only repeat matches rejected; the anti-chaining safeguard), then groups
  reads into connected-component clusters via a **parity union-find** that assigns each read
  an orientation relative to its cluster seed (longest read) and flags orientation
  conflicts. Singletons kept. `assemble()` returns `(clusters, edges)`.
- `core/feature_assembly.py` `cluster_consensus`: Engine B's **seed-anchored consensus** —
  each cluster member is oriented to the seed and locally aligned to it, and every aligned
  column votes its feature at the corresponding seed position (the seed votes for itself).
  Yields the per-position majority feature with support/coverage over the seed backbone.
- **`cluster` subcommand**: Engine B's CLI — clusters an overlay-annotations BED subset of
  reads into structural haplotypes (OLC), writing a clusters table, a per-cluster consensus
  BED (with per-position support/coverage), and a member layout TSV (orientation + seed
  flag). Substitution scorer + all overlap thresholds are options; `--min-length` keeps the
  long-read-only default. Rendering of the layouts is left to the plotting tier. This
  completes Engine B (aligner -> overlap graph -> clusters -> consensus -> CLI).
- Engine B **feature weighting** (anti-chaining) — `core/feature_assembly` per-feature
  `weight` scales both the match reward and the overlap-length criterion so overlaps must
  rest on distinctive structure; `idf_weights` (frequency) and
  `FeatureHierarchy.interspersed_repeat_features` (the `Interspersed_Repeat` subtree). The
  `cluster --weight-method` defaults to `repeat-mask` (zero interspersed repeats + `nonrepeat`).
  Validated on real v2 data: repeat-mask cuts a 223/250-read mega-cluster to 127 (singletons
  11 -> 112); the residual is a structurally-coherent telomeric group whose sub-division is
  left to community detection (open). See `docs/audit/rearrangement_detection.md` §10.
- Engine B **chromosome identity** (the bigger anti-chaining lever) — a
  `chromosome-telomere-satellite` overlay preset emits two-layer `chromosome:structural`
  labels (e.g. `chr13:canonical_telomere`), enabled by a new **explicit-list composite emit**
  (`{composite: [chromosome, subtelomeric]}`) in the resolution spec. The aligner gains
  `chromosome_aware_substitution` (and the `cluster --cross-chromosome-penalty` option): the
  structural layer is hierarchy-graded; the chromosome layer adds a soft per-bp penalty unless
  two positions are the *same specific chromosome* (different/ambiguous chromosomes penalized),
  so reads cluster by chromosome and translocation reads still bridge. Validated: under
  `telomere-satellite` the subset chained almost entirely (249/250), and the chromosome
  composite splits it into
  clean per-chromosome clusters (largest 18) with multi-chromosome clusters surfacing as
  candidate translocations. See `docs/audit/rearrangement_detection.md` §10.
- **`genome-weights` subcommand** + `core/genome_weights.py`: per-feature **information-content
  weights from the annotated CHM13 reference** (one C4 BED per featureset). Tallies each
  feature's genome bp, takes its fraction `p` of its featureset's partition, and writes
  `-ln(p)` scaled to `(0, 1]` by the rarest feature across all featuresets (ubiquitous → ~0,
  rare → 1). New `cluster --weight-method genome-freq --genome-weights <tsv>` applies them.
  This is the principled answer to *structural* chaining (read-frequency `idf` couldn't
  separate the drivers): genome-wide, `q_arm`/`p_arm` cover 57%/26% of the genome → weight
  0.027/0.064 (crushed), while `canonical_telomere` covers 0.003% → weight 0.51 (kept
  informative). **Also fixes a latent bug:** `_weight`/`idf_weights` now key on the
  **structural layer** of `chromosome:structural` labels, so weighting actually applies on
  composite overlays (previously `repeat-mask`/`idf` silently no-op'd there). Validated on the
  U2OS w1001/τ0 binned subset: clusters become chromosome-coherent and recurrent
  candidate translocations surface (chr4+chr22 across 3 reads; chr12+chr9 across 2), where
  before reads were glued by shared `p_arm` with a meaningless consensus.
- Engine B **distinctive-overlap edge criterion** (`cluster --min-distinctive-bp`,
  `--distinctive-weight`) — the final anti-arm-chaining lever. An edge now also requires a
  minimum bp of matched *distinctive* features (weight ≥ `--distinctive-weight`, default 0.15),
  so an overlap explained only by filler (a shared chromosome arm) is rejected even though its
  weighted overlap is large (huge arm block × tiny weight). Genome-frequency weighting alone
  left ~150 kb arm blocks clearing the weighted-overlap floor; this closes that gap. Validated
  on the U2OS w1001/τ0 subset at `--min-distinctive-bp 1000`: the chr4 arm-star dissolves to
  singletons and the all-`p_arm` chr5 read drops out of the chr5 cluster, leaving three
  structurally-justified clusters (chr12+chr9; chr20; chr5+`canonical_telomere`). Off by default
  (0); recommended with a weighting method.
- Engine B **scales to a whole sample** — three levers so `cluster` no longer needs the
  O(N²) all-vs-all alignment that made 4005 reads intractable (43 min, unfinished):
  (1) a **blocking index** (`--block-min-bp`) that buckets reads by the *specific chromosome*
  leaf of their composite labels and only aligns within-bucket pairs (translocation reads join
  every chromosome they span, so the signal is kept; ambiguous `autosome`/`categorized` content
  is excluded); (2) **scorer memoization** (the substitution scorer is called per DP cell but has
  only ~hundreds of distinct feature pairs); (3) **process parallelism** (`--workers/-j`) over the
  candidate pairs via `fork`, so workers inherit the reads + scorer copy-on-write (no pickling)
  and the result is identical to serial. Full U2OS (4005 reads, w1001/τ0 overlay) now clusters in
  **~100 s on 8 cores** (≈6.4× speedup), vs >8 min serial / 43 min unoptimized. See the
  "running on a whole sample" runbook in `docs/audit/rearrangement_detection.md` §13.
- Engine B **community detection** (`cluster --communities`, default on) — in-house weighted
  **label propagation** subdivides each connected component so a sparse bridge (a noisy
  multi-chromosome read, or a lone translocation) no longer transitively merges distinct groups.
  Orientation parities still come from the component-wide union-find, so they stay consistent
  within each sub-community. On the full U2OS run this dissolved a **1058-read mega-cluster** (it
  spanned every chromosome — diagnosed as bridge-chaining through 52 noisy ≥3-chromosome reads,
  whose chromosome layer is scattered slivers, not the large blocks a real complex rearrangement
  would show) into **clean per-chromosome haplotype groups** (some chromosomes resolve into
  *multiple* haplotype communities, e.g. chr18 → 148 + 83) and **24 recurrent translocation
  candidates** (chr4+chr22 ×43, chr18+chr19 ×34, chr13+chr11 ×33, chr1+chr21 ×23, …). 4005 reads →
  2544 clusters (75 multi-read), no mega-cluster.
- Engine B **consensus-coordinate layout** (`feature_assembly.consensus_layout`, replacing the
  seed-anchored consensus + single-offset layout). Each read is placed in the cluster's consensus
  frame via a piecewise-linear map through its read→seed **alignment columns** (slope-1 outside
  the anchors), so **matched features stack vertically** instead of drifting (a chr20 cluster's
  shared `HSat3` now lands at identical coordinates in both reads, vs ~20 kb off before), and a
  read's overhang **extends the frame** so the consensus spans the **union** of all reads (not
  just the seed). `cluster`'s `layout.tsv` is now **per-segment** (`cluster_id, read_id, is_seed,
  reversed, start, end, feature` in consensus coords) and the consensus BED spans the union;
  `cluster-plot` draws straight from those coordinates (no `--overlay` needed) with length
  filtering off by default so a gap unambiguously means "didn't align".
- Engine B **filler-aware distinctiveness** — a `FeatureHierarchy.filler_features` set (telomere +
  arm + ct + non-* / novel) is now what "distinctive" means, replacing the genome-frequency weight
  threshold for clustering. **Telomere is genome-rare (high weight) yet read-set-ubiquitous** (every
  chromosome end), so the weight couldn't flag it; the explicit set can. Two effects: `cluster
  --min-interesting-bp` (default 2000) **drops "boring" reads** (telomere/arm-only, no satellite /
  ITS / TAR1 / rDNA content) before clustering, and the distinctive-overlap edge criterion counts
  only non-filler shared content, so reads can no longer **chain through shared telomere**. On full
  U2OS this dissolved the residual 323-read, 23-chromosome `cluster_0` (it was 53% boring reads) —
  every cluster is now chromosome-coherent: 39 single-chromosome haplotype groups + 14 clean
  translocation candidates (chr4+chr22 ×34, chr13+chr11 ×33, chr14+chr22 ×28, chr21+chr1 ×23, …).
  `cluster-plot` titles now list **all** major chromosomes (`chr4+chr22`), so translocations are
  labeled as such.
- **`cluster-plot` subcommand** + `core/cluster_plot.py` + `core/io/colors.py`: the package's
  single **read-renderer** (collapsing the legacy `plot-reads`/`cluster-plot`/`telogator-reads-viz`).
  Renders **one cluster (`--cluster-id`) or all clusters stacked in one SVG** (omit `--cluster-id`;
  `--min-cluster-size` filters, each panel titled with its dominant chromosome) — each read a row of
  feature-colored bars, oriented and offset into the seed frame (from `cluster`'s `layout.tsv`, which
  now carries a seed-relative `offset`/`length` via `feature_assembly.cluster_layout`), with the
  consensus track and a shared legend.
  Colors from the database `colors.tsv` (committed to `tests/data/`; structural layer of
  `chr:feature` labels), with a deterministic auto-palette fallback; `--min-segment-bp` denoises.
  Self-contained raw SVG (no plotting deps) — the `karyoplot.svg` push-down is deferred, as are the
  animation/video (D7) and the old `cluster_analysis`-specific visuals (dendrogram/enrichment bubbles).

### Notes

- Migration of the analysis scripts into the package (with bug fixes, v2-only
  feature vocabulary, and `karyoplot` push-down) is in progress; see
  `docs/audit/DECISIONS.md`.
