# `build-feature-matrix` — metrics reference (for review)

What every column in the per-sequence feature matrix means, how it's computed (with code refs in
`KaryoScope_sequence_annotate.py`), and which constants are unexplained "voodoo" numbers we need to
justify + parameterize. **Please annotate where you (dis)agree.**

The matrix has **one row per `seq_id`**; it is the direct input to `cluster`. Per-feature columns are
emitted for every `(featureset, feature)` pair.

---

## A. Proposed column schema (resolves your parsability concern)

Today columns are `{featureset}_{metric}__{feature}` — the single `_` between featureset and metric is
**ambiguous** because both contain underscores. Proposed unambiguous scheme using `__` as the *only*
delimiter (no featureset/metric/feature token contains a double underscore):

| Kind | Pattern | Example |
|------|---------|---------|
| per-(featureset, feature) | `{featureset}__{metric}__{feature}` | `region__frac__bSat`, `region__dmax__bSat` |
| per-featureset total | `{featureset}__total_bp` | `region__total_bp` |
| interspersion (global) | `interspersion__{type}` | `interspersion__total`, `interspersion__tel_sat` |
| alignment (optional) | `align__{metric}` | `align__primary_mapq`, `align__read_length` |
| key column | `seq_id` | (currently `sequence`) |

**This is a shared contract**: `cluster` and `cluster-annotate` parse these column names back
(`cluster-annotate` builds `{pfx}__{col_kind}__{name}` lookups). So the delimiter change is made **once**
and applied across `build-feature-matrix → cluster → cluster-annotate` together.

---

## B. Coverage metrics (per featureset × feature)

| metric | definition | code | notes |
|--------|-----------|------|-------|
| `bp` | total bp of this feature on the sequence | `compute_per_read_feature_bp` L451-460 | `groupby(seq,feature).length.sum()` |
| `frac` | `bp / total_bp` — fraction of the sequence's annotated length that is this feature | `compute_read_feature_fractions` L199-210 | denominator = sum of *all* feature lengths on the seq |
| `total_bp` | total annotated bp on the sequence (the `frac` denominator) | L205/459 | **assumes** intervals within a featureset partition the sequence (no overlap), else inflated — ties to `overlay-annotations` full-tiling assumption |

## C. Local-density metrics (per featureset × feature) — sliding `window_size` (default 1000 bp)

`compute_per_read_window_densities_bulk` L314-388. For each feature on each sequence, build a per-bp 0/1
coverage array, then slide a `window_size`-bp window; **density of a window = (#covered bp)/window_size**.

| metric (token) | dict key | meaning |
|----------------|----------|---------|
| `dmax` | `max` | densest window — how *concentrated* the feature is |
| `dmin` | `min` | least dense window |
| `dmedian` | `median` | typical local density |
| `dfirst` | `first` | density of the first window (5′ end) |
| `dlast` | `last` | density of the last window (3′ end) |
| `dterminal` | `terminal` | `max(dfirst, dlast)` — "is the feature at *an* end?" (e.g. telomere) |
| `dterminal_min` | `terminal_min` | `min(dfirst, dlast)` |

Sequences shorter than `window_size` get their whole-sequence fraction broadcast to all density fields
(L344-358). **`window_size = 1000` is an unexplained default** → make it a documented CLI parameter.

## D. Contiguity metric

| metric | definition | code |
|--------|-----------|------|
| `max_block_bp` | longest *contiguous* run of the feature, **bridging gaps ≤ `BLOCK_GAP_TOL`** | `_max_block_length` L236-253 |

**What "gap" means (your question):** the coverage array is per-*feature*, so a "gap" is a stretch where
this feature is absent and a *different* feature is annotated. Even though every bp is annotated by
*something*, feature X's own track has gaps wherever X isn't the feature. `BLOCK_GAP_TOL = 100` bridges
runs of X separated by ≤100 bp of non-X, so a feature briefly interrupted still counts as one block.
**`100` is an unexplained heuristic** → document + parameterize.

## E. Interspersion metrics (global, per sequence) — transitions per kb

`compute_per_read_interspersion_bulk` L395-444 + `classify_bed_feature` L270-307. Sort intervals by start;
classify each into a category; count category *changes* between adjacent intervals, normalized per kb.

| metric | counts transitions between… |
|--------|------------------------------|
| `interspersion__total` | any two different categories |
| `interspersion__can_ncan` | canonical ↔ noncanonical telomere |
| `interspersion__tel_sat` | telomere ↔ satellite |
| `interspersion__arm_tel` | arm ↔ {telomere / ITS_TAR1} |

Categories come from the hierarchy-derived vocab (satellite / canonical / noncanonical / ITS_TAR1 / ct /
arm / other). **This is the only place the `:` composite separator matters** — see §G.

---

## F. The adaptive-thresholds sidecar — what it is and where it's used

**Definition** (`compute_adaptive_thresholds` L213-226): per feature,
`threshold = clamp(median(nonzero frac) / 3, 0.1%, 5%)`; features never seen nonzero get the 0.1% floor.
Written to `{output}.adaptive_thresholds.tsv` (plain text).

**Where it's used (audit correction):** it **is** consumed downstream — by `cluster-annotate`
(`load_adaptive_thresholds`, L92; `--adaptive-thresholds`, auto-derived from the matrix path, L1445-1456).
`cluster-annotate` uses each feature's threshold to compute **`readpct`** columns = "% of a cluster's
sequences where `frac__feature` exceeds the adaptive threshold" (its docstring L13, L155). So the threshold
is a **per-feature presence cutoff**: "how much of feature X must a sequence have to count as having it."

**Voodoo:** the `/3` factor and the `[0.1%, 5%]` clamp have **no documented rationale**. We need to
determine and document the rationale, and expose them as CLI parameters (D6).

---

## G. The `:` composite label contract (your "treated specially?" question)

Features **never contain `:`** in any current/intended database — `:` is reserved as the **composite
separator** produced by `overlay-annotations` (e.g. `DJ_TAR1`… actually composites use `_`; the two-layer
`a:b` form is the *overlay/composite* join output). The **only** place `build-feature-matrix` parses `:`
is `classify_bed_feature` (for interspersion): it splits `layer1:layer2` and applies a priority to pick one
category. Everywhere else a feature string is atomic (a composite `a:b` would simply be its own
`frac`/`bp`/density columns).

Consequence: interspersion behaves differently depending on whether the input is **raw single-featureset
BEDs** (`layer1` only) or **`overlay-annotations` composite output** (`layer1:layer2`). That's the shared
contract worth documenting + jointly testing. (For user-built databases: we document that `:` is reserved.)

---

## H. Alignment columns (optional; `--readnames-dir`) — what, source, downstream

**Source today (telogator-specific):** `load_readnames` L467-496 reads
`{dir}/{sample}/telogator/{sample}.readnames.txt` → `sequencing_approach`; `load_stats` L499-614 reads
`{dir}/{sample}/telogator/aligned/{sample}.{reference}.stats.tsv` → per-read alignment QC:
`read_length`, `primary_mapq`, `primary_de`, `primary_align_len`, `primary_align_fraction`,
`total_align_len`, `total_align_fraction`, `n_alignments`, `n_secondary`, `n_supplementary`, `max_mapq`,
`mean_de`, `is_mapped`.

**What they are:** per-sequence alignment quality (mapping quality, alignment fraction, divergence `de`,
multi-mapping counts). The "total aligned bases" sums `align_len` over **non-secondary** alignments
(primary + supplementary), which is a deliberate choice (inline rationale L543-547).

**Downstream use (traced):** consumed **only by `cluster-diagnostics`**, purely for **QC plots**
(`read_length`, `primary_align_fraction`, `primary_mapq`, `primary_align_len`, `n_alignments` by cluster).
They are **not used** by `cluster` (clustering), `cluster-annotate` (labels), `cluster-plot`, or
`select-representatives`. So they are optional diagnostic columns, not part of the core signal.

**Genericization (your ask):** decouple from telogator's directory layout — accept a generic
**per-alignment stats TSV** (documented required columns: `seq_id/readname, is_primary,
is_not_supplementary, align_len, read_len, mapq, de, align_fraction, is_mapped`) and an optional
seq→`sequencing_approach` TSV, by explicit path(s). Any aligner's per-alignment stats then work, not just
telogator's. **Open decision:** keep this as an optional join inside `build-feature-matrix`, or move it to
`cluster-diagnostics` (its only consumer)?

---

## Voodoo constants to justify + parameterize (D6)
- `window_size = 1000` (density window)
- `BLOCK_GAP_TOL = 100` (block gap bridging)
- adaptive threshold `/ 3` and clamp `[0.1%, 5%]`

## Open questions for maintainer
1. Do you agree with each metric definition above (esp. `frac` denominator, `dterminal`, interspersion
   transition typing)?
2. Confirm the `{featureset}__{metric}__{feature}` schema + namespaces (`interspersion__`, `align__`).
3. Alignment join: keep in `build-feature-matrix` (generic stats TSV) or move to `cluster-diagnostics`?
4. Rationale for the four voodoo constants (so we can document the defaults).
