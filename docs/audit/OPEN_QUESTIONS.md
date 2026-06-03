# Open questions register (Phase 2)

Everything still needing a **maintainer / coauthor** decision before Phase 3 (scaffold) and Phase 4
(migrate). Cross-cutting conventions (C1–C5) and decisions D1–D8, M*, F*, CA*, RP*, CC*, CL1–4 are already
**resolved** in `DECISIONS.md`; this file lists only what's still **open**, grouped by who/what it needs.
Each item cites its decision ID and the doc with the detail.

---

## A. Clustering statistics — for coauthor review (highest stakes)
Full rigor + proposed fixes in **`clustering_methods.md` §4, §7**. These can change scientific conclusions,
not just code.
- **CL-open #1 🔴 Circularity / double-dipping** — k is chosen partly to maximize group enrichment, then the
  same clusters are tested for enrichment. Make k-selection unsupervised by default + run enrichment as a
  separate step?
- **CL-open #2 🔴 Per-sample stats invalid** — min-p over samples + cluster-only FDR. Replace with full
  sample×cluster FDR or a per-cluster omnibus test?
- **CL-open #3 🟠 Compositional normalization** — replace the round-and-Fisher size-factor hack with a
  proportion/offset model? Keep or drop `telomeric` vs `total`?
- **CL-open #4 🟠 Three "enriched" definitions** (odds>1.5 / raw p<0.05 / FDR q) — unify on one.
- **CL-open #5 🟠 Default matrix** — z-scoring the sparse edge block inflates rare transitions; change the
  default representation? Keep the edge block at all, or abundance-only?
- **CL-open #6 🟡 Default k-selector** — `composite-knee` is unstable; switch default to silhouette/Calinski?
- **CL-open #7 🟠 `--exclude-features` default** removes `canonical_telomere*` (the dominant telomere signal)
  from clustering by default — intended?

## B. Scope / naming consolidations (the "generalize over-specific scope" calls)
- **telogator-reads-viz** (item 11) — it renders *any* per-sequence feature BED as vertical bars and is a
  **subset of `plot-reads`**. **Merge into `plot-reads` as a vertical preset, or keep a separate
  generally-named tool** that shares the renderer? (If separate: pick the general name.)
- **Representative selection** (CA5 / RP-open / B6) — `cluster-annotate` and `select-representatives`
  implement *different* algorithms for the same goal. Unify into `core/representatives.py` — **which strategy
  is authoritative**, and should `centroid_distance` actually drive selection (its name implies it does, but
  no code uses it)?
- **`cluster-plot` layout engines** (item 10) — unify the two near-duplicate vertical/horizontal engines, or
  keep separate? (Real effort either way.)
- **`animate` fixed-zoom mode** (item 12) — still used, or has adaptive superseded it? Dropping it removes
  most of the duplication.

## C. Biology to parameterize / relabel (needs scientific intent)
- **Type I ALT relabel** (CA-open #1) — sample-name-specific labeling baked into `cluster-annotate`. Extract
  to a config-driven post-step, or fully parameterize?
- **Auto-label thresholds** (CA-open #2) — ship a documented default preset (human/CHM13) + override, or
  require explicit config?
- **Translocation label color** (item 15) — red/blue hard-coded to `chr2_chr13`; source from a palette/config?

## D. Hidden-constant rationale to supply (we'll parameterize with these as defaults either way)
- **`build-feature-matrix`** (F4) — rationale for `window_size=1000`, `BLOCK_GAP_TOL=100`, adaptive threshold
  `/3` + clamp `[0.1%,5%]`.
- **`cluster`** (CL3) — `min_k=40`, composite weights `0.5/0.1/0.4`, purity `1.0/0.8`, `reduce_dims=500`,
  `early_stopping=150`, odds cutoff `1.5`, etc.

## E. Pipeline / contract questions
- **`find-translocation-reads`** (item 13) — drop unused `_bp`/`_pct` columns? support `chrX/chrY`
  translocations? merge overlaps before summing coverage? read length from an authoritative source vs
  `max(end)`?
- **`cluster-translocation-reads`** (item 14) — make the frozen cluster recipe configurable? is a `primary`
  group required? in-process vs subprocess after migration? consume `find`'s TSV instead of re-globbing?
- **`visualize-translocation-reads`** (item 15) — reconcile the two BED naming conventions / two length
  sources between modes; is the transposed `chromosome.smoothed.{trans}` name intentional?
- **`draw-legend`** (item 8) — push the 4 input-shaping helpers into `karyoplot.svg.legend`, or keep
  analysis-local? (Lean: karyoplot.)
- **`karyoplot` push-down scope** (decision #3) — confirm OK to add to `karyoplot`: a new `svg.dendrogram`,
  `video` module (D7), bubble/matrix drawers, unified legend, PIL raster path, and to route Fisher/FDR through
  `mpl.statistics`. (D4.5 already approved editing `karyoplot.core.colors` for `colors.tsv`.)

## F. Dependencies / tooling confirmations
- **Package deps:** `karyoscope-analysis` will depend on the `karyoscope` package (for `paths`/DB resolution,
  C5/D4.5) and `karyoplot`. Confirm.
- **v2 example data / fixtures** — maintainer to provide a `KS_human_CHM13_v2` example dataset for tests
  (D4.1 consequence; in progress).
- **v1→v2 migration utility** (D4.3) — confirm still deferred (build later only if old v1 result files must
  be reread).
- **DB-file inconsistency** (D4.6) — resolved (`repeat` added to `hierarchy.tsv`); no action.

---

### Status
Phase 1 audit: complete (16 per-script docs). Phase 2 scoping: complete (Review items 1–15 in
`DECISIONS.md`; structure mode dropped). Remaining before Phase 3: resolve the items above (A–F). No code in
`scripts/` has been modified.
