# KaryoScope-analysis — Script Audit (Phase 1 synthesis)

This directory holds a per-script audit of every tool in `scripts/`, produced before
restructuring the repo into an installable package modeled after
[`KaryoScope`](../../../KaryoScope) (core engine) and
[`KaryoScope-plotlib`](../../../KaryoScope-plotlib) (`karyoplot`, shared plotting library).

**Goal of the audit:** let the maintainer confirm the *scope*, *purpose*, *decisions*, and
*assumptions* of each script, surface bugs/dead code/duplication, and propose where each piece
of logic should live in the new layout — **before** any code is changed.

One file per script. This README is the cross-cutting synthesis. Each per-script file follows the
same 11-section template (purpose, CLI surface, inputs/outputs, control flow, design decisions,
assumptions, dependencies, proposed home, smells/risks, testability, open questions).

---

## 1. Inventory (16 scripts, ~24,762 lines)

| Script | Lines | Purpose (1-line) | Proposed subcommand |
|---|---:|---|---|
| `_feature_vocab.py` | 109 | Shared v1/v2 satellite/arm/telomere vocabulary (library module) | *(→ `core/feature_vocab.py`)* |
| `KaryoScope_draw_legend.py` | 195 | Standalone legend SVG from a `.colors.txt` file | `draw-legend` |
| `KaryoScope_find_translocation_reads.py` | 327 | Discover translocation reads; per-chrom bp/pct coverage TSV | `find-translocation-reads` |
| `KaryoScope_cluster_translocation_reads.py` | 519 | Orchestrate translocation clustering (shells to merge/cluster/plot) | `cluster-translocation-reads` |
| `KaryoScope_select_representatives.py` | 519 | Pick representative reads per cluster for `cluster-plot` | `select-representatives` |
| `KaryoScope_compare_clusterings.py` | 577 | Compare two clusterings of the same reads (ARI/NMI, concordance) | `compare-clusterings` |
| `KaryoScope_visualize_translocation_reads.py` | 734 | Per-read stacked-track SVG/PNG for translocation reads | `visualize-translocation-reads` |
| `KaryoScope_telogator_reads_viz.py` | 856 | Telomeric read feature-bar viz (Telogator output) | `telogator-reads-viz` |
| `KaryoScope_sequence_annotate.py` | 1065 | Per-read wide feature-matrix TSV from smoothed BEDs | `annotate-sequences` / `seq-features` |
| `create_panning_animation.py` | 1195 | PNG → panning MP4 (16:9 slides) | `animate` |
| `KaryoScope_merge_beds.py` | 1280 | Per-read BED overlay / priority-feature merge | `merge-beds` |
| `KaryoScope_cluster_annotate.py` | 1719 | Aggregate per-read annotations → per-cluster + auto-label | `cluster-annotate` |
| `KaryoScope_cluster_diagnostics.py` | 1731 | Post-hoc per-cluster metric & comparison figures | `cluster-diagnostics` |
| `KaryoScope_cluster_analysis.py` | 3103 | Feature matrix → Ward clustering → enrichment (the core step) | `cluster` |
| `KaryoScope_plot_reads.py` | 3620 | Stacked read feature-bar figures (+ animation hook) | `plot-reads` |
| `KaryoScope_cluster_plot.py` | 7213 | Cluster-representative figure (dendrograms / matrices / bubbles) | `cluster-plot` |

---

## 2. Pipeline / how the tools relate

```
                 (KaryoScope engine produces *.smoothed.features.bed.gz)
                                      │
        ┌─────────────────────────────┼─────────────────────────────────┐
        ▼                             ▼                                   ▼
  merge-beds                  annotate-sequences                  find-translocation-reads
  (overlay/priority BED)      (wide per-read matrix TSV)          (per-chrom coverage TSV)
        │                             │                                   │
        └──────────────┬──────────────┘                                   ▼
                       ▼                                        visualize-translocation-reads
                   cluster   ──►  NPZ + TSVs                    cluster-translocation-reads
            (feature matrix, Ward,                              (orchestrator: shells out to
             k-selection, enrichment)                            merge-beds/cluster/cluster-plot)
                       │
        ┌──────────────┼───────────────┬────────────────┐
        ▼              ▼               ▼                ▼
  cluster-annotate  cluster-plot   cluster-diagnostics  compare-clusterings
  (per-cluster      (the big        (post-hoc figures)  (two clusterings)
   labels)           figure)
        │
        ▼
  select-representatives ──► reads-file consumed by cluster-plot

  Independent viz: plot-reads, telogator-reads-viz ──► animate (panning MP4)
```

Key coupling facts (from the audits):
- `find-translocation-reads` and `cluster-translocation-reads` are **siblings**, not a pipe — each
  re-globs the results dir with its own regex (cluster's adds a `featureset` group). `find`'s TSV is
  consumed only by `visualize`'s `--input-tsv` mode.
- `cluster-translocation-reads` currently **shells out** (`subprocess`) to `merge_beds`,
  `cluster_analysis`, `cluster_plot` via a `--scripts-dir`. Post-packaging these become in-process calls.
- `plot-reads` imports `create_panning_animation` via a `sys.path.insert` hack — becomes a real import.
- `cluster-annotate` and `sequence-annotate` import `_feature_vocab` via `sys.path.insert`.

---

## 3. Cross-cutting duplication map (the heart of the karyoplot push-down)

These patterns are reimplemented across many scripts. Consolidating them is the bulk of the
"push shared logic into `karyoplot`" work (decision #3).

| Shared pattern | Reimplemented in | Proposed home |
|---|---|---|
| `TeeLogger` (stdout+file tee) | ~9 scripts (cluster_plot, cluster_annotate, merge_beds, sequence_annotate, all 3 translocation tools, …) | `karyoplot.core.io` or analysis `core/logging.py` |
| Scale-bar (pick round bp + draw) | plot_reads, telogator (×2), visualize_translocation; **cluster_plot already uses** `karyoplot.core.coords`/`text` | `karyoplot.core.coords.pick_round_scale_bp` + `text.format_genomic_distance` + an SVG drawer |
| Legend layout + render | plot_reads (full stack), draw_legend, visualize_translocation | `karyoplot.svg.legend` / `mpl.legend` |
| `hex_to_rgb` / `hex_to_rgba` | plot_reads, telogator | `karyoplot.core.colors` |
| Read feature-bar rasterization | cluster_plot, plot_reads (SVG path already via `karyoplot.svg.reads`; **PIL path still local**) | extend `karyoplot.svg.reads` with the PIL path |
| Stats: Fisher / FDR / Mann-Whitney / rank-biserial / significance stars | cluster_analysis, cluster_diagnostics, cluster_plot, compare_clusterings | `karyoplot.mpl.statistics` (partly exists; reconcile 3- vs 4-tier stars) |
| Plot style / dark-mode rcParams / dual-background save | cluster_analysis, cluster_diagnostics, cluster_annotate, compare_clusterings | `karyoplot.mpl.style` (exists) |
| BED loaders | cluster_analysis, select_representatives, telogator, merge_beds, sequence_annotate | `karyoplot.core.io.load_bed` (+ analysis `core/io` for specialized readers) |
| Results-layout discovery regex (`sample.data_type.replicate.db…featureset`) | all 3 translocation tools, plot_reads, cluster_plot | analysis `core/io/result_layout.py` |
| `svg_to_png` / rsvg-convert | visualize_translocation, telogator, create_panning_animation; cluster_plot delegates to `karyoplot.svg.export` | `karyoplot.svg.export` (exists) |
| ffmpeg / video encode | create_panning_animation, plot_reads (via sibling) | **no home yet** — new `karyoplot.video`? (open question) |
| Font loading (`pil_font` / `register_fonts`) | create_panning_animation has local `_load_font`; others delegate | `karyoplot.core.fonts` |
| Feature-vocab literals | `_feature_vocab` (2 users), plot_reads (own set), merge_beds (own **v1-only** set) | single source: analysis `core/feature_vocab.py` |

---

## 4. Flagged bugs / correctness concerns (⚠ verify before acting)

These were surfaced by reading, not by running. They need a quick runtime confirmation, but each is a
concrete decision point for "do we change the code?":

| # | Script | Finding | Likely effect | Confidence |
|---|---|---|---|---|
| B1 | compare_clusterings | imports `adjusted_rand_index` (real name is `adjusted_rand_score`); `ImportError` is caught silently | **ARI has likely never been computed** in any report (NMI is fine) | High |
| B2 | cluster_analysis | structure mode `groupby(['read', …])` but the column is `sequence` (~L355) | structure mode likely broken for multi-read chromosomes | Medium |
| B3 | merge_beds | priority-feature sets are **v1-only literals**, not sourced from `_feature_vocab` | v2 names (`arm`, `bSat`, …) silently misrouted | Medium-High |
| B4 | visualize_translocation | if `rsvg-convert` is missing, the SVG is deleted in a `finally` | run produces **no output**, silently | High |
| B5 | cluster_analysis | three inconsistent significance rules across the pipeline (k-loop odds>1.5, raw p<0.05, FDR q<thresh); per-sample FDR applied only over per-cluster min-p | inconsistent / under-corrected enrichment calls | Medium |
| B6 | select_representatives | `centroid_distance` is computed but never drives selection; required `--cluster-analysis` arg is unused | flag/name misleads; option does nothing | High |

Plus pervasive **silent-failure** patterns: broad `except Exception: pass`, malformed BED lines skipped
without warning, unseeded jitter (non-reproducible), `--format` honored only on some outputs.

---

## 5. Cross-cutting risk themes

1. **All 16 scripts use `argparse`**, not the gold-standard `click` CLI. Converting to a single click
   group with subcommands is the central structural change.
2. **Import-time execution / no `main()` guard** (notably `merge_beds`) blocks importing for tests —
   must be wrapped before it can be unit-tested.
3. **`sys.path.insert` hacks** to import sibling scripts / `_feature_vocab` — disappear once packaged.
4. **Hard-coded single-experiment biology**: `cluster_diagnostics` (`'E6E7'`/`'primary'`),
   `cluster_annotate` (Type I ALT relabel keyed on sample names), `merge_beds` (v1-only priority sets),
   hard-coded `…/1/KaryoScope/{database}/` BED paths in several tools. Limits generality.
5. **Dead code / unimplemented flags**: documented flags that do nothing (`cluster_plot`
   `--enrichment-normalization total`, `--orient-telomere-top`, …; `plot_reads --animate-crop-ratio`;
   `create_panning_animation --uniform-zoom`); dead functions in several files.
6. **Feature-vocab fragmentation**: three separate vocabularies, only one (`_feature_vocab`) v1/v2-aware.

---

## 6. Proposed target layout (to refine during migration)

Modeled on `KaryoScope` (src layout, click CLI, `commands/` + `core/` + `core/io/`):

```
src/karyoscope_analysis/
  __init__.py  _version.py  __main__.py
  cli.py                       # click group; registers all subcommands
  commands/                    # thin click wrappers, one per tool (arg parsing only)
    merge_beds.py  annotate_sequences.py
    cluster.py  cluster_diagnostics.py  cluster_annotate.py
    select_representatives.py  compare_clusterings.py
    cluster_plot.py  plot_reads.py  draw_legend.py  animate.py
    telogator_reads_viz.py
    find_translocation_reads.py  cluster_translocation_reads.py  visualize_translocation_reads.py
  core/
    feature_vocab.py           # ex-_feature_vocab.py (single source of truth)
    logging.py                 # TeeLogger (unless pushed to karyoplot)
    merge_beds.py              # interval/overlay/priority logic
    seq_features.py            # per-read matrix builders
    clustering/                # matrix.py, cluster.py, enrichment.py, kselect.py, structure.py
    representatives.py         # shared by cluster-annotate + select-representatives
    feature_importance.py
    translocations/            # discovery, coverage, layout shared by the 3 trans tools
    plot_reads/                # decomposition of the 3.6k-line script
    cluster_plot/              # decomposition of the 7.2k-line script
    io/
      bed.py  results.py  result_layout.py  colors.py  reads_bed.py  ...
tests/                         # one test file per core module + integration + golden
examples/  (data/)             # example inputs; decide vs tests/data
```

`karyoplot` push-downs (separate repo, decision #3): scale-bar drawer, legend, stats reconciliation,
style, color helpers, PIL rasterization, SVG `dendrogram` module (new), video/ffmpeg module (new?).

---

## 7. Consolidated decision points for the maintainer

Per-script questions live in each audit file (§11). The cross-cutting decisions that shape the whole
migration:

- **D1 — Output stability.** Your history validates refactors as *byte-identical* (`bench: 18
  byte-identical`). Converging scale-bar/legend/stats onto `karyoplot` will change output bytes. Is
  byte-identical output a hard requirement (→ keep some rendering local / extract conservatively), or
  is visual-equivalent acceptable for some tools?
- **D2 — Bug policy.** For B1–B6 above: fix as we migrate, or preserve current behavior and decide
  per-bug after you review?
- **D3 — CLI naming & splits.** `annotate-sequences` vs collision with engine's `annotate`
  (rename to `seq-features`?); split `cluster` structure-mode into its own `cluster-structure`?
  Final subcommand names.
- **D4 — Feature-vocab consolidation.** One v1/v2-aware source in `core/feature_vocab.py`, and migrate
  `plot_reads` + `merge_beds` off their private literals (fixes B3)?
- **D5 — Dead code / unimplemented flags.** Remove, or implement, or keep as-is?
- **D6 — Hard-coded biology.** Generalize/parametrize the single-experiment assumptions
  (E6E7/primary, Type I ALT), or keep for now?
- **D7 — Video home.** New `karyoplot.video` module vs analysis-local for the ffmpeg encoder?
- **D8 — Build backend.** `hatchling` (match `KaryoScope`) vs `setuptools` (match `karyoplot`)?

---

## Status

- Phase 1 (audit): **complete** — 16 per-script files + this synthesis.
- Phase 2 (scoping/review): **complete** — every script scoped; decisions recorded in **`DECISIONS.md`**
  (Review items 1–15 + conventions C1–C5 + cross-cutting D1–D8). Deep dives: **`KNOWN_ISSUES.md`** (v1→v2
  vocab errors + bugs for lab review), **`feature_matrix_metrics.md`**, **`clustering_methods.md`** (rigorous
  stats critique). **`OPEN_QUESTIONS.md`** is the consolidated list of what still needs maintainer/coauthor
  input before building.
- Phases 3–6 (scaffold → migrate → test → docs): not started; gated on the `OPEN_QUESTIONS.md` items.

> Note: some subcommand names/splits were refined during scoping (authoritative in `DECISIONS.md`): `merge-beds
> → overlay-annotations`; `annotate-sequences → build-feature-matrix`; `cluster` → **`build-matrix` +
> `cluster` + `test-enrichment`** (with `cluster`/`cluster-plot` **structure mode dropped**); `telogator-reads-viz`
> may merge into `plot-reads`. `_feature_vocab` becomes a hierarchy-derived `core/feature_vocab.py`.
