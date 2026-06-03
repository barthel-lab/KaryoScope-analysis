# Reorganization decisions

Living record of decisions made during the audit/restructuring of KaryoScope-analysis.
Cross-cutting decisions first; per-script scope/naming decisions appended as we review.

## Cross-cutting decisions (resolved 2026-06-03)

| ID | Decision |
|---|---|
| **D1 — Output stability** | **Correctness over byte-identical.** We do *not* require byte-identical output. Expect substantial output changes as bugs are fixed. Validate via tests + visual review, not byte-diff. The old `bench` byte-diff workflow is retired as an acceptance gate (may still be used opportunistically as a "did anything change?" signal). |
| **D2 — Bug policy** | **Fix bugs as we migrate.** B1–B6 (and any newly found) are real and get fixed during the relevant script's migration, with a test that locks the corrected behavior. |
| **D3 — CLI naming** | **Names are derived per-script from the true I/O contract**, not from the script's historical origin. Many names are over-specific (e.g. `telogator_reads_viz` really renders reads from *any* BED). Final subcommand names are decided one-by-one during the script review and recorded below. |
| **D4 — Feature vocabulary** | **Single v2-only source of truth** at `core/feature_vocab.py`. Drop all v1 names, the v1↔v2 alias tables, and the v2→v1 fallback logic. All consumers (cluster_annotate, sequence_annotate, merge_beds, plot_reads) use this one module. **Accepted consequence:** old v1 annotation TSVs are no longer understood; v1 example data/colors must be replaced with v2 fixtures. |
| **D5 — Dead code** | **Remove all of it** (dead functions, unimplemented/no-op flags). |
| **D6 — Hard-coded biology** | **Remove/generalize all of it.** Single-experiment assumptions (`E6E7`/`primary` enrichment, Type I ALT relabel, v1-only priority sets, hard-coded `…/1/KaryoScope/{db}/` paths) become parameters or config. |
| **D7 — Video home** | **New `karyoplot.video` module** for the ffmpeg encoder + panning-frame generation (sibling to `karyoplot.svg.export`). `animate` becomes a thin CLI over it; `plot-reads` calls it in-process (removes the `sys.path.insert` hack). |
| **D8 — Build backend** | **`hatchling`**, matching KaryoScope (the model). Path-based dynamic version via `_version.py`; clean sdist/wheel targets; no `setup.py`. (plotlib uses setuptools; that minor inconsistency is acceptable.) |

## Confirmed bugs to fix (maintainer reviewed; all true)

- **B1** `compare-clusterings`: `adjusted_rand_index` → `adjusted_rand_score`; stop swallowing the ImportError. ARI must actually compute.
- **B2** `cluster` structure mode: `groupby` on `sequence` not `read`.
- **B3** `merge-beds`: priority feature sets must come from the v2 vocab, not v1-only literals.
- **B4** `visualize-translocation`: do not delete the SVG when PNG conversion is unavailable; fail loudly or keep the SVG.
- **B5** `cluster`: unify the significance rules and FDR scope across the pipeline.
- **B6** `select-representatives`: make `centroid_distance` actually drive selection (or remove it); wire or drop the unused required `--cluster-analysis` arg.

## Review protocol (Phase 2)

Script-by-script "scoping card", in dependency order. For each script we confirm:
1. **What it actually does** (generalized beyond its origin)
2. **Input contract** — file types + essential format (vs incidental source)
3. **Output contract** — files + formats
4. **Correct, general name** (the D3 decision for this script)
5. Bugs to fix, dead code to cut, hard-coded biology to parameterize, vocab impact
6. Proposed home / decomposition (commands/ + core/ + core/io/, + karyoplot push-downs)

Decisions recorded below as we go.

## Cross-cutting conventions (apply repo-wide)

- **C1 — Sequence terminology.** BED column 1 is a **`seq_id`** — a generic sequence identifier (a read,
  contig, chromosome, …), not specifically a sequencing read. Use "sequence"/`seq_id` in code, docs, and
  CLI help, **not** "read", except in tools that are genuinely read-specific. Rename existing `read`
  variables/columns to `seq_id` during migration.
- **C2 — `novel` is the only out-of-taxonomy feature allowed.** `novel` = k-mers absent from the database;
  it is the single feature not present in `colors.tsv`/`hierarchy.tsv`/`features.tsv`. The vocab loader and
  all feature-consuming logic must accept `novel` (pass-through). **Any other feature not in the taxonomy is
  an error → the program terminates** (no silent pass-through).
- **C3 — bgzip output by default.** All KaryoScope-derived pipelines write **bgzip**-compressed outputs by
  default (consistent + indexable), matching the engine's BGZF inputs. (Tabix indexing only where genome
  coordinates make it meaningful; per-`seq_id` BEDs may not benefit.)
- **C4 — Annotation-BED input invariant (validated, not assumed).** Input annotation BEDs are **sorted by
  `(seq_id, start)`** and the intervals **partition** each `(seq_id, featureset)` — i.e. **non-overlapping
  and gapless** (every bp annotated exactly once, possibly `novel`). This is **validated once on read**
  (`core/io/bed.py`); any violation (overlap, gap, or out-of-order) is an **error** (C2). Downstream per-seq
  code **relies on this and never re-sorts** (drops the per-seq `sort_values` in interspersion, etc.). This
  subsumes `overlay-annotations`' coverage check (M4) and `build-feature-matrix`'s `total_bp` no-overlap
  assumption.
- **C5 — One result-layout module; no hard-coded paths.** Multiple scripts independently hard-code the
  telogator/KaryoScope on-disk layout (e.g. `{prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.
  telogator.1.{database}.{featureset}.{smoothness}.features.bed[.gz]`) and a sample/featureset **discovery
  regex**. Centralize both in **one** `core/io/result_layout.py` (path builders + discovery), make the layout
  **explicit/configurable** (resolve the DB dir via `karyoscope.paths`, D4.5), and **generalize off
  "telogator"** so it works for any sequence source (C1). No per-script path literals.

## Per-script decisions

### Review item 1 — feature vocabulary (foundation)

Definitive v2 sources (per maintainer): `KaryoScope-databases/KS_human_CHM13_v2/{features,hierarchy,colors}.tsv`.
Definitive v1→v2 rename map: `~/Downloads/karyoscope_feature_name_mapping.tsv`. The `hierarchy.tsv`
encodes the full v2 taxonomy as `(feature_set, child, parent)` rows; `colors.tsv` is `(feature_set, feature, color)`.

**Discrepancies found in `scripts/_feature_vocab.py` `SATELLITE_V1_TO_V2` (vs definitive map):**
1. **WRONG:** `hor_multigroup1 → hor`. Definitive: `hor_multigroup1 → alpha_hor`; v2 `hor` actually
   comes from v1 `inactive_specific`. → a real semantic error (selects the wrong column).
2. **MISSING:** no `inactive_specific → hor`, no `divergent_specific → dhor` (both listed in `SATELLITE_V1`
   but absent from the map).
3. **INCONSISTENT v1 key forms:** drops `_specific` (`bsat`, `gsat`, `censat`, `hsat1A/B/2/3`, `active`,
   `monomeric`) but keeps `_multigroup1` (`asat_multigroup1`, `hor_multigroup1`, …). Definitive v1 names
   carry their suffix → real v1 columns may never match.
4. **TYPO:** `'censet'` in `SATELLITE_V1` (should be `censat`).
5. **INCOMPLETE scope:** only satellites/arms are mapped. The definitive map covers all six feature sets
   (dozens of renames, e.g. `noncentromeric_specific → rDNA`, `repeat_multigroup1 → repeat`,
   `autosome_multigroup1 → autosome`) — none of which the module handles.

**plotlib (`karyoplot`) vs definitive v2 — FLAGGED, pending maintainer confirm before changing plotlib:**
- `karyoplot.core.colors` only reads the legacy per-featureset `{db}.{featureset}.colors.txt` format and
  manipulates `_specific`/`_multigroup1` suffixes (v1 artifacts). It has **no loader** for the definitive
  single `colors.tsv` (`feature_set, feature, color`).
- `karyoplot.svg.legend` strips `_specific`/`_multigroup1` (v1 baggage; inert in pure v2).
- `karyoplot.core.chromosomes` hard-codes acrocentric q_arm coordinates — genome geometry, not vocab; noted, likely fine.
→ In a v2-only world plotlib needs a `colors.tsv` reader, and the suffix logic becomes dead.

**Recommendation (pending decisions D4.1–D4.4 below):**
- Drop v1 from all analysis logic; v2-only.
- Replace hard-coded vocab sets with **derivation from `hierarchy.tsv`** (+ palette from `colors.tsv`) →
  one source of truth, satisfies D6 (no hard-coded biology) and D4 (single place).
- If old v1 result files must stay readable: a one-shot migration utility built from the definitive map
  (kept as a data asset, not threaded through analysis logic). Else skip.

**Decisions:**
- **D4.1 — RESOLVED:** Drop v1 from all analysis logic; v2-only.
- **D4.2 — RESOLVED:** Derive vocab from `hierarchy.tsv` (+ palette from `colors.tsv`) rather than
  hard-coded sets. Single source of truth; removes hard-coded biology per D6. `core/feature_vocab.py`
  becomes a loader over the hierarchy exposing e.g. `descendants(feature_set, root)`, with the satellite
  set = centromeric subtree minus `ct` (D4.4).
- **D4.3 — DEFERRED:** A v1→v2 one-shot migration utility for old result files may be built later, separate
  from the analysis path. Not now.
- **D4.4 — RESOLVED:** "Satellite" = the `centromeric` subtree **excluding `ct`** (preserve current logic).
  Rationale (maintainer): `ct` = "centromeric transition" = gaps *between* annotated centromeric satellites
  — they occur within centromeres but are not satellite sequences. Concretely: satellite set =
  `{centromeric} ∪ descendants("region","centromeric") \ {ct}`.
- **D4.5 — RESOLVED:** Update `karyoplot` to read the definitive `colors.tsv`. **Path resolution:** reuse
  KaryoScope's convention — `karyoscope.paths.default_db_root()` precedence (`--db` flag → `KARYOSCOPE_DB`
  env var → `~/.karyoscope/db/`) + database name → `<db_root>/<database>/colors.tsv` (and `hierarchy.tsv`).
  Keep `karyoplot` a *pure loader* (takes a file path); the analysis layer resolves the db path and passes
  it in. Implies `karyoscope-analysis` depends on the `karyoscope` package for path resolution.

**New findings recorded (see KNOWN_ISSUES.md):**
- Mapping file is **complete** on the v2 side (every v2 feature in hierarchy/colors/features is mapped).
- One db inconsistency: feature **`repeat`** was in `colors.tsv` + the mapping but **absent from
  `hierarchy.tsv`**. → **D4.6 RESOLVED:** maintainer updated `hierarchy.tsv`; `repeat` is now the repeat-set
  root (`repeat → categorized`, with `Interspersed_Repeat`/`Satellite`/`Noninterspersed` reparented under
  it). Cross-check now fully consistent (hierarchy = colors = mapping = 91).
- `plot_reads` and `merge_beds` carry their own **v1-only** vocab literals (V6/V7 in KNOWN_ISSUES) →
  must be migrated to the single v2 source.

### `karyoplot.core.chromosomes` (generality note)
- Provides chromosome ordering (`chrom_sort_key`), the human `ACROCENTRIC` set, `TELOMERIC_MOTIFS`
  (TTAGGG etc.), and a pluggable `Reference` (lengths + acrocentric `q_arm_starts`) registered as
  `CHM13_v2`. The acrocentric `q_arm_starts` are for drawing acrocentric ideograms (short arm = rDNA).
- **Currently NOT imported by any analysis script** — it's plotlib-internal/aspirational.
- Designed to be reference-pluggable (`register_reference()`), BUT `ACROCENTRIC` and `TELOMERIC_MOTIFS`
  are module-level **human** constants. For non-human genomes these would need generalizing. Low priority
  (unused by analysis today); tracked as a plotlib generalization item under D6.

### Review item 2 — `merge_beds` → `overlay-annotations`

**True scope:** a per-sequence **annotation-overlay** tool. Input = N≥2 annotation BEDs (4-col
`seq_id, start, end, feature`) over the *same* sequences; output = one resolved annotation BED. "Annotation"
= a 4-col BED. The "merge" name is dropped (it's an overlay/intersection, not a union).

**Decisions:**
- **M1 — Name: RESOLVED** — subcommand `overlay-annotations` (module `overlay_annotations`).
- **M2 — Resolution model: RESOLVED (YAML + jsonschema).** Every mode is a function
  `(features-per-featureset) → resolved_feature`, expressed **explicitly in a spec file**, not in if/then.
  To avoid the combinatorial product blow-up, the spec = **(a) `precedence`** (the *single* list of
  featuresets, in default-winner order — there is no separate `inputs` key) + **(b) an ordered list of
  `rules`** whose `when` patterns match feature *names*, *lists*, or *hierarchy classes* (`@arm` =
  descendants of `arm`), first match wins, unmatched positions fall to `precedence`.
  - `emit` forms: `<featureset>` (pass through that featureset's matched feature) · `{literal: X}` (force X)
    · `composite` (**join winners in `precedence` order** — deterministic, no accidental ordering).
  - **`when` keys are always written in `precedence` order** (enforced by a spec lint), so a rule reads
    left-to-right in priority order and the comment matches expectation without mental reordering. E.g. for
    precedence `[acrocentric, subtelomeric]`: `when: {acrocentric: [DJ, PHR, rDNA], subtelomeric: [TAR1, ITS]}`.
  - The 4 legacy modes ship as named **preset spec files** (`composite`/overlay, telomere-satellite,
    priority, chromosome-acrocentric, telomere-acrocentric); users may supply their own.
  - Presets are **ported 1:1 from the current code's resolvers** (`apply_conditional_region_repeat_rules`,
    `_resolve_telomere_acro_feature`, etc.), including v1→v2 name translation (e.g. `noncentromeric → rDNA`),
    and **reviewed with the maintainer** before finalizing. Confirmed faithful mapping example
    (telomere-acrocentric): precedence `[acrocentric, subtelomeric]` (acro default; composite → `DJ_TAR1`),
    + telomere-override / TAR1-ITS-composite / TAR1-ITS-keep rules.
- **M3 — Spec validation:** every feature/`@class` referenced in a spec must exist in `hierarchy.tsv` (else
  config error). At runtime every observed feature-tuple must resolve (via rule or precedence) else error.
  Honors C2 (`novel` allowed; other unknowns terminate).
- **M4 — Overlay semantics:** inputs are assumed to **fully tile** every sequence (no unannotated positions).
  Resolve across the full tiling and keep all positions; **add a coverage check** that each input BED tiles
  each `seq_id` with no gaps — a gap is unexpected → error/warn. (Replaces the current inner-join that
  silently drops singly-annotated positions.)
- **M5 — Implementation:** single interval implementation; **delete the ~5 duplicated pyranges/pandas
  fallbacks**. A sorted-boundary sweep-line per `seq_id` over a full tiling may remove the need for pyranges
  entirely (decide at implementation, favoring simplicity).
- **M6 — Output:** bgzip by default (C3); single write path (remove the two divergent ones).
- **Bugs/cleanup folded in:** V7 (v1 literals → hierarchy-derived classes), add `main()`/importable module,
  delete dead `import pyranges as pr`, coalesce adjacent same-feature intervals consistently, no silent
  fall-through on unknown labels (C2).
- **Proposed home:** `commands/overlay_annotations.py` (thin click) · `core/intervals.py` (pure overlay /
  subtract / coalesce) · `core/annotation_resolution.py` (spec loader + resolver, hierarchy-validated) ·
  `core/io/bed.py` (shared bgzip-aware read/write) · preset specs under package data.

### Review item 3 — `KaryoScope_sequence_annotate` → `build-feature-matrix`

**True scope:** builds a **wide per-`seq_id` feature matrix** (one row per sequence) from already-annotated
BEDs — the direct input to `cluster`. Not an annotator (the engine produces the BEDs). Full metric
reference + proposed column schema: `feature_matrix_metrics.md`.

**Decisions:**
- **F1 — Name: RESOLVED** — subcommand `build-feature-matrix` (active verb, matches `overlay-annotations`).
- **F2 — Column schema: RESOLVED (pending sign-off)** — `{featureset}__{metric}__{feature}` using `__` as the
  sole delimiter (unambiguous; no token contains `__`); `{featureset}__total_bp`, `interspersion__{type}`;
  key column `seq_id`. (Alignment columns move to `cluster-diagnostics` per F6.) **Shared contract** between
  `build-feature-matrix` (producer) and its consumers **`cluster-annotate` + `cluster-diagnostics`** (which
  parse these names) — changed once, together. NOTE: `cluster` does **not** consume this TSV; it builds its
  own edge/abundance adjacency matrix directly from the BEDs (see Review item 4).
- **F3 — Malformed input: RESOLVED** — always **terminate** on malformed BED (bad field count or
  non-integer coords), with line number (no silent skip). Matches C2.
- **F4 — Voodoo constants → parameterize (D6):** `window_size=1000`, `BLOCK_GAP_TOL=100`, adaptive threshold
  `/3` + clamp `[0.1%,5%]` become documented CLI parameters. **Rationale still TBD** (maintainer to supply;
  tracked as an open doc task).
- **F5 — Adaptive-thresholds sidecar:** **IS consumed** by `cluster-annotate` (computes `readpct` =
  % of a cluster's sequences exceeding the per-feature presence threshold). Audit's "advisory only" was
  wrong. Keep the sidecar; document its meaning.
- **F6 — Alignment columns: RESOLVED** — **move the alignment-stats join out of `build-feature-matrix`
  into `cluster-diagnostics`** (its only consumer), keeping the feature matrix clean. When implemented there,
  **genericize** away from the telogator directory layout → a generic per-alignment **stats TSV**
  (documented schema) + optional `seq_id→sequencing_approach` TSV, by explicit path (any aligner works).
- **F7 — total_bp / partition: RESOLVED** — `total_bp` correctness relies on the C4 partition invariant,
  which is now **validated on read** (error on overlap/gap). No double-counting possible.
- **Cleanup folded in:** `read`→`seq_id` (C1); shared `TeeLogger`/BED reader; click defaults (drop
  `_argparse_defaults`); fix `load_stats` brittle `== True`; canonical feature ordering from
  `colors.tsv`/hierarchy (D4.5/D4.2) not per-featureset `.colors.txt`.
- **Proposed home:** `commands/build_feature_matrix.py` · `core/seq_features.py` (pure metrics) ·
  `core/io/bed.py` · `core/io/alignment_stats.py` (generic stats reader).

### Review item 4 — `KaryoScope_cluster_analysis` → split into `build-matrix` + `cluster` + `test-enrichment`

Full rigorous methods audit: **`clustering_methods.md`**.

**Resolved:**
- **CL1 — Structure mode: DROP.** `--analysis-mode structure` is broken (B2, confirmed: `groupby(['read',…])`
  at L355 vs `sequence` column) and bolted-on. Remove it (D5). A coauthor can revive it later from git
  history if needed.
- **CL2 — Decomposition: RESOLVED (pending detailed design).** Split the monolith into three subcommands —
  **`build-matrix`** (BEDs → adjacency matrix), **`cluster`** (matrix → linkage + k-selection + cut →
  assignments), **`test-enrichment`** (assignments + metadata → Fisher/FDR enrichment) — plus the existing
  downstream `cluster-diagnostics`. Decoupling testing from clustering is also the structural fix for the
  circularity (CL-open #2).
- **CL3 — Expose all hidden constants (D6):** every constant in `clustering_methods.md` §5 becomes a
  documented CLI parameter (span window, exclude-features, min/max-k, composite weights, thresholds,
  reduce-dims, silhouette-sample-size, early-stopping, odds cutoff). Length window made explicit (no silent
  default filtering).
- **CL4 — `main()` / importable / no `sys.exit` in helpers; fail loudly (replace bare `except:` at L1927/1946,
  per C2); drop duplicate NPZ `read_names`.**

**OPEN — to work through together (skeptical review; see `clustering_methods.md` §7):**
- **CL-open #1 🔴 Circularity:** k-selection's `composite` mixes silhouette with group enrichment/purity, then
  the same clusters are tested for enrichment → optimistically biased results. Proposal: unsupervised
  k-selection by default; enrichment composite demoted to explicit exploratory; enrichment as a separate step.
- **CL-open #2 🔴 Per-sample stats:** min-p over samples + cluster-only FDR is invalid (max-statistic; min-p
  not a valid p-value). Replace with full sample×cluster FDR or a per-cluster omnibus test.
- **CL-open #3 🟠 Compositional normalization:** round-and-Fisher size-factor hack → proportion/offset model;
  keep/drop `telomeric` vs `total`.
- **CL-open #4 🟠 Three inconsistent "enriched" definitions** (odds>1.5 / raw p<0.05 / FDR q) — unify.
- **CL-open #5 🟠 Default matrix:** z-scoring sparse edge counts inflates rare transitions; reconsider default
  (presence? abundance-only? keep edge block at all?).
- **CL-open #6 🟡 Default k-selector:** `composite-knee` is unstable (depends on min/max-k + early-stop +
  smoothing window) → consider silhouette/Calinski default.
- **CL-open #7 🟠 `--exclude-features` default `canonical_telomere*`** removes the dominant telomere signal
  from clustering by default — intended?

### Review item 5 — `KaryoScope_cluster_annotate` → `cluster-annotate`

**True scope:** Step-2 aggregator. Reads cluster assignments + `cluster_analysis.tsv` + the
`build-feature-matrix` TSV (F2 schema) + adaptive thresholds → **per-cluster annotation TSV**; optional
structural **auto-label**, **representative selection**, and a **feature-importance** PDF. Recomputes nothing
from BEDs (pure aggregation/heuristics).

**Resolved (apply conventions):**
- **CA1 — Name:** `cluster-annotate` (annotates clusters; distinct from the engine's `annotate`).
- **CA2 — Vocab:** drop `_feature_vocab` `sys.path` hack; use hierarchy-derived v2 vocab (D4.1/D4.2); this
  fixes **V1** (`hor→alpha_hor`) and **V8** (0-as-absent) inherited via `lookup_satellite_col`.
- **CA3 — `acrocentric_dmax__rDNA` (L785):** stop hard-coding the featureset name; resolve rDNA's featureset
  from the hierarchy/config; error (not silent 0) if absent (C2).
- **CA4 — De-hardcode biology (D6):** the ~11 auto-label thresholds (L746-756: `CAN_ECTR=70`, `NCAN_ECTR=10`,
  `DMAX_HIGH=35`, `SAT_DOMINANT=80`, `ALT_BLOCK_BP=6000`, …) become documented CLI params / a config preset.
- **CA5 — Consolidate representatives:** the representative-selection code here overlaps
  `select-representatives` → single `core/representatives.py` (see RP-open below).
- **CA6 — Push-downs:** `TeeLogger`→shared `core/io`; `plot_feature_importance` + hard-coded Barthel palette
  (L1082-1096)→`karyoplot.mpl` + `karyoplot.core.colors`; brand rcParams→`karyoplot.mpl.style`.
- **CA7 — Dead code (D5):** `entity_columns` recompute (L1647-1652), convoluted `_fmt`, duplicated read-span
  detection (×3).
- **Home:** `commands/cluster_annotate.py` (thin) · `core/cluster_annotate.py` (aggregation + `auto_label`) ·
  `core/representatives.py` · `core/feature_importance.py` · `core/io/` readers. `auto_label_cluster` is the
  top unit-test target (table-driven per decision-tree branch).

**OPEN:**
- **CA-open #1 🟠 Type I ALT relabel (L1573-1581):** sample-name-specific business logic ("prepend 'Type I
  ALT' if Σ `{sample}_pct` over `--alt-samples` > `--alt-threshold`") embedded in a generic tool, with
  brittle string surgery. Per D6: extract to an explicit, config-driven relabeling step (or fully
  parameterize the label + excluded-label set)? Scientific labeling call — flag for review.
- **CA-open #2:** auto-label thresholds — ship a documented default preset (human/CHM13) and allow override,
  or require explicit config? (Leaning: default preset + override.)

### Review item 6 — `KaryoScope_select_representatives` → `select-representatives`

**True scope:** pick representative sequences per cluster for `cluster-plot --reads-file`, balancing
feature-match and length (length-tiered + feature-group + `_top`-feature modes).

**Resolved:**
- **RP1 — Name:** `select-representatives`.
- **RP2 — B6 fixes:** remove the **unused required `--cluster-analysis`** arg; resolve the
  **`centroid_distance` contradiction** (the column is carried but never drives selection though the name
  implies it) — see RP-open. Fix the `_get_length` arg-order confusion (L47-48).
- **RP3 — Silent failures → loud (C2):** `load_all_features`' bare `except: pass` (L125-126) becomes an
  error/warning; no silent zero-feature results.
- **RP4 — Paths/featuresets (C5):** the hard-coded telogator BED path (L107-110) and featuresets
  `['region','subtelomeric','repeat']` (L412) move to the shared `result_layout` / become configurable.
- **Home:** `commands/select_representatives.py` (thin) · shared `core/representatives.py`.

**OPEN:**
- **RP-open 🟠 Unify representative selection:** `cluster-annotate` (strategies `annotation`|`centroid`) and
  `select-representatives` (length-tiered / feature-group / `_top`) implement **different** algorithms for
  the same goal. Consolidate into one `core/representatives.py` — **which strategy set is authoritative**, and
  should `centroid_distance` actually drive a strategy? Scientific call — flag for review.

### Review item 7 — `KaryoScope_compare_clusterings` → `compare-clusterings`

**True scope:** compare two clusterings of the **same** sequences (inner-join on `seq_id`) → ARI/NMI,
cluster×cluster crosstab, auto-label flow, concordance report + PDF.

**Resolved:**
- **CC1 — Name:** `compare-clusterings`.
- **CC2 — B1 fix (high impact):** `adjusted_rand_index` → **`adjusted_rand_score`** (L84). ARI has **never**
  been computed (the wrong name raises `ImportError`, silently caught) — fix + **narrow the bare
  `except ImportError`** so a real missing-sklearn differs from a typo (C2). A unit test (identical labels →
  ARI=1) would have caught it — add it.
- **CC3 — Dead code (D5):** unused `numpy`/`scipy.stats` imports; unused `generate_report` params; dead
  `annot_matrix` (L509); dead `load_clustering(analysis_file=…)` param; rename misleading
  `plot_enrichment_sankey` (it's a heatmap).
- **CC4 — De-hardcode biology (D6):** the `Pre`/`Post` enrichment color coding (L189-196) is study-specific →
  generalize/parameterize.
- **CC5 — Robust inputs (C2):** replace fragile `str.replace` sidecar discovery (L49,59) with explicit CLI
  args / robust derivation; warn when the two clusterings have disjoint label vocabularies (label
  retention/purity is meaningless otherwise).
- **CC6 — Push-downs:** ARI/NMI → `karyoplot.mpl.statistics` (or core stats); concordance plots → **new**
  `karyoplot.mpl.clustering` (NOT `mpl.comparison` — different semantics); style → `karyoplot.mpl.style`.
- **Home:** `commands/compare_clusterings.py` (thin) · `core/clustering_comparison.py` (pure metrics) · plots
  to `karyoplot.mpl.clustering`.

---

## Tier 4 — plotting (Review items 8–12)

This tier carries the bulk of the **`karyoplot` push-down** (decision #3). Cross-tier consolidation backlog
(extends `README.md` §3): **scale-bar** (4 separate implementations across plot-reads/telogator/animate →
one `karyoplot` helper, converged on `core.coords.pick_round_scale_bp` + `core.text.format_genomic_distance`);
**legend stack** (plot-reads' `compute_legend_layout`/draw/composite + draw-legend's helpers →
`karyoplot.svg.legend`); **PIL/numpy raster** (`_draw_rect_rgba` + raster bar loops → `karyoplot.svg.reads`,
the PIL counterpart to the already-extracted `rasterize_features`); **dendrogram drawing** (→ NEW
`karyoplot.svg.dendrogram`); **matrix/bubble/grid** (→ `karyoplot.svg` + `karyoplot.core.colors` gradient);
**`hex_to_rgb(a)`** (delete dups → `karyoplot.core.colors`); **ffmpeg/video** (→ NEW `karyoplot.video`, D7);
**heatmap tracks** (→ `karyoplot.svg.tracks`/`mpl.heatmap`); **Fisher/FDR** (→ `karyoplot.mpl.statistics`,
for p-value consistency with `cluster`). Also **centralize shared render defaults** (`ratio` 1/300,
`top_margin` 80, `left_margin` 60) so plot-reads/telogator/animate can't drift.

### Review item 8 — `KaryoScope_draw_legend` → `draw-legend`
- **Scope:** legend SVG from a colors file. Rendering already delegated to
  `karyoplot.svg.legend.make_legend_drawing`; ~90 lines of input-shaping helpers remain.
- **Resolved:** name `draw-legend`; `print`→`click.echo`; colors-file reader → `core/io/colors` (or read the
  definitive `colors.tsv`, D4.5); adopt `karyoplot` `DEFAULT_THEME` (drop the back-compat inline `Theme`).
- **OPEN (minor):** push the 4 helpers (`load_colors_file`, `filter_items`, `parse_groups`, `group_items`)
  **into `karyoplot.svg.legend`** (next to `merge_by_color`/`strip_label_suffixes`) so the subcommand becomes
  trivial, vs keep them in analysis `core/legend.py`? (Lean: push to karyoplot.) Resolve the include/exclude
  exact-match vs `--groups` prefix-match asymmetry.

### Review item 9 — `KaryoScope_plot_reads` → `plot-reads`
- **Scope:** stacked per-sequence feature-bar figures (vertical/horizontal, SVG/PNG) + heatmap track +
  two-tier labels + auto-legend + animation hook. One of two big read renderers (shares the rasterization
  stack with `cluster-plot`).
- **Resolved:** name `plot-reads`. Decompose: `commands/plot_reads.py` (thin) · `core/plot_reads/`
  (`build_render_config`, `prepare_reads`, `apply_viewport_ratio` from the ~570-line `main()`) ·
  `core/io/reads_bed.py` (loaders) · `core/orientation.py`. Push-downs per the backlog above (legend is the
  single biggest extraction). `read`→`seq_id` (C1). Hardcoded human chrom/satellite/telomere vocab
  (L2334-2355) → shared `core/feature_vocab` (hierarchy/DB-driven, D4). Hardcoded telogator path → C5.
  Replace the `sys.path.insert` animation import with a real import (ties to item 12).
- **Dead code (D5):** `--animate-crop-ratio` (parsed, never read); `font_size=14` fallbacks (always 11);
  `compute_viewport_params` (L2284, unused); `_estimate_heatmap_legend_height` (unused); dead `draw_ctx`
  reassign; triplicated `compute_group_spans`/`_label`.
- **OPEN:** merge `--max-length` vs `--max-read-length` into one? · unify the legend with `karyoplot` vs keep
  the section-aware grid local? · converge scale-bar onto karyoplot helpers (changes SVG bytes — fine, D1).

### Review item 10 — `KaryoScope_cluster_plot` → `cluster-plot`
- **Scope:** the cluster-representative figure (per-read feature bars + dendrograms + sample×cluster matrices
  + enrichment bubbles/grids + legends), horizontal & vertical layouts. 7,213 lines; `main()` ≈ 1,490 with
  **two near-duplicate layout engines**.
- **Resolved:** name `cluster-plot`. **MUST split** → `commands/cluster_plot.py` (thin) ·
  `core/cluster_plot/{io,selection,dendrogram,colors,density,layout,render_horizontal,render_vertical}.py` ·
  `core/cluster_plot/draw/` (thin over karyoplot). Push-downs: dendrogram drawing → NEW
  `karyoplot.svg.dendrogram`; bubbles/grids/legends → `karyoplot.svg.legend`; matrix gradient →
  `karyoplot.core.colors` + svg matrix drawer; **Fisher/FDR → `karyoplot.mpl.statistics`** (consistency).
  Keep the cleaner inline above-cut subtree extraction; **drop the Bio.Phylo Newick round-trip** (two paths,
  one problem). Write derived `*.FIRE_LINKER.bed` to output/temp, **not** the input dir (side-effect in a
  plot command). Hardcoded path → C5.
- **CP-structural — DROP:** the auto-detected **structural mode** plots `cluster`'s structure-mode output,
  which we dropped (CL1) → drop `plot_structural_mode` too (D5).
- **Dead code / unimplemented flags (D5):** unused `import subprocess`; likely-dead feature-scoring selection
  functions (`parse_top_features`, `score_read_features`, `_select_by_strategy`, …) — **confirm not
  externally called** before deleting; flags parsed-but-never-used (`--enrichment-normalization
  telomeric/total`, `--total-reads-file`, `--orient-telomere-top`, `--show-clade-id/-count`) → remove (or
  implement). `--min-feature-width` default 0.5 vs help 1.0 — pick one.
- **OPEN:** unify the two layout engines vs keep separate (real effort)? · scope of karyoplot drawing
  push-down (dendrogram/bubble/matrix are the strongest candidates) · confirm the dead selection functions.

### Review item 11 — `KaryoScope_telogator_reads_viz` → generalize (the over-specific-name case)
- **True scope:** renders **any** per-sequence feature BED as side-by-side **vertical feature bars** (SVG +
  PNG), sorted by length, grouped by sample. The "telogator" framing is incidental (C1 / D3) — input is just
  a feature BED. It is a **simpler subset of `plot-reads`** (which already does vertical/horizontal feature
  bars, with more features).
- **Resolved:** drop the telogator-specific framing; `read`→`seq_id`; delete local `hex_to_rgb` (use
  `karyoplot.core.colors`); collapse the duplicated BED loader (2 copies) and the duplicated SVG/PNG renderer;
  consolidate its 2 scale-bars into the shared karyoplot helper; set `Image.MAX_IMAGE_PIXELS=None`; **fail
  loudly** on missing input / no reads (currently exits 0 — violates C2); hardcoded telogator path → C5.
- **OPEN 🟠 (scoping):** is this a **distinct tool or redundant with `plot-reads`?** Both render per-sequence
  feature bars. Options: **(a)** merge into `plot-reads` as a vertical preset (and reuse
  `karyoplot.svg.reads`); **(b)** keep a separate, generally-named subcommand that shares the renderer. Needs
  your call. If kept separate, pick a general name (not "telogator").

### Review item 12 — `create_panning_animation` → `animate`
- **Scope:** PNG → panning MP4 (fixed-zoom loop + adaptive content-aware zoom). Already a de-facto library
  (`plot-reads` imports it via `sys.path` hack).
- **Resolved:** name `animate`. Decompose → `commands/animate.py` (thin) ·
  `core/animation/{profile,adaptive,fixed,encoder}.py`. Replace the `plot-reads` `sys.path` import with a real
  package import. Push-downs: `_load_font` → delete (`karyoplot.core.fonts.pil_font`);
  `svg_to_png`/`ensure_png`/`get_svg_dimensions` → `karyoplot.svg.export`; **ffmpeg encoder → NEW
  `karyoplot.video`** (D7), de-duplicating the 4 inline encoder arg-lists; scale-bar → shared helper. Add an
  ffmpeg/rsvg **presence check** (fail loudly, C2); `.convert("RGB")` the main input; fix temp-file leak in
  fixed mode.
- **Dead code (D5):** `compute_strip_zoom_levels` (unused); unused imports `math`/`io`; `--uniform-zoom` is a
  no-op; `--crop-ratio` ignored in adaptive mode.
- **OPEN:** is **fixed-zoom mode** still used, or has **adaptive** superseded it? Dropping fixed-zoom would
  remove most of the duplication (`create_horizontal/vertical_panning`, `--crop-ratio`,
  `--vertical-only-zoom`, static `--scale-bar`). Centralize the shared `ratio`/margin defaults (vs telogator /
  plot-reads).

---

## Tier 5 — translocation trio (Review items 13–15)

Three tools fed by the engine's `*.translocations.bed.gz`. **`find` and `cluster` are siblings** (each
independently re-globs the results dir), **not** a linear pipe; `find`'s TSV feeds only `visualize`'s
`--input-tsv` mode. None do breakpoint detection — translocation reads are detected upstream by the engine;
these only catalog/cluster/draw them. Shared across all three: `TeeLogger` (→ `core/logging`), the
discovery regex (→ C5 `result_layout`), and `seq_id` (C1). **File-naming wrinkle to reconcile in
`result_layout`:** stages 1&2 read `*.presmoothed.translocations.bed.gz`; stage 3 reads
`*.smoothed.…translocations.bed.gz` (TSV mode) or `*.smoothed.features.bed.gz` (direct mode).

### Review item 13 — `KaryoScope_find_translocation_reads` → `find-translocation-reads`
- **Scope:** discover chromosome translocation BEDs → per-`seq_id` length + per-target-chrom bp/pct coverage
  TSV. Stdlib-only; lightest of the trio.
- **Resolved:** name `find-translocation-reads`; discovery regex → C5; BED read → `karyoplot.core.io`;
  `TeeLogger` → shared; `seq_id` (C1); fail loudly (C2).
- **OPEN:** drop the `_bp`/`_pct` columns (written but unused downstream)? · widen `chr\d+_chr\d+` to allow
  `chrX/chrY` translocations? · read length from an authoritative source vs `max(end)` inference? · merge
  overlapping intervals before summing coverage (today `pct` can exceed 100)?

### Review item 14 — `KaryoScope_cluster_translocation_reads` → `cluster-translocation-reads`
- **Scope:** orchestrator — discover/group/merge translocation BEDs per group, then run clustering + plotting.
  Today **shells out** to `merge_beds` / `cluster_analysis` / `cluster_plot` via `--scripts-dir`.
- **Resolved:** name `cluster-translocation-reads`. **CRITICAL: replace the 3 subprocess calls with
  in-process calls** to the migrated `overlay-annotations` + (`build-matrix`→`cluster`→`test-enrichment`) +
  `cluster-plot` — drops `--scripts-dir`. Discovery regex → C5; `concatenate_beds`/`create_samples_tsv` →
  `core/io` + `karyoplot.core.sample_metadata`; `TeeLogger` → shared; `seq_id`. The frozen cluster recipe
  (~12 hard-coded flags) and hard-coded `--control-group primary` → exposed defaults/`--group-config` (D6).
- **OPEN:** make the frozen recipe configurable (vs fixed for reproducibility)? · is a group literally named
  `primary` a hard requirement? · in-process vs keep subprocess isolation (matplotlib backend/memory)? ·
  should it consume `find`'s TSV (confirmed reads) instead of re-globbing? · is the `≥2 samples` gate right?

### Review item 15 — `KaryoScope_visualize_translocation_reads` → `visualize-translocation-reads`
- **Scope:** the only drawing tool of the trio — per-read stacked multi-featureset SVG/PNG (TSV-batch or
  direct-spec). Already uses `karyoplot.core.text`/`colors`.
- **Resolved:** name `visualize-translocation-reads`. **Fix B4** (missing `rsvg-convert` → SVG deleted in
  `finally` → silent no-output): switch to `karyoplot.svg.export.svg_to_png` (keeps SVG, graceful, C2).
  Push-downs (high value): `compute_scale_bar_bp` → `karyoplot.core.coords.pick_round_scale_bp`;
  `draw_scale_bar` → karyoplot; `draw_legend` → `karyoplot.svg.legend`; BED → `karyoplot.core.io`;
  `natural_sort_key` → `karyoplot.core.chromosomes`. `seq_id` (C1). Fix `featureset_colors[fs]` KeyError (use
  `.get`).
- **De-hardcode biology (D6):** the red/blue label color keyed on `"chr2_chr13"` (L403) → from a
  translocation palette/config; hard-coded `font_family="Basic Sans"` → `--font-family`.
- **OPEN:** reconcile the **two BED naming conventions** + two length sources between modes (and is the
  transposed `chromosome.smoothed.{trans}` name intentional?) · is the `:`-delimited read-spec (breaks on
  read IDs containing `:`) acceptable?
