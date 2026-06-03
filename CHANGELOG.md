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
  `.gz` (gzip/bgzip) read transparently.
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
- **`overlay-annotations` subcommand** (replaces `KaryoScope_merge_beds.py`): reads one
  C4-validated annotation BED per featureset, refines the tracks per `seq_id`, resolves
  each segment via a preset / custom spec / the default basic overlay, and writes one
  resolved annotation BED. Validates input features against the hierarchy (C2). First
  fully-migrated tool, wired into the `karyoscope-analysis` CLI.
- `core/seq_features.py`: pure per-`seq_id` feature metrics for `build-feature-matrix`
  — coverage (`bp`/`frac`/`total_bp`), 1-bp-step sliding-window density
  (`dmax`/…/`dterminal`), `max_block_bp` (gap-bridged), hierarchy-derived
  interspersion, and the adaptive-threshold computation. Every magic constant
  (window size, block-gap tolerance, threshold factor/bounds) is a parameter (F4).
- **`build-feature-matrix` subcommand** (replaces the matrix-building part of
  `KaryoScope_sequence_annotate.py`): reads one annotation BED per featureset and
  emits the wide per-`seq_id` matrix (`{featureset}__{metric}__{feature}` schema, F2)
  + an adaptive-threshold sidecar (F5). Constants are CLI options (F4); alignment-QC
  columns intentionally omitted (they move to `cluster-diagnostics`, F6). Second
  fully-migrated tool; completes the data-foundation tier.

### Notes

- Migration of the analysis scripts into the package (with bug fixes, v2-only
  feature vocabulary, and `karyoplot` push-down) is in progress; see
  `docs/audit/DECISIONS.md`.
