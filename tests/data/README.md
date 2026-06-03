# Test fixtures

Small, committed fixtures for the test suite (mirrors the gold-standard
`KaryoScope/tests/data/` convention). Keep everything here tiny.

Planned (added during Phase 4/5 migration, see `docs/audit/DECISIONS.md`):

- A minimal **`KS_human_CHM13_v2`** example dataset (a handful of sequences across
  the six featuresets) to drive unit + golden tests, since the package is
  standardizing on v2 (decision D4.1). The legacy v1 example data currently under
  the repo's top-level `data/` will be replaced by these v2 fixtures.
- Tiny annotation BEDs, a sample-metadata TSV, and a feature-matrix `.npz` for the
  clustering/annotation tools.

Fixtures are intentionally checked in (see the `!tests/data/**` rule in
`.gitignore`).
