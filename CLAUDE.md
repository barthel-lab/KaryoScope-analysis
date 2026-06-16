# KaryoScope-analysis — orientation for Claude

Reorganizing this repo (originally a flat `scripts/` of legacy `KaryoScope_*.py` tools) into an
installable package, `karyoscope_analysis`, exposing a unified `karyoscope-analysis` CLI. Modeled on
the sibling `KaryoScope` engine and `KaryoScope-plotlib` (`karyoplot`) library.

## Branch & workflow
- All work is on **`package-reorg`** (off `main`; `main` and the legacy top-level `scripts/` are left
  untouched). Pushed to `origin/package-reorg` (`git@github.com:barthel-lab/KaryoScope-analysis.git`).
- **Commit and push only when the user asks.** If on `main`, branch first.
- v2-only: the vocab/DB is the definitive `KaryoScope-databases/KS_human_CHM13_v2`
  (`hierarchy.tsv`/`colors.tsv`); v1 was dropped.

## Setup
```bash
git checkout package-reorg
python3 -m venv .venv
# Sibling packages (not on PyPI) must be installed editable FIRST — `karyoscope`
# provides the canonical DB parsers (core.io.{colors,hierarchy}); `karyoplot` the renderers.
.venv/bin/pip install -e ../KaryoScope -e ../KaryoScope-plotlib
.venv/bin/pip install -e '.[dev]'          # base + pytest; then: .venv/bin/pip install ruff
.venv/bin/python -m pytest -q              # ~227 pass; integration tests deselected by default
.venv/bin/ruff check src tests             # line-length 100 (config in pyproject.toml)
```
Integration tests (marked `@pytest.mark.integration`) run on the full committed BEDs and are
deselected by default via `addopts`; run them with `-m integration`.

## Read these first (the durable design record)
- `docs/audit/DECISIONS.md` — conventions (C1–C5), cross-cutting decisions (D1–D8), per-script review.
- `docs/audit/OPEN_QUESTIONS.md` — what still needs maintainer/coauthor input (esp. **Group A** =
  Engine A clustering statistics).
- `docs/audit/rearrangement_detection.md` — the authoritative Engine A/B model, spec, and validation;
  **§11–13** cover `bin-annotations`, `genome-weights`, and the **whole-sample RUNBOOK (§13)**.
- `docs/audit/{README,KNOWN_ISSUES,clustering_methods,feature_matrix_metrics}.md`.
- `CHANGELOG.md` `[Unreleased]`, `CONTRIBUTING.md`, and the roadmap comment in
  `src/karyoscope_analysis/cli.py` (~line 121).

## CLI today
`bin-annotations`, `overlay-annotations`, `build-feature-matrix`, `detect-rearrangements`,
`genome-weights`, `cluster`, `cluster-plot`, `version`. The whole-sample clustering pipeline is
wrapped as one command: `scripts/run_cluster_pipeline.sh --sample S --prefix P --db DB`
(bin → overlay → cluster → cluster-plot; see §13).

## Data (what's in git vs not)
- **Committed** (so it travels with the repo): the telogator featureset BEDs
  `data/raw_bed/{HeLa,IMR90,U2OS}.telogator.1.KS_human_CHM13(_v2).{region,subtelomeric,chromosome,
  gene,repeat,acrocentric}.smoothed.features.bed.gz`; the genome weights
  `data/chm13v2_feature_weights.tsv`; small test fixtures in `tests/data/`.
- **NOT in git** (copy to the machine only to *run the pipeline*, not to develop/test):
  the DB `hierarchy.tsv`/`colors.tsv` (from the separate `KaryoScope-databases/KS_human_CHM13_v2`
  repo); the 475 MB CHM13 reference BEDs `data/raw_bed/chm13v2.0.*` (only needed to *recompute*
  genome weights — the weights are already committed).
- `.venv/` and `plots_preview/` are gitignored. `pytest` needs none of the external data.

## State: done / next
- **Done:** package scaffold (hatchling/ruff/pytest/CI); the data-foundation tier
  (`core/feature_vocab`, `core/intervals`, `core/io/bed`, `overlay-annotations`,
  `build-feature-matrix`, `bin-annotations`); Engine A (`detect-rearrangements`) and Engine B
  (`cluster` — chromosome-aware overlap-layout-consensus) + `genome-weights`; and the
  consensus-coordinate layout / `cluster-plot` (extensively refined: backbone selection,
  read orientation, concordance refinement). 227 tests pass, ruff clean.
- **Next:** broaden tests (unit/integration/golden); update docs for the new structure; and the open
  items in `docs/audit/OPEN_QUESTIONS.md` (Engine A Group-A statistics sign-off; read-dedup; the
  deferred Engine B `cluster_5` under-split). The legacy `scripts/KaryoScope_*.py` not yet folded in
  are tracked in `DECISIONS.md` / the `cli.py` roadmap.
