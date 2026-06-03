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

### Notes

- Migration of the analysis scripts into the package (with bug fixes, v2-only
  feature vocabulary, and `karyoplot` push-down) is in progress; see
  `docs/audit/DECISIONS.md`.
