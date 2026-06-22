# Contributing to KaryoScope-analysis

Thanks for your interest in improving KaryoScope-analysis! This package provides
the clustering, annotation, and visualization tools that sit on top of the core
[KaryoScope](https://github.com/barthel-lab/KaryoScope) engine and the shared
[KaryoScope-plotlib](https://github.com/barthel-lab/KaryoScope-plotlib)
(`karyoplot`) plotting library.

## Development setup

This package depends on two sibling KaryoScope-ecosystem packages that are **not
on PyPI**: `karyoplot` (KaryoScope-plotlib) and `karyoscope` (KaryoScope). Install
them editable from their sibling checkouts first, then this package:

```bash
# From the parent directory that holds all three repos:
python3 -m venv KaryoScope-analysis/.venv
source KaryoScope-analysis/.venv/bin/activate
python -m pip install --upgrade pip

# Sibling packages (editable):
pip install -e KaryoScope-plotlib        # provides `karyoplot`
pip install -e KaryoScope                # provides `karyoscope`

# This package (with dev + docs extras):
pip install -e "KaryoScope-analysis[dev,docs]"

# Enable the lint/format hooks:
cd KaryoScope-analysis && pre-commit install
```

The package targets **Python 3.10+** and uses a `src/` layout
(`src/karyoscope_analysis/`), matching the core engine.

## Running the tests

```bash
python -m pytest                 # unit tests (integration deselected by default)
python -m pytest --cov=karyoscope_analysis --cov-report=term-missing
python -m pytest -m integration  # tests needing sibling pkgs / ffmpeg / rsvg-convert
```

## Linting and formatting

Ruff (lint + format) is pinned in `.pre-commit-config.yaml` — the single source of
truth, used by both developers and CI:

```bash
pre-commit run --all-files
```

## Architecture

- `commands/` — thin click subcommands (option parsing only).
- `core/` — pure, importable logic (clustering, enrichment, intervals, feature
  matrices, representative selection, …).
- `core/io/` — readers/writers (bgzip-aware BED, result-layout paths, colors,
  hierarchy, clustering artifacts).
- Shared rendering/IO that belongs to the whole ecosystem is pushed **down** into
  `karyoplot`.

The **design record** for the in-progress reorganization lives in
[`docs/audit/`](docs/audit/) — start with `docs/audit/README.md` and
`docs/audit/DECISIONS.md`. Open scoping/statistics questions are tracked in
`docs/audit/OPEN_QUESTIONS.md`.

## Pull requests

- Branch off `main`; keep PRs focused.
- Make sure `pre-commit run --all-files` and `pytest` pass.
- Add or update tests for behavior changes; for figure-producing code, prefer a
  golden-output test.
- Note user-facing changes in `CHANGELOG.md` under `[Unreleased]`.

By contributing, you agree that your contributions are licensed under the
project's GPL-3.0-or-later license.
