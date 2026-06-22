"""Shared pytest fixtures for the karyoscope-analysis test suite."""

from __future__ import annotations

import sys
from pathlib import Path

# Friendly diagnostic for a common environment mistake — especially on macOS
# where multiple Python installations coexist. If pytest runs from a different
# Python than the one where `pip install -e .[dev]` was run, show this instead
# of a wall of ModuleNotFoundError tracebacks. Must run BEFORE the imports below.
try:
    import karyoscope_analysis  # noqa: F401
except ImportError as _import_err:  # pragma: no cover
    _msg = (
        "\n"
        "    karyoscope_analysis is not importable from this Python interpreter.\n"
        "\n"
        f"    Python being used:  {sys.executable}\n"
        "\n"
        "    This usually means pytest is running from a different Python\n"
        "    environment than the one where 'pip install -e .[dev]' was run.\n"
        "\n"
        "    To fix, either:\n"
        "      (a) Run pytest through the right Python (recommended):\n"
        "              python -m pytest\n"
        "      (b) Or install the package in this Python:\n"
        f"              {sys.executable} -m pip install -e '.[dev]'\n"
    )
    raise RuntimeError(_msg) from _import_err

import pytest
from click.testing import CliRunner

#: ``tests/`` directory and its committed data subdirectories.
TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
V2_SUBSET_DIR = DATA_DIR / "v2_subset"
#: Repo root and the full (large, not part of the default test run) raw BEDs.
REPO_ROOT = TESTS_DIR.parent
RAW_BED_DIR = REPO_ROOT / "data" / "raw_bed"

#: The six KaryoScope featuresets, in the canonical order.
FEATURESETS = ("acrocentric", "chromosome", "gene", "region", "repeat", "subtelomeric")


@pytest.fixture
def cli_runner() -> CliRunner:
    """A click CliRunner for invoking subcommands in tests."""
    return CliRunner()


@pytest.fixture(scope="session")
def hierarchy_tsv() -> Path:
    """Path to the committed v2 hierarchy.tsv fixture."""
    return DATA_DIR / "hierarchy.tsv"


@pytest.fixture(scope="session")
def v2_subset_beds() -> dict[str, Path]:
    """``{featureset: path}`` for the tiny committed v2 HeLa subset fixtures."""
    return {fs: V2_SUBSET_DIR / f"HeLa.v2.{fs}.bed.gz" for fs in FEATURESETS}


def raw_bed_paths(sample: str) -> dict[str, Path] | None:
    """``{featureset: path}`` for a sample's full v2 raw BEDs, or ``None`` if absent.

    The full BEDs are large and not required for the default test run; integration
    tests use this and skip when the data isn't present locally.
    """
    paths = {
        fs: RAW_BED_DIR / f"{sample}.telogator.1.KS_human_CHM13_v2.{fs}.smoothed.features.bed.gz"
        for fs in FEATURESETS
    }
    return paths if all(p.is_file() for p in paths.values()) else None


@pytest.fixture(scope="session")
def raw_bed_lookup():
    """Return the :func:`raw_bed_paths` helper (returns ``None`` when data is absent)."""
    return raw_bed_paths
