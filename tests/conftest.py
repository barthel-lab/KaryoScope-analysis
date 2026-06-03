"""Shared pytest fixtures for the karyoscope-analysis test suite."""

from __future__ import annotations

import sys

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


@pytest.fixture
def cli_runner() -> CliRunner:
    """A click CliRunner for invoking subcommands in tests."""
    return CliRunner()
