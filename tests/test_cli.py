"""Smoke tests for the karyoscope-analysis CLI skeleton.

These exercise the click group wiring (version, help, verbosity handling) so the
entry point can't silently break as subcommands are added during migration.
"""

from __future__ import annotations

from karyoscope_analysis import __version__
from karyoscope_analysis.cli import main


def test_version_flag(cli_runner):
    result = cli_runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help(cli_runner):
    result = cli_runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "version" in result.output  # the registered subcommand is listed


def test_version_subcommand(cli_runner):
    result = cli_runner.invoke(main, ["version"])
    assert result.exit_code == 0
    assert __version__ in result.output
    assert "Python dependencies" in result.output


def test_quiet_and_verbose_conflict(cli_runner):
    result = cli_runner.invoke(main, ["-q", "-v", "version"])
    assert result.exit_code != 0
