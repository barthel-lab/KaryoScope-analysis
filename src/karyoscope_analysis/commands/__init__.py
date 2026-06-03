"""KaryoScope-analysis subcommands.

Each module in this package contributes one click command, exported as ``cmd``.
The :mod:`karyoscope_analysis.cli` module aggregates them into the top-level
``karyoscope-analysis`` group. Command modules should stay thin: parse options,
then delegate to pure logic in :mod:`karyoscope_analysis.core`.
"""
