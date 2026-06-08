"""Top-level command-line interface for KaryoScope-analysis.

This module wires the analysis subcommands into a single ``karyoscope-analysis``
entry point. Each subcommand lives in its own module under
``karyoscope_analysis.commands`` so that they can grow independently, mirroring
the structure of the core ``karyoscope`` engine CLI.

Logging vs. program output
==========================

Two channels coexist deliberately:

* **Program output** — what a command is doing for the user. Emitted via
  ``click.echo``. Always visible regardless of verbosity.
* **Logging / diagnostics** — behind-the-scenes information for developers and
  power users debugging issues. Emitted via the standard ``logging`` module.
  Hidden by default; opt in with ``-v`` or ``-vv``.

The verbosity flags only affect the logging channel.
"""

from __future__ import annotations

import logging
import sys

import click

from karyoscope_analysis._version import __version__
from karyoscope_analysis.commands import (
    bin_annotations,
    build_feature_matrix,
    cluster,
    cluster_plot,
    detect_rearrangements,
    genome_weights,
    overlay_annotations,
    version,
)

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "max_content_width": 100,
}

#: Maps a verbosity integer (negative for quiet, 0 default, positive for verbose)
#: to a stdlib logging level. Anything beyond ``2`` is clamped to DEBUG.
_VERBOSITY_TO_LEVEL = {
    -1: logging.ERROR,
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def _configure_logging(verbosity: int) -> None:
    """Install a stderr log handler with a format suited to a CLI tool.

    Replaces any pre-existing handlers on the root logger so that repeated
    calls (e.g., from tests) don't compound output.
    """
    level = _VERBOSITY_TO_LEVEL.get(verbosity, logging.DEBUG if verbosity > 0 else logging.ERROR)
    handler = logging.StreamHandler(stream=sys.stderr)
    if level <= logging.DEBUG:
        fmt = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
    elif level <= logging.INFO:
        fmt = "%(asctime)s %(levelname)s: %(message)s"
    else:
        fmt = "%(levelname)s: %(message)s"
    handler.setFormatter(logging.Formatter(fmt, datefmt="%H:%M:%S"))
    root = logging.getLogger()
    root.handlers = [handler]
    root.setLevel(level)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    __version__,
    "-V",
    "--version",
    message="karyoscope-analysis %(version)s",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase logging verbosity. Repeat for more (-v=info, -vv=debug).",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Decrease logging verbosity to errors only. Conflicts with -v.",
)
def main(verbose: int, quiet: bool) -> None:
    """KaryoScope-analysis: cluster, annotate, and visualize KaryoScope sequence annotations.

    Run a subcommand with ``--help`` to see its options, e.g.:

    \b
        karyoscope-analysis version

    For more information, see https://github.com/barthel-lab/KaryoScope-analysis.
    """
    if quiet and verbose:
        raise click.UsageError("--quiet and --verbose cannot be combined.")
    verbosity = -1 if quiet else verbose
    _configure_logging(verbosity)


# Register subcommands. The order here determines the order in `--help`.
main.add_command(bin_annotations.cmd, name="bin-annotations")
main.add_command(overlay_annotations.cmd, name="overlay-annotations")
main.add_command(build_feature_matrix.cmd, name="build-feature-matrix")
main.add_command(detect_rearrangements.cmd, name="detect-rearrangements")
main.add_command(genome_weights.cmd, name="genome-weights")
main.add_command(cluster.cmd, name="cluster")
main.add_command(cluster_plot.cmd, name="cluster-plot")
main.add_command(version.cmd, name="version")

# --- Roadmap (Phase 4 migration; see docs/audit/DECISIONS.md) ---
# Data foundation:   bin-annotations ✓ (hierarchy-aware mode filter; denoise pre-overlay),
#                    overlay-annotations ✓, build-feature-matrix ✓
# Rearrangements:    detect-rearrangements ✓ (Engine A; differential colocalization)
# Clustering:        cluster ✓ (Engine B; OLC clustering + consensus),
#                    genome-weights ✓ (reference-genome information-content feature weights)
# Plotting:          cluster-plot ✓ (read-renderer; SVG). Deferred: animation/video (D7),
#                    Engine A bubble/matrix, karyoplot.svg push-down.
# Translocations:    find-/cluster-/visualize-translocation-reads


if __name__ == "__main__":
    main()
