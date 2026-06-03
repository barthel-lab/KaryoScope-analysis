"""``karyoscope-analysis version`` — print version and environment info.

This is the detailed form of ``karyoscope-analysis --version``. It reports the
package version, the Python interpreter, the install location, the presence and
versions of key Python dependencies and the sibling KaryoScope-ecosystem packages
(``karyoplot``, ``karyoscope``), and the external command-line tools used for
rendering (``rsvg-convert``, ``ffmpeg``). The output is suitable for pasting into
bug reports.
"""

from __future__ import annotations

import platform
import shutil
import sys
from importlib.metadata import PackageNotFoundError, version

import click

from karyoscope_analysis._version import __version__

#: Python dependencies whose versions we report.
_PYTHON_DEPS: tuple[str, ...] = (
    "click",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "drawsvg",
    "pillow",
    "pyyaml",
    "jsonschema",
)

#: Sibling KaryoScope-ecosystem packages (installed editable; see CONTRIBUTING.md).
_ECOSYSTEM_PACKAGES: tuple[str, ...] = (
    "karyoplot",
    "karyoscope",
)

#: External command-line tools used for rendering.
_EXTERNAL_TOOLS: tuple[str, ...] = (
    "rsvg-convert",
    "ffmpeg",
)


def _dep_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


@click.command(
    help="Print karyoscope-analysis version and environment info (useful for bug reports).",
    no_args_is_help=False,
)
def cmd() -> None:
    """Show detailed version and environment information."""
    click.echo(f"karyoscope-analysis {__version__}")
    click.echo(f"  Python: {sys.version.split()[0]} ({sys.executable})")
    click.echo(f"  Platform: {platform.platform()}")

    click.echo("\nPython dependencies:")
    for name in _PYTHON_DEPS:
        click.echo(f"  {name}: {_dep_version(name)}")

    click.echo("\nKaryoScope-ecosystem packages:")
    for name in _ECOSYSTEM_PACKAGES:
        click.echo(f"  {name}: {_dep_version(name)}")

    click.echo("\nExternal tools:")
    for tool in _EXTERNAL_TOOLS:
        path = shutil.which(tool)
        click.echo(f"  {tool}: {path if path else 'not found on PATH'}")
