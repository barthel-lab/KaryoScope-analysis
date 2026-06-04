"""Read the database ``colors.tsv`` (the v2 feature palette).

``colors.tsv`` is ``feature_set <TAB> feature <TAB> color`` (decision D4.5 / D4.2). For
rendering we usually want a feature -> color map collapsed across featuresets (feature names
are effectively unique across featuresets in the v2 vocabulary; on the rare collision the
last row wins). The path is resolved by the caller (eventually via ``karyoscope.paths``);
this stays a pure file loader.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def _open_text(path: Path) -> Iterator:
    fh = gzip.open(path, "rt") if path.suffix == ".gz" else path.open()  # noqa: SIM115 (closed below)
    try:
        yield fh
    finally:
        fh.close()


def load_colors(path: str | Path) -> dict[str, str]:
    """Read ``colors.tsv`` into ``{feature: hex_color}`` (collapsed across featuresets)."""
    path = Path(path)
    colors: dict[str, str] = {}
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3 or fields[0] == "feature_set":  # skip header
                continue
            _feature_set, feature, color = fields[0], fields[1], fields[2]
            colors[feature] = color
    return colors
