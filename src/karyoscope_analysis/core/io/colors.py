"""Read the database ``colors.tsv`` (the v2 feature palette).

``colors.tsv`` is ``feature_set <TAB> feature <TAB> color`` (decision D4.5 / D4.2).
Parsing is delegated to :func:`karyoscope.core.io.colors.parse_colors`, the canonical
loader in the KaryoScope engine, so the three-column format is parsed and validated in
exactly one place across the ecosystem (no duplicate parser here, in ``karyoplot``, or in
the engine). The path is resolved by the caller (eventually via ``karyoscope.paths``).

For rendering we usually want a ``feature -> color`` map collapsed across featuresets
(feature names are effectively unique across featuresets in the v2 vocabulary; on the rare
collision the row from the later featureset block wins, matching the file order).
"""

from __future__ import annotations

from pathlib import Path

from karyoscope.core.io.colors import parse_colors


def load_colors_by_featureset(path: str | Path) -> dict[str, dict[str, str]]:
    """Read ``colors.tsv`` into the canonical ``{feature_set: {feature: hex}}`` mapping.

    Thin pass-through to :func:`karyoscope.core.io.colors.parse_colors`; kept here so
    analysis code has a single import site for DB colors and callers needing the
    featureset-aware structure (e.g. legend grouping) don't reach across packages.
    """
    return parse_colors(Path(path))


def load_colors(path: str | Path) -> dict[str, str]:
    """Read ``colors.tsv`` into ``{feature: hex_color}`` (collapsed across featuresets).

    Featuresets appear in contiguous blocks in the file and :func:`parse_colors` preserves
    that order, so updating in featureset order reproduces the file's last-row-wins behaviour
    on the rare cross-featureset name collision.
    """
    colors: dict[str, str] = {}
    for fs_colors in load_colors_by_featureset(path).values():
        colors.update(fs_colors)
    return colors
