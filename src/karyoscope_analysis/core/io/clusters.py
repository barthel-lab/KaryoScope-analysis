"""Readers for Engine B clustering outputs (``layout.tsv``).

Shared so the clustering-downstream tools (``test-enrichment``, ``compare-clusterings``, …) parse
``cluster``'s sidecars in exactly one place.
"""

from __future__ import annotations

from pathlib import Path


def read_layout_assignments(path: str | Path) -> dict[str, str]:
    """Map ``read_id -> cluster_id`` from a ``layout.tsv`` (all of a read's segments share it)."""
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"empty layout file: {path}")
    header = lines[0].split("\t")
    try:
        ci, ri = header.index("cluster_id"), header.index("read_id")
    except ValueError as e:
        raise ValueError(
            f"{path}: expected 'cluster_id' and 'read_id' columns, got {header}"
        ) from e
    out: dict[str, str] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) > max(ci, ri):
            out[fields[ri]] = fields[ci]
    return out
