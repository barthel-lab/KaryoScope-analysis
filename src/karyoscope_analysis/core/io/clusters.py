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


def _header_and_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"empty file: {path}")
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:] if line]


def read_clusters_table(path: str | Path) -> tuple[dict[str, int], dict[str, int]]:
    """Read ``clusters.tsv`` into ``({cluster_id: size}, {cluster_id: width})``."""
    path = Path(path)
    header, rows = _header_and_rows(path)
    try:
        ci, si, wi = header.index("cluster_id"), header.index("size"), header.index("width")
    except ValueError as e:
        raise ValueError(
            f"{path}: expected 'cluster_id', 'size', 'width' columns, got {header}"
        ) from e
    sizes: dict[str, int] = {}
    widths: dict[str, int] = {}
    for r in rows:
        if len(r) > max(ci, si, wi):
            sizes[r[ci]] = int(r[si])
            widths[r[ci]] = int(r[wi])
    return sizes, widths


def read_consensus_segments(path: str | Path) -> dict[str, list[tuple[int, int, str]]]:
    """Read ``consensus.bed`` into ``{cluster_id: [(start, end, feature), ...]}``."""
    path = Path(path)
    header, rows = _header_and_rows(path)
    try:
        ci, sti, eni, fti = (
            header.index("cluster_id"),
            header.index("start"),
            header.index("end"),
            header.index("feature"),
        )
    except ValueError as e:
        raise ValueError(
            f"{path}: expected 'cluster_id', 'start', 'end', 'feature' columns, got {header}"
        ) from e
    out: dict[str, list[tuple[int, int, str]]] = {}
    for r in rows:
        if len(r) > max(ci, sti, eni, fti):
            out.setdefault(r[ci], []).append((int(r[sti]), int(r[eni]), r[fti]))
    return out
