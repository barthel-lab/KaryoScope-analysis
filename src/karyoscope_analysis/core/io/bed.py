"""Annotation-BED reading/writing with the C4 input invariant.

An *annotation* is a 4-column BED: ``seq_id``, ``start``, ``end``, ``feature``
(tab-separated), where ``seq_id`` (column 0) is a generic sequence identifier — a
read, contig, or chromosome (convention C1) — not a genome chromosome, and
coordinates are in that sequence's own ``[start, end)`` frame.

The C4 invariant (validated on read; a violation is an error, convention C2):

* Rows are grouped by ``seq_id`` (all rows for a sequence are contiguous).
* Within each sequence the intervals **partition** the sequence: sorted, with no
  gaps and no overlaps (each interval's ``start`` equals the previous ``end``).

Downstream code relies on this and never re-sorts. ``.gz`` inputs (plain gzip or
bgzip/BGZF, which is gzip-compatible to read) are handled transparently.

Feature-value validity (C2: only ``novel`` may be out-of-taxonomy) is enforced
separately by consumers via :mod:`karyoscope_analysis.core.feature_vocab`, so this
reader stays purely structural.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path

#: A parsed annotation row.
BedRow = tuple[str, int, int, str]
#: A valued interval within one sequence (start, end, feature).
Interval = tuple[int, int, str]


@contextmanager
def _open_text(path: Path, mode: str = "r") -> Iterator:
    """Open ``path`` as text, transparently (de)compressing ``.gz`` files."""
    if path.suffix == ".gz":
        fh = gzip.open(path, mode + "t", newline="")  # noqa: SIM115 (closed in finally below)
    else:
        fh = path.open(mode, newline="")
    try:
        yield fh
    finally:
        fh.close()


def iter_bed_rows(path: str | Path) -> Iterator[BedRow]:
    """Yield ``(seq_id, start, end, feature)`` rows, validating each row's shape.

    Raises ``ValueError`` (with the file path and 1-based line number) for a row
    with fewer than 4 tab-separated fields, non-integer coordinates, or a
    non-positive interval (``end <= start``). Columns beyond the 4th are ignored.
    Blank lines and ``#`` comment lines are skipped.
    """
    path = Path(path)
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{lineno}: expected >=4 tab-separated fields "
                    f"(seq_id, start, end, feature), got {len(fields)}"
                )
            seq_id, start_s, end_s, feature = fields[0], fields[1], fields[2], fields[3]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{lineno}: non-integer coordinates {start_s!r}, {end_s!r}"
                ) from exc
            if end <= start:
                raise ValueError(
                    f"{path}:{lineno}: non-positive interval [{start}, {end}) for {seq_id!r}"
                )
            yield seq_id, start, end, feature


def iter_annotation_rows(path: str | Path, *, validate: bool = True) -> Iterator[BedRow]:
    """Stream ``(seq_id, start, end, feature)`` rows, enforcing the C4 invariant on the fly.

    The streaming counterpart of :func:`read_annotation_bed`: it never holds more than
    one row in memory, so several BEDs can be overlaid in a single pass without loading
    any of them (see :func:`karyoscope_analysis.core.overlay_annotations.overlay_streams`).

    Args:
        path: BED path (``.gz`` handled transparently).
        validate: enforce C4 (contiguous ``seq_id`` grouping + a gapless,
            non-overlapping partition per sequence).

    Raises:
        ValueError: on a malformed row (see :func:`iter_bed_rows`) or, when
            ``validate`` is set, a C4 violation (non-contiguous ``seq_id``, a gap,
            or an overlap).
    """
    path = Path(path)
    seen: set[str] = set()
    prev_seq: str | None = None
    prev_end: int | None = None

    for seq_id, start, end, feature in iter_bed_rows(path):
        if seq_id != prev_seq:
            if validate and seq_id in seen:
                raise ValueError(
                    f"{path}: rows for seq_id {seq_id!r} are not contiguous; the file "
                    f"must be grouped/sorted by seq_id (C4 invariant)"
                )
            seen.add(seq_id)
            prev_seq = seq_id
            prev_end = None
        elif validate and prev_end is not None:
            if start < prev_end:
                raise ValueError(
                    f"{path}: overlap in seq_id {seq_id!r}: interval starts at {start} "
                    f"but previous ended at {prev_end} (C4 requires a partition)"
                )
            if start > prev_end:
                raise ValueError(
                    f"{path}: gap in seq_id {seq_id!r}: interval starts at {start} but "
                    f"previous ended at {prev_end} (C4 requires a gapless tiling)"
                )
        yield seq_id, start, end, feature
        prev_end = end


def read_annotation_bed(path: str | Path, *, validate: bool = True) -> dict[str, list[Interval]]:
    """Read a 4-column annotation BED grouped by ``seq_id``, preserving file order.

    Eagerly materializes the whole file as ``{seq_id: [(start, end, feature), ...]}``.
    For overlaying several BEDs without loading them, prefer the streaming
    :func:`iter_annotation_rows`.

    Args:
        path: BED path (``.gz`` handled transparently).
        validate: enforce the C4 invariant (contiguous ``seq_id`` grouping + a
            gapless, non-overlapping partition per sequence). Leave on unless a
            caller explicitly needs non-partitioned input.

    Returns:
        ``{seq_id: [(start, end, feature), ...]}`` in the order rows appear.

    Raises:
        ValueError: on a malformed row (see :func:`iter_bed_rows`) or, when
            ``validate`` is set, a C4 violation (non-contiguous ``seq_id``, a gap,
            or an overlap).
    """
    groups: dict[str, list[Interval]] = {}
    for seq_id, start, end, feature in iter_annotation_rows(path, validate=validate):
        groups.setdefault(seq_id, []).append((start, end, feature))
    return groups


def write_annotation_bed(path: str | Path, rows: Iterable[BedRow]) -> None:
    """Write 4-column annotation rows to ``path`` (gzip if the suffix is ``.gz``).

    NOTE (C3): KaryoScope-derived pipelines should ultimately emit **bgzip** output.
    Plain ``.gz`` here is gzip (not block-gzip), which is fine for per-``seq_id``
    annotation BEDs (column 0 is a sequence id, not a genome coordinate, so tabix
    indexing does not apply). True BGZF output — via the ``bgzip`` tool or a
    library — will be wired in for genome-coordinate outputs where indexing matters.
    """
    path = Path(path)
    with _open_text(path, "w") as fh:
        for seq_id, start, end, feature in rows:
            fh.write(f"{seq_id}\t{start}\t{end}\t{feature}\n")
