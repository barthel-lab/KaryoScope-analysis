"""``karyoscope-analysis bin-annotations`` — hierarchical mode-filter of one featureset BED.

A denoise step run *before* ``overlay-annotations``: replaces each base's feature with the
locally dominant feature in a centered rolling window (a hierarchy-aware majority vote; see
``core/annotation_binning.py``), collapsing tiny fragmented segments. Input and output are
both C4 annotation BEDs of the same length, so the rest of the pipeline is unchanged.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
import tempfile
from collections.abc import Iterator
from pathlib import Path

import click

from karyoscope_analysis.core import annotation_binning as binning
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, Interval, iter_annotation_rows

logger = logging.getLogger(__name__)

#: Rows-per-chunk target for the parallel worker pool. A chunk always ends at a
#: sequence boundary, so the real size is "at least this many rows, rounded up to
#: the next whole sequence" -- binning is per-sequence, so no sequence is ever split.
DEFAULT_CHUNK_SIZE = 50_000


def _binned_rows(
    input_path: Path,
    tree: binning.BinTree,
    hierarchy: FeatureHierarchy,
    feature_set: str,
    *,
    window: int,
    step: int,
    majority_fraction: float,
    scope: str,
    novel_min_fraction: float,
) -> Iterator[BedRow]:
    """Stream the input BED one sequence at a time, yielding its binned rows.

    Sequences are contiguous (C4), so only one sequence's intervals are held at once. Every
    feature is validated against ``feature_set`` (C2) as it is read.
    """
    cur_seq: str | None = None
    buf: list[Interval] = []

    def flush() -> Iterator[BedRow]:
        if cur_seq is None:
            return
        for s, e, f in binning.bin_sequence(
            buf,
            tree,
            window=window,
            step=step,
            majority_fraction=majority_fraction,
            scope=scope,
            novel_min_fraction=novel_min_fraction,
        ):
            yield cur_seq, s, e, f

    for seq_id, start, end, feature in iter_annotation_rows(input_path):
        hierarchy.require_valid_feature(feature, feature_set)
        if seq_id != cur_seq:
            yield from flush()
            cur_seq, buf = seq_id, []
        buf.append((start, end, feature))
    yield from flush()


def _write_streaming(output: Path, rows: Iterator[BedRow]) -> int:
    """Stream ``rows`` to ``output`` atomically (temp file + replace); gzip if ``.gz``."""
    import gzip

    out_dir = str(output.parent) or "."
    fd, tmp_name = tempfile.mkstemp(dir=out_dir, prefix=f"{output.name}.", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    opener = gzip.open if output.suffix == ".gz" else open
    count = 0
    try:
        with opener(tmp, "wt", newline="") as fh:
            for seq_id, start, end, feature in rows:
                fh.write(f"{seq_id}\t{start}\t{end}\t{feature}\n")
                count += 1
        tmp.replace(output)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return count


def _write_lines(output: Path, lines: Iterator[str]) -> int:
    """Stream pre-formatted output lines to ``output`` atomically; gzip if ``.gz``."""
    import gzip

    out_dir = str(output.parent) or "."
    fd, tmp_name = tempfile.mkstemp(dir=out_dir, prefix=f"{output.name}.", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    opener = gzip.open if output.suffix == ".gz" else open
    count = 0
    try:
        with opener(tmp, "wt", newline="") as fh:
            for line in lines:
                fh.write(line)
                count += 1
        tmp.replace(output)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return count


# --- parallel path (opt-in via --threads) --------------------------------------------
# Binning is per-sequence independent, so whole sequences fan out to a process pool and
# reassemble in input order. This mirrors ``karyoscope.core.bin``'s proven Pool dispatch.
# Validation stays in the parent process (``_validated_rows``) as the single source of
# truth, so the parallel output is byte-identical to the single-threaded path.

_WORKER: dict = {}


def _worker_init(
    hierarchy_path: str,
    feature_set: str,
    window: int,
    step: int,
    majority_fraction: float,
    scope: str,
    novel_min_fraction: float,
) -> None:
    """Build this worker's binning tree once (workers receive already-validated rows)."""
    hierarchy = FeatureHierarchy.from_tsv(Path(hierarchy_path))
    _WORKER.update(
        tree=binning.BinTree.from_hierarchy(hierarchy, feature_set),
        window=window,
        step=step,
        majority_fraction=majority_fraction,
        scope=scope,
        novel_min_fraction=novel_min_fraction,
    )


def _bin_seq_chunk(chunk: list[tuple[str, list[Interval]]]) -> list[str]:
    """Bin one chunk of whole sequences into output lines (input order preserved)."""
    tree = _WORKER["tree"]
    out: list[str] = []
    for seq_id, intervals in chunk:
        for s, e, f in binning.bin_sequence(
            intervals,
            tree,
            window=_WORKER["window"],
            step=_WORKER["step"],
            majority_fraction=_WORKER["majority_fraction"],
            scope=_WORKER["scope"],
            novel_min_fraction=_WORKER["novel_min_fraction"],
        ):
            out.append(f"{seq_id}\t{s}\t{e}\t{f}\n")
    return out


def _validated_rows(
    input_path: Path, hierarchy: FeatureHierarchy, feature_set: str
) -> Iterator[BedRow]:
    """``iter_annotation_rows`` (C4) plus per-row feature validation (C2).

    The single validation path shared by the serial and parallel drivers.
    """
    for seq_id, start, end, feature in iter_annotation_rows(input_path):
        hierarchy.require_valid_feature(feature, feature_set)
        yield seq_id, start, end, feature


def _chunk_sequences(
    rows: Iterator[BedRow], chunk_size: int
) -> Iterator[list[tuple[str, list[Interval]]]]:
    """Group rows into chunks of >= ``chunk_size`` rows, never splitting a sequence."""
    chunk: list[tuple[str, list[Interval]]] = []
    chunk_rows = 0
    cur: str | None = None
    buf: list[Interval] = []
    for seq_id, start, end, feature in rows:
        if seq_id != cur:
            if cur is not None:
                chunk.append((cur, buf))
                chunk_rows += len(buf)
                if chunk_rows >= chunk_size:
                    yield chunk
                    chunk, chunk_rows = [], 0
            cur, buf = seq_id, []
        buf.append((start, end, feature))
    if cur is not None:
        chunk.append((cur, buf))
    if chunk:
        yield chunk


def _parallel_binned_lines(
    input_path: Path,
    hierarchy: FeatureHierarchy,
    hierarchy_path: Path,
    feature_set: str,
    *,
    window: int,
    step: int,
    majority_fraction: float,
    scope: str,
    novel_min_fraction: float,
    threads: int,
    chunk_size: int,
) -> Iterator[str]:
    """Yield formatted binned lines via a ``spawn`` pool, preserving input order."""
    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=threads,
        initializer=_worker_init,
        initargs=(
            str(hierarchy_path),
            feature_set,
            window,
            step,
            majority_fraction,
            scope,
            novel_min_fraction,
        ),
    ) as pool:
        chunks = _chunk_sequences(_validated_rows(input_path, hierarchy, feature_set), chunk_size)
        for out_lines in pool.imap(_bin_seq_chunk, chunks):
            yield from out_lines


@click.command(
    name="bin-annotations",
    help="Mode-filter a featureset BED (hierarchy-aware rolling window) to denoise it.",
)
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Annotation BED for one featureset (.gz handled transparently).",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database hierarchy.tsv (defines the feature tree + validates features).",
)
@click.option(
    "--feature-set",
    required=True,
    help="The featureset these annotations belong to (e.g. region, subtelomeric, chromosome).",
)
@click.option(
    "--window",
    default=binning.DEFAULT_WINDOW,
    show_default=True,
    type=int,
    help="Rolling-window size in bp (centered on each base; clipped at sequence ends).",
)
@click.option(
    "--step",
    default=binning.DEFAULT_STEP,
    show_default=True,
    type=int,
    help="Stride between window centers in bp. 1 = evaluate every base (exact, slowest). "
    ">1 strides the window for an O(intervals) speed/coarseness trade-off: output boundaries "
    "snap to the step grid and the result is no longer reverse-complement invariant. Keep "
    "step well below the smallest feature you need to localise.",
)
@click.option(
    "--majority-fraction",
    default=binning.DEFAULT_MAJORITY,
    show_default=True,
    type=float,
    help="Majority bar (tau) for descending into a child. Lower = more specific/aggressive; "
    "0 = always descend to a specific leaf (no internal/ambiguous labels).",
)
@click.option(
    "--threshold-scope",
    type=click.Choice(binning.THRESHOLD_SCOPES),
    default="node",
    show_default=True,
    help="Denominator for the majority bar: 'node' (bp at the current node; conditional "
    "majority, more specific) or 'window' (whole-window bp; conservative).",
)
@click.option(
    "--novel-min-fraction",
    default=binning.DEFAULT_NOVEL_MIN,
    show_default=True,
    type=float,
    help="Minimum window fraction for 'novel' to win a window. novel is an index property "
    "(shared across featuresets), so this absolute gate keeps the binned-novel extent "
    "featureset-independent -> overlaying yields 'novel:novel', not 'chrN:novel' mixes.",
)
@click.option(
    "--threads",
    "-t",
    default=1,
    show_default=True,
    type=int,
    help="Worker processes for binning. 1 = single-threaded (default). 0 = auto "
    "(os.cpu_count()). Binning is per-sequence independent, so the parallel output is "
    "byte-identical to single-threaded (input order preserved).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output annotation BED (.gz for gzip).",
)
def cmd(
    input_path: Path,
    hierarchy_path: Path,
    feature_set: str,
    window: int,
    step: int,
    majority_fraction: float,
    threshold_scope: str,
    novel_min_fraction: float,
    threads: int,
    output: Path,
) -> None:
    """Mode-filter a featureset BED and write the denoised C4 BED."""
    if window < 1:
        raise click.BadParameter("--window must be >= 1")
    if step < 1:
        raise click.BadParameter("--step must be >= 1")
    if threads < 0:
        raise click.BadParameter("--threads must be >= 0 (0 = auto)")
    if not 0.0 <= majority_fraction <= 1.0:
        raise click.BadParameter("--majority-fraction must be in [0, 1]")
    if not 0.0 <= novel_min_fraction <= 1.0:
        raise click.BadParameter("--novel-min-fraction must be in [0, 1]")

    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)
    if feature_set not in hierarchy.feature_sets():
        raise click.BadParameter(
            f"feature set {feature_set!r} is not in the hierarchy "
            f"(have: {', '.join(sorted(hierarchy.feature_sets()))})"
        )
    tree = binning.BinTree.from_hierarchy(hierarchy, feature_set)
    pool_size = (os.cpu_count() or 1) if threads == 0 else threads

    try:
        if pool_size > 1:
            # Parallel: validate + chunk whole sequences in the parent, bin in the pool.
            count = _write_lines(
                output,
                _parallel_binned_lines(
                    input_path,
                    hierarchy,
                    hierarchy_path,
                    feature_set,
                    window=window,
                    step=step,
                    majority_fraction=majority_fraction,
                    scope=threshold_scope,
                    novel_min_fraction=novel_min_fraction,
                    threads=pool_size,
                    chunk_size=DEFAULT_CHUNK_SIZE,
                ),
            )
        else:
            rows = _binned_rows(
                input_path,
                tree,
                hierarchy,
                feature_set,
                window=window,
                step=step,
                majority_fraction=majority_fraction,
                scope=threshold_scope,
                novel_min_fraction=novel_min_fraction,
            )
            count = _write_streaming(output, rows)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {count} binned intervals to {output}")
