"""``karyoscope-analysis overlay-annotations`` — combine annotation BEDs into one.

Reads one annotation BED per featureset and resolves overlapping annotations to a
single feature per position, using either a built-in preset, a custom spec, or the
default basic overlay (join all featuresets with ``--separator``). Replaces the
legacy ``KaryoScope_merge_beds.py``.
"""

from __future__ import annotations

import gzip
import multiprocessing as mp
import os
import tempfile
from collections.abc import Iterator
from pathlib import Path

import click

from karyoscope_analysis.core import overlay_annotations as core
from karyoscope_analysis.core.annotation_resolution import (
    SpecError,
    builtin_preset_names,
    load_builtin_preset,
    load_spec,
    load_spec_file,
)
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, iter_annotation_rows


def _parse_beds(bed_specs: tuple[str, ...]) -> dict[str, Path]:
    """Parse ``FEATURESET=PATH`` specs into an ordered ``{featureset: path}`` map."""
    beds: dict[str, Path] = {}
    for spec in bed_specs:
        if "=" not in spec:
            raise click.BadParameter(f"--bed must be FEATURESET=PATH, got {spec!r}")
        feature_set, path = spec.split("=", 1)
        if feature_set in beds:
            raise click.BadParameter(f"--bed featureset {feature_set!r} given more than once")
        beds[feature_set] = Path(path)
    return beds


def _write_streaming(output: Path, rows: Iterator[BedRow]) -> int:
    """Stream ``rows`` to ``output``, writing atomically (temp file + replace).

    Keeps memory flat (one row at a time) while guaranteeing the output file appears
    only if the whole overlay succeeds — a mid-stream error leaves no partial file.
    """
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


#: Rows-per-chunk target for the parallel worker pool (summed across all tracks). Chunks end
#: at a sequence boundary so no sequence is split. Matches bin-annotations' default.
DEFAULT_CHUNK_SIZE = 50_000


def _build_spec(preset, spec_path, separator, precedence, hierarchy):
    """Build the ResolutionSpec (preset / custom spec file / default composite overlay).

    Factored out so the parent (early validation + serial path) and each pool worker build
    the identical spec from the same inputs, rather than pickling a spec across processes.
    """
    if preset:
        return load_builtin_preset(preset, hierarchy)
    if spec_path is not None:
        return load_spec_file(spec_path, hierarchy)
    return load_spec(
        {
            "name": "overlay",
            "precedence": list(precedence),
            "rules": [{"emit": {"composite": "all", "sep": separator}}],
        },
        hierarchy,
    )


def _write_lines(output: Path, lines: Iterator[str]) -> int:
    """Stream pre-formatted output lines to ``output`` atomically; gzip if ``.gz``."""
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
# Overlay resolution is per-sequence independent (the k-way sweep coalesces within a seq_id
# and flushes at its boundary), so whole sequences -- each carrying every track's intervals
# -- fan out to a spawn Pool and reassemble in input order. Mirrors bin-annotations. The
# spec is rebuilt per worker (not pickled), so the output is byte-identical to serial.

_WORKER: dict = {}


def _iter_seq_groups(streams, featuresets):
    """Regroup the k concurrent row streams into per-sequence groups.

    Yields ``(seq_id, {featureset: [rows...]})`` per sequence, advancing all tracks in
    lockstep and enforcing that they present the same seq_ids in the same order -- the same
    cross-track invariant :func:`overlay_streams` checks globally.
    """
    cur = {fs: next(streams[fs], None) for fs in featuresets}
    while True:
        live = [fs for fs in featuresets if cur[fs] is not None]
        if not live:
            return
        if len(live) != len(featuresets):
            dead = [fs for fs in featuresets if cur[fs] is None]
            raise ValueError(
                f"featureset(s) {dead} ran out of sequences before the others; "
                f"featuresets must cover the same seq_ids in the same order"
            )
        seq_id = cur[featuresets[0]][0]
        for fs in featuresets[1:]:
            if cur[fs][0] != seq_id:
                raise ValueError(
                    f"featureset {fs!r} is at seq_id {cur[fs][0]!r} but expected {seq_id!r}; "
                    f"featuresets must list seq_ids in the same order"
                )
        group: dict[str, list[BedRow]] = {}
        for fs in featuresets:
            rows: list[BedRow] = []
            while cur[fs] is not None and cur[fs][0] == seq_id:
                rows.append(cur[fs])
                cur[fs] = next(streams[fs], None)
            group[fs] = rows
        yield seq_id, group


def _chunk_groups(groups, chunk_size):
    """Batch per-sequence groups into chunks of >= chunk_size total rows (whole seqs)."""
    chunk: list = []
    n = 0
    for item in groups:
        chunk.append(item)
        n += sum(len(rows) for rows in item[1].values())
        if n >= chunk_size:
            yield chunk
            chunk, n = [], 0
    if chunk:
        yield chunk


def _worker_init(hierarchy_path, featuresets, preset, spec_path, separator):
    hierarchy = FeatureHierarchy.from_tsv(Path(hierarchy_path))
    _WORKER["hierarchy"] = hierarchy
    _WORKER["featuresets"] = list(featuresets)
    _WORKER["spec"] = _build_spec(
        preset, Path(spec_path) if spec_path else None, separator, featuresets, hierarchy
    )


def _overlay_seq_chunk(chunk):
    """Overlay one chunk of whole-sequence groups -> output lines (input order preserved)."""
    spec = _WORKER["spec"]
    hierarchy = _WORKER["hierarchy"]
    featuresets = _WORKER["featuresets"]
    out: list[str] = []
    for _seq_id, group in chunk:
        streams = {fs: iter(group[fs]) for fs in featuresets}
        for sid, start, end, feature in core.overlay_streams(streams, spec, hierarchy):
            out.append(f"{sid}\t{start}\t{end}\t{feature}\n")
    return out


def _parallel_overlay_lines(bed_paths, hierarchy_path, preset, spec_path, separator, threads):
    """Yield resolved overlay lines via a spawn Pool, preserving input order."""
    featuresets = list(bed_paths)
    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=threads,
        initializer=_worker_init,
        initargs=(
            str(hierarchy_path),
            featuresets,
            preset,
            str(spec_path) if spec_path else None,
            separator,
        ),
    ) as pool:
        streams = {fs: iter_annotation_rows(path) for fs, path in bed_paths.items()}
        groups = _iter_seq_groups(streams, featuresets)
        for out_lines in pool.imap(_overlay_seq_chunk, _chunk_groups(groups, DEFAULT_CHUNK_SIZE)):
            yield from out_lines


@click.command(
    name="overlay-annotations",
    help="Combine per-featureset annotation BEDs into one resolved annotation BED.",
)
@click.option(
    "--bed",
    "bed_specs",
    multiple=True,
    required=True,
    metavar="FEATURESET=PATH",
    help="Annotation BED for one featureset (repeat once per featureset).",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database hierarchy.tsv (for spec validation + feature checks).",
)
@click.option(
    "--preset",
    default=None,
    help=f"Built-in resolution preset: {', '.join(builtin_preset_names())}.",
)
@click.option(
    "--spec",
    "spec_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Custom resolution spec (YAML). Mutually exclusive with --preset.",
)
@click.option(
    "--separator",
    default=":",
    show_default=True,
    help="Separator for the default overlay mode (when neither --preset nor --spec is given).",
)
@click.option(
    "--threads",
    "-t",
    default=1,
    show_default=True,
    type=int,
    help="Worker processes for the overlay sweep. 1 = single-threaded (default). 0 = auto "
    "(os.cpu_count()). Overlay is per-sequence independent, so the parallel output is "
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
    bed_specs: tuple[str, ...],
    hierarchy_path: Path,
    preset: str | None,
    spec_path: Path | None,
    separator: str,
    threads: int,
    output: Path,
) -> None:
    """Overlay annotation BEDs and resolve to one feature per position."""
    if preset and spec_path:
        raise click.UsageError("--preset and --spec are mutually exclusive.")

    bed_paths = _parse_beds(bed_specs)
    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)

    if threads < 0:
        raise click.BadParameter("--threads must be >= 0 (0 = auto)")
    try:
        spec = _build_spec(preset, spec_path, separator, list(bed_paths), hierarchy)
    except SpecError as exc:
        raise click.ClickException(str(exc)) from exc

    pool_size = (os.cpu_count() or 1) if threads == 0 else threads
    try:
        if pool_size > 1:
            # Parallel: regroup the tracks by sequence, overlay each sequence in the pool.
            count = _write_lines(
                output,
                _parallel_overlay_lines(
                    bed_paths, hierarchy_path, preset, spec_path, separator, pool_size
                ),
            )
        else:
            # Serial: stream every input BED concurrently — only the current interval of each.
            streams = {fs: iter_annotation_rows(path) for fs, path in bed_paths.items()}
            count = _write_streaming(output, core.overlay_streams(streams, spec, hierarchy))
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {count} resolved intervals to {output}")
