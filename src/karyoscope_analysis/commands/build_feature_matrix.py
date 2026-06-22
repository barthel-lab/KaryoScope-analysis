"""``karyoscope-analysis build-feature-matrix`` — wide per-sequence feature matrix.

Reads one annotation BED per featureset and emits a wide per-``seq_id`` TSV of
coverage / density / contiguity / interspersion metrics (the input to
``cluster``/``cluster-annotate``), plus an adaptive-threshold sidecar. Replaces the
matrix-building part of the legacy ``KaryoScope_sequence_annotate.py`` (alignment-QC
columns move to ``cluster-diagnostics``; decision F6).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import click

from karyoscope_analysis.core import build_feature_matrix as core
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import iter_annotation_rows


def _parse_beds(bed_specs: tuple[str, ...]) -> dict[str, Path]:
    beds: dict[str, Path] = {}
    for spec in bed_specs:
        if "=" not in spec:
            raise click.BadParameter(f"--bed must be FEATURESET=PATH, got {spec!r}")
        feature_set, path = spec.split("=", 1)
        if feature_set in beds:
            raise click.BadParameter(f"--bed featureset {feature_set!r} given more than once")
        beds[feature_set] = Path(path)
    return beds


def _fmt(value: float) -> str:
    """Compact, deterministic cell formatting: ints stay int-like, floats round-trip."""
    return str(int(value)) if isinstance(value, int) else repr(value)


def _write_matrix(path: Path, matrix: core.FeatureMatrix) -> None:
    header = ["seq_id", *matrix.columns]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for seq_id in matrix.seq_ids():
            row = matrix.rows[seq_id]
            cells = [seq_id, *(_fmt(row.get(col, 0)) for col in matrix.columns)]
            fh.write("\t".join(cells) + "\n")


def _write_thresholds(path: Path, thresholds: list[tuple[str, str, float]]) -> None:
    with path.open("w", newline="") as fh:
        fh.write("featureset\tfeature\tthreshold\n")
        for feature_set, feature, threshold in thresholds:
            fh.write(f"{feature_set}\t{feature}\t{threshold!r}\n")


def _default_thresholds_path(output: Path) -> Path:
    """``foo.tsv``/``foo.tsv.gz`` -> ``foo.adaptive_thresholds.tsv`` (always plain TSV)."""
    name = output.name
    for suffix in (".tsv.gz", ".tsv", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return output.with_name(f"{name}.adaptive_thresholds.tsv")


@click.command(
    name="build-feature-matrix",
    help="Build the wide per-sequence feature matrix (+ adaptive-threshold sidecar).",
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
    help="Database hierarchy.tsv (feature validation + interspersion categories).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output wide matrix TSV (.gz for gzip).",
)
@click.option(
    "--thresholds-output",
    "thresholds_path",
    default=None,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Adaptive-thresholds sidecar TSV (default: <output>.adaptive_thresholds.tsv).",
)
@click.option(
    "--interspersion-featureset",
    default=None,
    help="Featureset to compute interspersion over (typically an overlay composite). "
    "Omitted by default.",
)
@click.option("--window-size", default=core.sf.DEFAULT_WINDOW_SIZE, show_default=True, type=int)
@click.option("--block-gap-tol", default=core.sf.DEFAULT_BLOCK_GAP_TOL, show_default=True, type=int)
@click.option(
    "--threshold-factor", default=core.sf.DEFAULT_THRESHOLD_FACTOR, show_default=True, type=float
)
@click.option(
    "--threshold-min", default=core.sf.DEFAULT_THRESHOLD_MIN, show_default=True, type=float
)
@click.option(
    "--threshold-max", default=core.sf.DEFAULT_THRESHOLD_MAX, show_default=True, type=float
)
def cmd(
    bed_specs: tuple[str, ...],
    hierarchy_path: Path,
    output: Path,
    thresholds_path: Path | None,
    interspersion_featureset: str | None,
    window_size: int,
    block_gap_tol: int,
    threshold_factor: float,
    threshold_min: float,
    threshold_max: float,
) -> None:
    """Build the wide per-sequence feature matrix."""
    bed_paths = _parse_beds(bed_specs)
    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)
    # Stream every input BED concurrently (lockstep by seq_id) — only one sequence's
    # intervals per featureset is held at a time, never the whole files.
    streams = {fs: iter_annotation_rows(p) for fs, p in bed_paths.items()}

    try:
        matrix = core.build_feature_matrix_streaming(
            streams,
            hierarchy,
            window_size=window_size,
            gap_tol=block_gap_tol,
            threshold_factor=threshold_factor,
            threshold_min=threshold_min,
            threshold_max=threshold_max,
            interspersion_featureset=interspersion_featureset,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    _write_matrix(output, matrix)
    sidecar = thresholds_path or _default_thresholds_path(output)
    _write_thresholds(sidecar, matrix.thresholds)
    click.echo(
        f"Wrote {len(matrix.seq_ids())} sequences x {len(matrix.columns)} columns to {output}\n"
        f"Wrote {len(matrix.thresholds)} adaptive thresholds to {sidecar}"
    )
