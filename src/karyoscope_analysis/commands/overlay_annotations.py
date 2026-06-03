"""``karyoscope-analysis overlay-annotations`` — combine annotation BEDs into one.

Reads one annotation BED per featureset and resolves overlapping annotations to a
single feature per position, using either a built-in preset, a custom spec, or the
default basic overlay (join all featuresets with ``--separator``). Replaces the
legacy ``KaryoScope_merge_beds.py``.
"""

from __future__ import annotations

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
from karyoscope_analysis.core.io.bed import read_annotation_bed, write_annotation_bed


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
    output: Path,
) -> None:
    """Overlay annotation BEDs and resolve to one feature per position."""
    if preset and spec_path:
        raise click.UsageError("--preset and --spec are mutually exclusive.")

    bed_paths = _parse_beds(bed_specs)
    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)

    try:
        if preset:
            spec = load_builtin_preset(preset, hierarchy)
        elif spec_path:
            spec = load_spec_file(spec_path, hierarchy)
        else:
            # Default: basic overlay — join all featuresets (in --bed order) with --separator.
            spec = load_spec(
                {
                    "name": "overlay",
                    "precedence": list(bed_paths),
                    "rules": [{"emit": {"composite": "all", "sep": separator}}],
                },
                hierarchy,
            )
    except SpecError as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        beds = {fs: read_annotation_bed(path) for fs, path in bed_paths.items()}
        rows = list(core.overlay_annotations(beds, spec, hierarchy))
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    write_annotation_bed(output, rows)
    click.echo(f"Wrote {len(rows)} resolved intervals to {output}")
