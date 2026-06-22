"""``karyoscope-analysis genome-weights`` — per-feature weights from the reference genome.

Tallies how much of the annotated CHM13 reference each feature covers (one C4 BED per
featureset) and writes information-content weights in ``(0, 1]`` (ubiquitous features → ~0,
rare/distinctive → 1; see ``core/genome_weights.py``). Engine B (`cluster --weight-method
genome-freq`) applies these to the structural layer so overlaps rest on distinctive content.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import genome_weights as gw
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import iter_annotation_rows


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
    name="genome-weights",
    help="Compute per-feature information-content weights from the reference genome.",
)
@click.option(
    "--bed",
    "bed_specs",
    multiple=True,
    required=True,
    metavar="FEATURESET=PATH",
    help="Reference annotation BED for one featureset (repeat once per featureset).",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database hierarchy.tsv (validates features; C2).",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output weights TSV (feature_set, feature, genome_bp, genome_fraction, info_content, weight).",
)
def cmd(bed_specs: tuple[str, ...], hierarchy_path: Path, output: Path) -> None:
    """Tally genome bp per feature and write information-content weights."""
    bed_paths = _parse_beds(bed_specs)
    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)

    streams = {fs: iter_annotation_rows(path) for fs, path in bed_paths.items()}
    try:
        bp_by_featureset = gw.tally_feature_bp(streams)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    for feature_set, per_feature in bp_by_featureset.items():
        for feature in per_feature:
            try:
                hierarchy.require_valid_feature(feature, feature_set)
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc

    weights = gw.compute_genome_weights(bp_by_featureset)
    with output.open("w", newline="") as fh:
        fh.write("\t".join(gw.WEIGHTS_HEADER) + "\n")
        for w in weights:
            fh.write(
                f"{w.feature_set}\t{w.feature}\t{w.genome_bp}\t"
                f"{w.fraction:.3e}\t{w.info_content:.4f}\t{w.weight:.4f}\n"
            )
    click.echo(
        f"Wrote {len(weights)} feature weights ({len(bp_by_featureset)} featuresets) to {output}"
    )
