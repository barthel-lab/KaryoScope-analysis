"""``karyoscope-analysis detect-rearrangements`` — differential colocalization test.

Compares feature **colocalization rates** between an experiment and a control overlay BED
(optionally against a normal reference = annotated CHM13 reads), reporting feature pairs
whose proximity is differentially enriched/depleted. The headline signal for recurrent
rearrangements. See ``docs/audit/rearrangement_detection.md`` (Engine A).

Inputs are ``overlay-annotations`` outputs (one coalesced feature-segment BED per sample).
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope_analysis.core import rearrangement as core

logger = logging.getLogger(__name__)

_COLUMNS = (
    "feature_a",
    "feature_b",
    "window",
    "exp_support",
    "exp_total",
    "exp_rate",
    "ctrl_support",
    "ctrl_total",
    "ctrl_rate",
    "ref_rate",
    "log2_ratio",
    "odds_ratio",
    "p_value",
    "q_value",
    "direction",
    "reference_abnormal",
    "passes",
)


def _row(call: core.RearrangementCall) -> list[str]:
    return [
        call.pair[0],
        call.pair[1],
        str(call.window),
        str(call.exp_support),
        str(call.exp_total),
        repr(call.exp_rate),
        str(call.ctrl_support),
        str(call.ctrl_total),
        repr(call.ctrl_rate),
        repr(call.ref_rate),
        repr(call.log2_ratio),
        repr(call.odds_ratio),
        repr(call.p_value),
        repr(call.q_value),
        call.direction,
        "1" if call.reference_abnormal else "0",
        "1" if call.passes else "0",
    ]


@click.command(
    name="detect-rearrangements",
    help="Differentially test feature colocalization rates (experiment vs control).",
)
@click.option(
    "--experiment",
    "experiment_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Experiment overlay-annotations BED.",
)
@click.option(
    "--control",
    "control_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Control overlay-annotations BED.",
)
@click.option(
    "--reference",
    "reference_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Normal reference BED (annotated CHM13 reads) — the artifact floor. Optional.",
)
@click.option(
    "--length-boundary",
    "length_boundaries",
    multiple=True,
    type=int,
    metavar="BP",
    help="Read-length bucket boundary in bp (repeat for more buckets). Default: one bucket.",
)
@click.option(
    "--window",
    "windows",
    multiple=True,
    type=int,
    metavar="BP",
    help=f"Colocalization window(s) in bp (repeatable). Default: {core.DEFAULT_WINDOWS}.",
)
@click.option(
    "--min-occurrence-bp",
    default=0,
    show_default=True,
    type=int,
    help="Ignore feature intervals shorter than this when locating a feature.",
)
@click.option(
    "--min-support",
    default=3,
    show_default=True,
    type=int,
    help="Minimum supporting reads (recurrence) in the higher condition.",
)
@click.option(
    "--min-log2-ratio",
    default=1.0,
    show_default=True,
    type=float,
    help="Minimum |log2 rate ratio| (effect size) to call.",
)
@click.option(
    "--fdr",
    "fdr_alpha",
    default=0.05,
    show_default=True,
    type=float,
    help="BH-FDR significance level.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TSV of differential colocalization results.",
)
def cmd(
    experiment_path: Path,
    control_path: Path,
    reference_path: Path | None,
    length_boundaries: tuple[int, ...],
    windows: tuple[int, ...],
    min_occurrence_bp: int,
    min_support: int,
    min_log2_ratio: float,
    fdr_alpha: float,
    output: Path,
) -> None:
    """Differentially test colocalization between experiment and control."""
    # Read independence is assumed, not enforced (see the design doc): recurrence counts
    # reads as independent molecules. Surface it rather than hide it.
    logger.warning(
        "assuming reads are independent molecules (no de-duplication is performed); "
        "support counts may be inflated if the upstream pipeline leaves duplicates."
    )

    boundaries = tuple(sorted(length_boundaries))
    win = tuple(windows) if windows else core.DEFAULT_WINDOWS

    def load(path: Path) -> core.SampleColocalization:
        return core.aggregate_bed(
            str(path), boundaries=boundaries, min_occurrence_bp=min_occurrence_bp
        )

    experiment = load(experiment_path)
    control = load(control_path)
    reference = load(reference_path) if reference_path else None

    calls = core.detect_rearrangements(
        experiment,
        control,
        reference,
        windows=win,
        min_support=min_support,
        min_log2_ratio=min_log2_ratio,
        fdr_alpha=fdr_alpha,
    )

    with output.open("w", newline="") as fh:
        fh.write("\t".join(_COLUMNS) + "\n")
        for call in calls:
            fh.write("\t".join(_row(call)) + "\n")

    n_pass = sum(1 for c in calls if c.passes)
    click.echo(f"Tested {len(calls)} (pair, window) combinations; {n_pass} passed. Wrote {output}")
