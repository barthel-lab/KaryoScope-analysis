"""``karyoscope-analysis cluster`` — overlap-layout-consensus clustering (Engine B).

Clusters an overlay-annotations BED subset of reads into structural haplotypes by feature-
sequence overlap (see ``docs/audit/rearrangement_detection.md``, Engine B), writing a
clusters table, a per-cluster consensus BED, and a member layout TSV. Reads are clustered
as given (subset upstream; in practice long reads only); rendering is the plotting tier.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope_analysis.core import feature_assembly as asm
from karyoscope_analysis.core.feature_align import (
    Segment,
    chromosome_aware_substitution,
    hierarchy_substitution,
    to_segments,
)
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.genome_weights import load_structural_weights
from karyoscope_analysis.core.io.bed import read_annotation_bed

logger = logging.getLogger(__name__)


def _load_reads(path: Path, *, min_length: int, min_occurrence_bp: int) -> dict[str, list[Segment]]:
    """Read an overlay BED into ``{read_id: segments}``, applying the length filters."""
    reads: dict[str, list[Segment]] = {}
    for seq_id, intervals in read_annotation_bed(path).items():
        segments = [
            (f, length) for f, length in to_segments(intervals) if length >= min_occurrence_bp
        ]
        if sum(length for _, length in segments) >= min_length and segments:
            reads[seq_id] = segments
    return reads


def _sidecar(output: Path, suffix: str) -> Path:
    name = output.name
    for ext in (".tsv.gz", ".tsv", ".gz"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return output.with_name(f"{name}.{suffix}")


@click.command(name="cluster", help="Cluster reads into structural haplotypes (OLC).")
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Overlay-annotations BED of the read subset to cluster.",
)
@click.option(
    "--hierarchy",
    "hierarchy_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database hierarchy.tsv (for the tiered substitution scorer).",
)
@click.option(
    "--min-length",
    default=0,
    show_default=True,
    type=int,
    help="Skip reads shorter than this total bp (use to keep long reads only).",
)
@click.option(
    "--min-occurrence-bp",
    default=0,
    show_default=True,
    type=int,
    help="Drop feature segments shorter than this (denoise) before aligning.",
)
@click.option(
    "--gap-factor",
    default=0.01,
    show_default=True,
    type=float,
    help="Per-bp gap cost: skipping a segment costs gap_factor x its length.",
)
@click.option("--match", default=1.0, show_default=True, type=float, help="Exact-match score.")
@click.option(
    "--partial",
    default=0.5,
    show_default=True,
    type=float,
    help="Same-group (sibling) substitution score.",
)
@click.option(
    "--mismatch", default=-1.0, show_default=True, type=float, help="Unrelated substitution score."
)
@click.option(
    "--min-overlap-bp",
    default=1000,
    show_default=True,
    type=int,
    help="Minimum overlap length (bp) for an edge.",
)
@click.option(
    "--min-identity",
    default=0.9,
    show_default=True,
    type=float,
    help="Minimum normalized identity for an edge.",
)
@click.option(
    "--min-jaccard",
    default=0.0,
    show_default=True,
    type=float,
    help="Feature-set Jaccard prefilter (0 = off).",
)
@click.option(
    "--min-distinctive-bp",
    default=0.0,
    show_default=True,
    type=float,
    help="Minimum bp of matched DISTINCTIVE features (weight >= --distinctive-weight) for an "
    "edge (0 = off). Rejects overlaps explained only by filler like a shared chromosome arm "
    "(anti-chaining); use with a weighting method.",
)
@click.option(
    "--distinctive-weight",
    default=0.15,
    show_default=True,
    type=float,
    help="Weight threshold above which a feature counts as distinctive (for --min-distinctive-bp).",
)
@click.option(
    "--block-min-bp",
    default=0.0,
    show_default=True,
    type=float,
    help="Blocking index: only align reads sharing a feature with at least this many bp in both "
    "(0 = off, all-vs-all). Lets clustering scale to whole samples by skipping the O(N^2) scan; "
    "with composite labels this compares only same-chromosome reads. Set near --min-overlap-bp.",
)
@click.option(
    "--workers",
    "-j",
    default=1,
    show_default=True,
    type=int,
    help="Align candidate read pairs across this many processes (exact; uses fork, so best on "
    "Linux/HPC). 1 = serial.",
)
@click.option(
    "--communities/--no-communities",
    default=True,
    show_default=True,
    help="Subdivide each connected component by label propagation, so a sparse bridge (a noisy "
    "multi-chromosome read, or a lone translocation) doesn't transitively merge distinct groups "
    "into one mega-cluster. --no-communities uses plain connected components.",
)
@click.option(
    "--weight-method",
    type=click.Choice(["repeat-mask", "idf", "genome-freq", "uniform"]),
    default="repeat-mask",
    show_default=True,
    help="Per-feature weighting (anti-chaining): 'repeat-mask' zeroes genome-wide "
    "interspersed-repeat features (LINE/SINE/... + nonrepeat); 'genome-freq' weights by "
    "information content from the reference genome (--genome-weights; down-weights ubiquitous "
    "features like arm, keeps rare ones like telomere); 'idf' down-weights by read frequency; "
    "'uniform' weights all equally. All operate on the structural layer of composite labels.",
)
@click.option(
    "--genome-weights",
    "genome_weights_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Weights TSV from `genome-weights` (required for --weight-method genome-freq).",
)
@click.option(
    "--weight-floor",
    default=0.1,
    show_default=True,
    type=float,
    help="Minimum weight for the most ubiquitous features (idf method).",
)
@click.option(
    "--cross-chromosome-penalty",
    default=-2.0,
    show_default=True,
    type=float,
    help="Per-bp penalty when two `chromosome:feature` labels name different specific "
    "chromosomes (soft, so translocation reads still bridge). Ignored for non-composite labels.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Clusters table TSV (consensus BED + layout TSV written as sidecars).",
)
def cmd(
    input_path: Path,
    hierarchy_path: Path,
    min_length: int,
    min_occurrence_bp: int,
    gap_factor: float,
    match: float,
    partial: float,
    mismatch: float,
    min_overlap_bp: int,
    min_identity: float,
    min_jaccard: float,
    min_distinctive_bp: float,
    distinctive_weight: float,
    block_min_bp: float,
    workers: int,
    communities: bool,
    weight_method: str,
    genome_weights_path: Path | None,
    weight_floor: float,
    cross_chromosome_penalty: float,
    output: Path,
) -> None:
    """Cluster reads into structural haplotypes and write clusters/consensus/layout."""
    reads = _load_reads(input_path, min_length=min_length, min_occurrence_bp=min_occurrence_bp)
    if not reads:
        raise click.ClickException("no reads left to cluster after filtering")

    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)
    # Structural scorer (hierarchy-tiered), wrapped to be chromosome-layer aware. The wrapper
    # is a no-op for plain (non-composite) labels, so this is safe for any overlay.
    sub_score = chromosome_aware_substitution(
        hierarchy_substitution(hierarchy, match=match, partial=partial, mismatch=mismatch),
        cross_chromosome_penalty=cross_chromosome_penalty,
    )
    if weight_method == "repeat-mask":
        masked = hierarchy.interspersed_repeat_features | {"nonrepeat"}
        weight: dict[str, float] | None = dict.fromkeys(masked, 0.0)
    elif weight_method == "idf":
        weight = asm.idf_weights(reads, floor=weight_floor)
    elif weight_method == "genome-freq":
        if genome_weights_path is None:
            raise click.UsageError("--weight-method genome-freq requires --genome-weights")
        weight = load_structural_weights(genome_weights_path)
    else:
        weight = None

    clusters, _edges = asm.assemble(
        reads,
        sub_score=sub_score,
        gap_factor=gap_factor,
        match_score=match,
        min_overlap_bp=min_overlap_bp,
        min_identity=min_identity,
        min_jaccard=min_jaccard,
        min_distinctive_bp=min_distinctive_bp,
        distinctive_weight=distinctive_weight,
        block_min_bp=block_min_bp,
        workers=workers,
        communities=communities,
        weight=weight,
    )
    layouts = [
        asm.consensus_layout(reads, c, sub_score=sub_score, gap_factor=gap_factor, weight=weight)
        for c in clusters
    ]

    clusters_tsv = output
    consensus_bed = _sidecar(output, "consensus.bed")
    layout_tsv = _sidecar(output, "layout.tsv")

    with clusters_tsv.open("w", newline="") as fh:
        fh.write("cluster_id\tsize\tseed\tconsensus_segments\twidth\torientation_conflict\n")
        for idx, (cluster, lo) in enumerate(zip(clusters, layouts, strict=True)):
            fh.write(
                f"cluster_{idx}\t{cluster.size}\t{cluster.seed}\t"
                f"{len(lo.consensus)}\t{lo.width}\t{int(cluster.orientation_conflict)}\n"
            )

    with consensus_bed.open("w", newline="") as fh:
        fh.write("cluster_id\tstart\tend\tfeature\tsupport\tcoverage\n")
        for idx, lo in enumerate(layouts):
            for p in lo.consensus:
                fh.write(
                    f"cluster_{idx}\t{p.start}\t{p.end}\t{p.feature}\t{p.support}\t{p.coverage}\n"
                )

    # per-segment layout in consensus coordinates (one row per read segment)
    with layout_tsv.open("w", newline="") as fh:
        fh.write("cluster_id\tread_id\tis_seed\treversed\tstart\tend\tfeature\n")
        for idx, lo in enumerate(layouts):
            for read in lo.placed:
                for s, e, feature in read.segments:
                    fh.write(
                        f"cluster_{idx}\t{read.read_id}\t{int(read.is_seed)}\t"
                        f"{int(read.reversed)}\t{s}\t{e}\t{feature}\n"
                    )

    n_multi = sum(1 for c in clusters if c.size > 1)
    click.echo(
        f"Clustered {len(reads)} reads into {len(clusters)} clusters "
        f"({n_multi} multi-read). Wrote {clusters_tsv}, {consensus_bed}, {layout_tsv}"
    )
