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
from karyoscope_analysis.core.feature_align import Segment, hierarchy_substitution, to_segments
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
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
    output: Path,
) -> None:
    """Cluster reads into structural haplotypes and write clusters/consensus/layout."""
    reads = _load_reads(input_path, min_length=min_length, min_occurrence_bp=min_occurrence_bp)
    if not reads:
        raise click.ClickException("no reads left to cluster after filtering")

    hierarchy = FeatureHierarchy.from_tsv(hierarchy_path)
    sub_score = hierarchy_substitution(hierarchy, match=match, partial=partial, mismatch=mismatch)

    clusters, _edges = asm.assemble(
        reads,
        sub_score=sub_score,
        gap_factor=gap_factor,
        match_score=match,
        min_overlap_bp=min_overlap_bp,
        min_identity=min_identity,
        min_jaccard=min_jaccard,
    )
    consensuses = [
        asm.cluster_consensus(reads, c, sub_score=sub_score, gap_factor=gap_factor)
        for c in clusters
    ]

    clusters_tsv = output
    consensus_bed = _sidecar(output, "consensus.bed")
    layout_tsv = _sidecar(output, "layout.tsv")

    with clusters_tsv.open("w", newline="") as fh:
        fh.write("cluster_id\tsize\tseed\tconsensus_segments\torientation_conflict\n")
        for idx, (cluster, cons) in enumerate(zip(clusters, consensuses, strict=True)):
            fh.write(
                f"cluster_{idx}\t{cluster.size}\t{cluster.seed}\t"
                f"{len(cons.positions)}\t{int(cluster.orientation_conflict)}\n"
            )

    with consensus_bed.open("w", newline="") as fh:
        fh.write("cluster_id\tstart\tend\tfeature\tsupport\tcoverage\n")
        for idx, cons in enumerate(consensuses):
            pos = 0
            for p in cons.positions:
                fh.write(
                    f"cluster_{idx}\t{pos}\t{pos + p.length}\t{p.feature}\t{p.support}\t{p.coverage}\n"
                )
                pos += p.length

    with layout_tsv.open("w", newline="") as fh:
        fh.write("cluster_id\tread_id\tis_seed\treversed\n")
        for idx, cluster in enumerate(clusters):
            for member in cluster.members:
                fh.write(
                    f"cluster_{idx}\t{member}\t{int(member == cluster.seed)}\t"
                    f"{int(cluster.reversed_relative_to_seed[member])}\n"
                )

    n_multi = sum(1 for c in clusters if c.size > 1)
    click.echo(
        f"Clustered {len(reads)} reads into {len(clusters)} clusters "
        f"({n_multi} multi-read). Wrote {clusters_tsv}, {consensus_bed}, {layout_tsv}"
    )
