"""Per-cluster representative structures (the cluster consensus is the representative).

Engine B already computes a per-cluster **consensus** — the recurrent structural haplotype shared
by the cluster's reads — so that consensus *is* the representative. (The legacy
``select_representatives`` picked the best read per cluster by a length/feature heuristic, with a
declared-but-unused ``centroid_distance``; that read-picking is obsolete now that the consensus
exists and ``cluster-plot`` renders it directly.)

This module turns ``cluster``'s ``clusters.tsv`` + ``consensus.bed`` into a catalog of
representative haplotypes — one row per cluster with a compact, human-readable consensus
*signature* (the ordered feature path), which the raw per-segment ``consensus.bed`` doesn't give.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

Segment = tuple[int, int, str]  # (start, end, feature)


@dataclass(frozen=True)
class Representative:
    """One cluster's representative structure (its consensus)."""

    cluster_id: str
    size: int  # number of reads in the cluster
    n_segments: int  # consensus segments
    width: int  # consensus span (consensus coordinate units)
    signature: str  # ordered feature path, consecutive duplicates collapsed


def consensus_signature(segments: Sequence[Segment], sep: str = " > ") -> str:
    """Compact feature path for a consensus: features in start order, consecutive dupes collapsed.

    e.g. ``[(0,5,'canonical_telomere'), (5,9,'canonical_telomere'), (9,20,'ITS'), (20,40,'bSat')]``
    -> ``"canonical_telomere > ITS > bSat"``.
    """
    path: list[str] = []
    for _s, _e, feature in sorted(segments, key=lambda seg: seg[0]):
        if not path or path[-1] != feature:
            path.append(feature)
    return sep.join(path)


def build_catalog(
    cluster_sizes: Mapping[str, int],
    cluster_widths: Mapping[str, int],
    consensus_by_cluster: Mapping[str, Sequence[Segment]],
    *,
    min_cluster_size: int = 2,
) -> list[Representative]:
    """Build the representative catalog, keeping clusters with at least ``min_cluster_size`` reads.

    Sorted by descending size (largest, most-supported haplotypes first), then cluster id.
    """
    reps: list[Representative] = []
    for cluster_id, size in cluster_sizes.items():
        if size < min_cluster_size:
            continue
        segs = list(consensus_by_cluster.get(cluster_id, []))
        reps.append(
            Representative(
                cluster_id=cluster_id,
                size=size,
                n_segments=len(segs),
                width=cluster_widths.get(cluster_id, 0),
                signature=consensus_signature(segs),
            )
        )
    reps.sort(key=lambda r: (-r.size, r.cluster_id))
    return reps


def catalog_tsv(reps: Sequence[Representative]) -> str:
    """Render the representative catalog as a TSV (one row per cluster)."""
    lines = ["cluster_id\tsize\tn_segments\twidth\tconsensus_signature"]
    for r in reps:
        lines.append(f"{r.cluster_id}\t{r.size}\t{r.n_segments}\t{r.width}\t{r.signature}")
    return "\n".join(lines) + "\n"
