"""Compare two clusterings of the same reads (label-invariant agreement + overlap).

Given two ``read_id -> cluster_id`` assignments (e.g. two ``cluster`` runs with different
parameters, or different featuresets), report how concordant they are: the Adjusted Rand Index and
Normalized Mutual Information (both permutation-invariant, since cluster ids are arbitrary between
runs), plus the cluster-to-cluster read overlap.

(The legacy ``compare_clusterings`` imported the wrong sklearn name — ``adjusted_rand_index`` rather
than ``adjusted_rand_score`` — and silently swallowed the ``ImportError``, so ARI was *never*
computed; this uses the correct name. The study-specific Pre/Post enrichment plots and fragile
sidecar discovery are dropped.)
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


@dataclass(frozen=True)
class ComparisonResult:
    """Agreement between two clusterings over their shared reads."""

    n_common: int  # reads present in both clusterings
    n_clusters_1: int
    n_clusters_2: int
    ari: float  # adjusted Rand index (1 = identical partitions, ~0 = independent)
    nmi: float  # normalized mutual information (1 = identical, 0 = independent)


def compare(labels1: Mapping[str, str], labels2: Mapping[str, str]) -> ComparisonResult:
    """Compute ARI/NMI over the reads shared by both clusterings.

    Raises ``ValueError`` if the two assignments share no reads.
    """
    common = sorted(set(labels1) & set(labels2))
    if not common:
        raise ValueError("the two clusterings share no read ids")
    l1 = [labels1[r] for r in common]
    l2 = [labels2[r] for r in common]
    return ComparisonResult(
        n_common=len(common),
        n_clusters_1=len(set(l1)),
        n_clusters_2=len(set(l2)),
        ari=float(adjusted_rand_score(l1, l2)),
        nmi=float(normalized_mutual_info_score(l1, l2)),
    )


def overlap_pairs(
    labels1: Mapping[str, str], labels2: Mapping[str, str]
) -> list[tuple[str, str, int]]:
    """Cluster-to-cluster shared-read counts ``(cluster_1, cluster_2, n_shared)``, descending.

    Long format (not a dense crosstab) so it stays compact when there are many clusters.
    """
    common = set(labels1) & set(labels2)
    pairs: Counter = Counter((labels1[r], labels2[r]) for r in common)
    return sorted(
        ((c1, c2, n) for (c1, c2), n in pairs.items()),
        key=lambda t: (-t[2], t[0], t[1]),
    )


def overlap_tsv(pairs: Sequence[tuple[str, str, int]], label1: str, label2: str) -> str:
    """Render the cluster-overlap pairs as a TSV."""
    lines = [f"cluster_{label1}\tcluster_{label2}\tn_shared"]
    lines += [f"{c1}\t{c2}\t{n}" for c1, c2, n in pairs]
    return "\n".join(lines) + "\n"


def report(result: ComparisonResult, label1: str, label2: str) -> str:
    """Human-readable concordance report."""
    if result.ari >= 0.75:
        agreement = "HIGH"
    elif result.ari >= 0.5:
        agreement = "MODERATE"
    else:
        agreement = "LOW"
    return (
        f"Clustering comparison: {label1} vs {label2}\n"
        f"  reads compared (in both): {result.n_common}\n"
        f"  clusters: {label1}={result.n_clusters_1}, {label2}={result.n_clusters_2}\n"
        f"  Adjusted Rand Index (ARI): {result.ari:.4f}\n"
        f"  Normalized Mutual Information (NMI): {result.nmi:.4f}\n"
        f"  agreement: {agreement}\n"
    )
