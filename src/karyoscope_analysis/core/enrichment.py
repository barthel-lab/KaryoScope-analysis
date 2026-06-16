"""Cluster enrichment across samples/groups (descriptive, compositional).

After pooling reads from several samples and clustering them once (Engine B), each cluster is a
recurrent structural haplotype that may be over-represented in one sample/group. This module
answers, per cluster, *which group concentrates on it* — the headline being structures that
define a line (e.g. ALT-associated architectures private to U2OS).

The model (see the design discussion):

- The **count unit is the read**, assigned to ``(cluster, group)``; ``group`` is a sample or a
  sample-group from a read-list.
- Comparison is **compositional**: each group's cluster sizes are normalized to that group's read
  total (``fraction``), which absorbs sequencing-depth differences. (bp-depth and per-chromosome-end
  copy-number normalization are planned refinements.)
- Enrichment of cluster *k* in group *g* is ``log2(fraction_{k,g} / pooled_fraction_k)`` — how
  over-/under-represented the cluster is in that group vs. the pool.
- **No biological replication** (typically one sample per line), so reads are *not* independent
  replicates of prevalence: this is an **exploratory, effect-size-first** report, not a hypothesis
  test. The ``enriched`` flag is an effect-size call, deliberately not a per-read p-value (which
  would mostly measure depth). A formal group contrast is a future addition for replicated designs.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class ClusterEnrichment:
    """Per-cluster enrichment across groups (descriptive)."""

    cluster_id: str
    n_total: int
    counts: dict[str, int]  # group -> reads in this cluster
    fractions: dict[str, float]  # group -> reads_in_cluster / group's total reads
    effects: dict[str, float]  # group -> log2(fraction / pooled cluster fraction)
    top_group: str  # the group most over-represented in this cluster (max effect)
    private: bool  # reads come from exactly one group
    enriched: bool  # exploratory call: top effect >= min_log2_effect


def compute_enrichment(
    read_to_cluster: Mapping[str, str],
    read_to_group: Mapping[str, str],
    *,
    min_log2_effect: float = 1.0,
    min_cluster_size: int = 2,
) -> tuple[list[ClusterEnrichment], list[str], dict[str, int]]:
    """Compute per-cluster, per-group enrichment.

    Args:
        read_to_cluster: ``{read_id: cluster_id}`` (e.g. from ``cluster``'s ``layout.tsv``).
        read_to_group: ``{read_id: group}`` (sample or sample-group, from a read-list).
        min_log2_effect: a cluster is flagged ``enriched`` when some group's effect (log2 of its
            fraction over the pooled cluster fraction) is at least this. Default 1.0 (= 2x).
        min_cluster_size: a cluster must have at least this many reads to be flagged ``enriched``.
            A single-read cluster is trivially "private" with a large fold-change but carries no
            compositional evidence, so it is never called enriched. Default 2.

    Returns:
        ``(results, groups, group_totals)`` — the per-cluster results (sorted by descending top
        effect, then size), the sorted group list, and each group's total clustered-read count.
        Only reads present in both maps are counted.
    """
    groups = sorted({read_to_group[r] for r in read_to_cluster if r in read_to_group})
    counts: dict[str, Counter] = defaultdict(Counter)
    group_totals: Counter = Counter()
    total = 0
    for read, cluster in read_to_cluster.items():
        group = read_to_group.get(read)
        if group is None:
            continue
        counts[cluster][group] += 1
        group_totals[group] += 1
        total += 1

    results: list[ClusterEnrichment] = []
    for cluster, group_counts in counts.items():
        n_total = sum(group_counts.values())
        pooled_fraction = n_total / total  # expected per-group fraction under independence
        fractions: dict[str, float] = {}
        effects: dict[str, float] = {}
        for g in groups:
            gt = group_totals[g]
            frac = group_counts.get(g, 0) / gt if gt else 0.0
            fractions[g] = frac
            # log2 fold-enrichment vs the pool; -inf floored to a large negative for absent groups.
            effects[g] = math.log2(frac / pooled_fraction) if frac > 0 else float("-inf")
        present = [g for g in groups if group_counts.get(g, 0) > 0]
        top_group = max(groups, key=lambda g: effects[g])
        results.append(
            ClusterEnrichment(
                cluster_id=cluster,
                n_total=n_total,
                counts={g: group_counts.get(g, 0) for g in groups},
                fractions=fractions,
                effects=effects,
                top_group=top_group,
                private=len(present) == 1,
                enriched=n_total >= min_cluster_size and max(effects.values()) >= min_log2_effect,
            )
        )

    results.sort(key=lambda r: (-r.effects[r.top_group], -r.n_total, r.cluster_id))
    return results, groups, dict(group_totals)


def _fmt(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.4g}"


def enrichment_tsv(results: Sequence[ClusterEnrichment], groups: Sequence[str]) -> str:
    """Render enrichment results as a tidy TSV (one row per cluster)."""
    header = (
        ["cluster_id", "n_total"]
        + [f"n_{g}" for g in groups]
        + [f"frac_{g}" for g in groups]
        + [f"log2fc_{g}" for g in groups]
        + ["top_group", "private", "enriched"]
    )
    lines = ["\t".join(header)]
    for r in results:
        row = (
            [r.cluster_id, str(r.n_total)]
            + [str(r.counts[g]) for g in groups]
            + [_fmt(r.fractions[g]) for g in groups]
            + [_fmt(r.effects[g]) for g in groups]
            + [r.top_group, str(int(r.private)), str(int(r.enriched))]
        )
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def summarize(results: Sequence[ClusterEnrichment], groups: Sequence[str]) -> str:
    """One-line-per-group summary of enriched clusters (for stdout).

    Counts only *callable* enriched clusters (the size gate already excluded singletons), and of
    those, how many are fully private to the group — the strongest, most interpretable hits.
    """
    out = []
    for g in groups:
        enriched = [r for r in results if r.enriched and r.top_group == g]
        n_private = sum(1 for r in enriched if r.private)
        out.append(f"  {g}: {len(enriched)} enriched cluster(s) ({n_private} fully private)")
    return "\n".join(out)
