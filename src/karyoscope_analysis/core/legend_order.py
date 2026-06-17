"""Featureset-aware legend ordering for multi-featureset (overlay) legends.

KaryoScope's :func:`karyoscope.core.karyotype.legend_sort_key` orders one featureset's legend
(chromosomes natural-sorted, ``categorized`` pinned, hierarchy.tsv order, ``novel`` last). Our
overlay legends mix featuresets, so we sort **featureset-first** (DB natural order — the order
featuresets appear in ``hierarchy.tsv``), then by KaryoScope's within-featureset key. That groups
e.g. the telomere types (``subtelomeric``) together, the centromeric satellites (``region``)
together, and the chromosomes together, with sensible ordering inside each group — reusing the
engine's logic rather than duplicating it.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.core.karyotype import legend_sort_key

#: Sort-key type: maps a (possibly composite ``chrom:structural``) feature name to a sort tuple.
SortKey = Callable[[str], tuple]


def _structural(name: str) -> str:
    """Structural layer of a composite ``chrom:structural`` label (``chr20:gSat`` -> ``gSat``)."""
    return name.split(":", 1)[1] if ":" in name else name


def feature_sort_key(hierarchy_path: str | Path) -> SortKey:
    """Build a legend sort key over all featuresets in a database ``hierarchy.tsv``.

    The returned callable maps a feature name (structural, or composite ``chrom:structural``) to
    ``(featureset_rank, *legend_sort_key(structural, that_featureset_order))``. Features whose
    featureset can't be determined (e.g. ``novel``) sort into a trailing group, after which
    KaryoScope's key still puts ``novel`` last and ``categorized`` first.
    """
    hierarchy = parse_hierarchy(Path(hierarchy_path))
    feature_sets = list(hierarchy.feature_sets())  # DB natural order
    fs_rank = {fs: i for i, fs in enumerate(feature_sets)}
    n_fs = len(feature_sets)

    feature_fs: dict[str, str] = {}
    fs_order: dict[str, list[str]] = {}
    for fs in feature_sets:
        order = [row.child for row in hierarchy.rows_in(fs)]
        fs_order[fs] = order
        for child in order:
            feature_fs.setdefault(child, fs)  # first featureset (DB order) wins on collision

    def key(name: str) -> tuple:
        struct = _structural(name)
        fs = feature_fs.get(struct)
        rank = fs_rank.get(fs, n_fs)  # unknown featureset (e.g. novel) -> trailing group
        return (rank, *legend_sort_key(struct, fs_order.get(fs, [])))

    return key
