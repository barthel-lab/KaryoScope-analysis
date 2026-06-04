"""v2 feature vocabulary derived from a KaryoScope database hierarchy.

The KaryoScope feature taxonomy is defined authoritatively by the database's
``hierarchy.tsv`` (rows of ``feature_set``, ``child``, ``parent``). This module
loads that hierarchy and derives the feature groupings the analysis tools need
(satellites, arm, ct, telomere types) directly from it — a single source of
truth, with no hard-coded biology (decisions D4.2, D4.4, D6).

KaryoScope-analysis is **v2-only**: legacy v1 names (``*_specific`` /
``*_multigroup1``) are intentionally not recognized (decision D4.1). ``novel`` is
the one feature that legitimately does not appear in the hierarchy — it marks
k-mers absent from the database — so it is always accepted; every *other*
out-of-taxonomy feature is an error (convention C2).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

#: The only feature value allowed to be absent from the hierarchy (convention C2).
NOVEL = "novel"


class FeatureHierarchy:
    """A parsed KaryoScope feature hierarchy (a per-feature-set parent/child tree)."""

    def __init__(self, edges: Iterable[tuple[str, str, str]]) -> None:
        # feature_set -> parent -> set of children
        self._children: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        # feature_set -> set of all node names (children + parents)
        self._nodes: dict[str, set[str]] = defaultdict(set)
        for feature_set, child, parent in edges:
            self._children[feature_set][parent].add(child)
            self._nodes[feature_set].add(child)
            self._nodes[feature_set].add(parent)

    @classmethod
    def from_tsv(cls, path: str | Path) -> FeatureHierarchy:
        """Load a hierarchy from a ``hierarchy.tsv`` (header: feature_set, child, parent)."""
        path = Path(path)
        required = {"feature_set", "child", "parent"}
        edges: list[tuple[str, str, str]] = []
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(
                    f"{path}: expected a TSV with columns {sorted(required)}, "
                    f"got {reader.fieldnames!r}"
                )
            for row in reader:
                edges.append((row["feature_set"], row["child"], row["parent"]))
        if not edges:
            raise ValueError(f"{path}: hierarchy is empty")
        return cls(edges)

    # ------------------------------------------------------------------ structure
    def feature_sets(self) -> frozenset[str]:
        """The set of feature-set names (e.g. region, subtelomeric, repeat, ...)."""
        return frozenset(self._nodes)

    def features(self, feature_set: str | None = None) -> frozenset[str]:
        """All feature names — within one feature set, or across all of them."""
        if feature_set is None:
            out: set[str] = set()
            for nodes in self._nodes.values():
                out |= nodes
            return frozenset(out)
        return frozenset(self._nodes.get(feature_set, set()))

    def children(self, feature_set: str, parent: str) -> frozenset[str]:
        """Direct children of ``parent`` in ``feature_set``."""
        return frozenset(self._children.get(feature_set, {}).get(parent, set()))

    def descendants(
        self, feature_set: str, root: str, *, include_root: bool = False
    ) -> frozenset[str]:
        """All features transitively under ``root`` in ``feature_set``."""
        seen: set[str] = set()
        stack = [root]
        while stack:
            node = stack.pop()
            for child in self._children.get(feature_set, {}).get(node, set()):
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        if include_root:
            seen.add(root)
        return frozenset(seen)

    def __contains__(self, feature: object) -> bool:
        return any(feature in nodes for nodes in self._nodes.values())

    def __iter__(self) -> Iterator[str]:
        return iter(self.features())

    # ------------------------------------------------------------------ C2 validity
    def is_valid_feature(self, feature: str, feature_set: str | None = None) -> bool:
        """True if ``feature`` is in the hierarchy (or is ``novel``)."""
        return feature == NOVEL or feature in self.features(feature_set)

    def require_valid_feature(self, feature: str, feature_set: str | None = None) -> None:
        """Raise ``ValueError`` for an out-of-taxonomy feature other than ``novel`` (C2)."""
        if not self.is_valid_feature(feature, feature_set):
            where = f" in feature set {feature_set!r}" if feature_set else ""
            raise ValueError(
                f"unknown feature {feature!r}{where}: not in the database hierarchy "
                f"(only {NOVEL!r} may be absent). KaryoScope-analysis is v2-only; "
                f"legacy v1 names are not recognized."
            )

    # -------------------------------------------------- derived analysis vocabulary
    def _require_node(self, feature_set: str, node: str) -> None:
        if node not in self._nodes.get(feature_set, set()):
            raise ValueError(
                f"hierarchy is missing the {node!r} node in feature set {feature_set!r}; "
                f"cannot derive the feature vocabulary"
            )

    @property
    def satellite_features(self) -> frozenset[str]:
        """Centromeric-satellite features: the ``centromeric`` subtree, excluding ``ct``.

        ``ct`` (centromere transition) sits under ``centromeric`` in the taxonomy but
        is not a satellite sequence — it marks the gaps *between* annotated centromeric
        satellites — so it is excluded (decision D4.4).
        """
        self._require_node("region", "centromeric")
        return self.descendants("region", "centromeric", include_root=True) - {"ct"}

    @property
    def arm_features(self) -> frozenset[str]:
        """Chromosome-arm features: ``arm`` and its subtree (``p_arm``, ``q_arm``)."""
        self._require_node("region", "arm")
        return self.descendants("region", "arm", include_root=True)

    @property
    def ct_features(self) -> frozenset[str]:
        """Centromere-transition features (``ct``)."""
        self._require_node("region", "ct")
        return frozenset({"ct"})

    @property
    def canonical_telomere(self) -> frozenset[str]:
        """Canonical-telomere features (subtelomeric featureset)."""
        self._require_node("subtelomeric", "canonical_telomere")
        return frozenset({"canonical_telomere"})

    @property
    def noncanonical_telomere(self) -> frozenset[str]:
        """Noncanonical-telomere features (subtelomeric featureset)."""
        self._require_node("subtelomeric", "noncanonical_telomere")
        return frozenset({"noncanonical_telomere"})

    @property
    def its_tar1(self) -> frozenset[str]:
        """Interstitial/associated telomeric repeats: ``ITS`` and ``TAR1``."""
        for node in ("ITS", "TAR1"):
            self._require_node("subtelomeric", node)
        return frozenset({"ITS", "TAR1"})

    @property
    def interspersed_repeat_features(self) -> frozenset[str]:
        """Genome-wide interspersed/transposable repeat classes (LINE/SINE/LTR/DNA/...).

        The ``Interspersed_Repeat`` subtree of the ``repeat`` featureset. These occur in
        nearly every long read and carry no structural-rearrangement signal, so Engine B
        clustering down-weights them (otherwise reads chain together through shared repeats).
        Empty if the hierarchy has no such node.
        """
        if "Interspersed_Repeat" not in self.features("repeat"):
            return frozenset()
        return self.descendants("repeat", "Interspersed_Repeat", include_root=True)
