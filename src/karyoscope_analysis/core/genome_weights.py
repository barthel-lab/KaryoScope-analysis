"""Genome-frequency feature weights for Engine B (anti-chaining by information content).

The Engine B aligner needs to know which shared features are *distinctive*. Frequency in the
read set (``idf_weights``) is a poor proxy — long telomere reads are dominated by the same
ubiquitous structural features (arm), so document-frequency can't separate the drivers. The
principled signal is **how much of the reference genome each feature covers**: a feature that
tiles most of the genome (chromosome ``arm``) carries little information when two reads share
it, whereas a feature confined to a tiny fraction (``canonical_telomere``, a specific
satellite) is highly distinctive.

Given the annotated CHM13 reference (one C4 BED per featureset, ``seq_id`` = chromosome), we
tally each feature's genome bp, take its fraction ``p`` of its featureset's partition, and use
the **information content** ``-ln(p)`` as the weight, scaled to ``(0, 1]`` by the most
informative (rarest) feature across all featuresets. So the rarest feature → 1, the most
ubiquitous → ~0, on a common scale across featuresets (each featureset partitions the same
genome, so ``p`` is comparable). The weight is applied to the *structural* layer of
``chromosome:structural`` composite labels by the clusterer.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

from karyoscope_analysis.core.io.bed import BedRow

#: Columns of the genome-weights TSV written by ``genome-weights``.
WEIGHTS_HEADER = (
    "feature_set",
    "feature",
    "genome_bp",
    "genome_fraction",
    "info_content",
    "weight",
)


@dataclass(frozen=True)
class FeatureWeight:
    """One feature's genome-frequency weight, with the quantities it was derived from."""

    feature_set: str
    feature: str
    genome_bp: int
    fraction: float  # genome_bp / featureset_total_bp
    info_content: float  # -ln(fraction), in nats
    weight: float  # info_content / max info_content over all features, in (0, 1]


def tally_feature_bp(named_streams: Mapping[str, Iterator[BedRow]]) -> dict[str, dict[str, int]]:
    """Sum genome bp per feature for each featureset stream (single pass, O(intervals))."""
    out: dict[str, dict[str, int]] = {}
    for feature_set, stream in named_streams.items():
        per_feature: dict[str, int] = {}
        for _seq_id, start, end, feature in stream:
            per_feature[feature] = per_feature.get(feature, 0) + (end - start)
        out[feature_set] = per_feature
    return out


def compute_genome_weights(
    bp_by_featureset: Mapping[str, Mapping[str, int]],
) -> list[FeatureWeight]:
    """Turn per-featureset feature bp into information-content weights scaled to ``(0, 1]``.

    ``weight = -ln(p) / max(-ln(p))`` where ``p`` is the feature's fraction of its featureset's
    total bp and the max is taken across every feature in every featureset (a common scale).
    Features covering their whole partition (``p == 1``) get weight 0; the rarest feature gets
    1. Rows are returned sorted by featureset then descending weight.
    """
    info: list[tuple[str, str, int, float, float]] = []
    for feature_set, per_feature in bp_by_featureset.items():
        total = sum(per_feature.values())
        if total <= 0:
            continue
        for feature, bp in per_feature.items():
            fraction = bp / total
            ic = -math.log(fraction) if 0.0 < fraction < 1.0 else 0.0
            info.append((feature_set, feature, bp, fraction, ic))

    ic_max = max((ic for *_, ic in info), default=1.0) or 1.0
    weights = [
        FeatureWeight(fs, feat, bp, frac, ic, ic / ic_max) for fs, feat, bp, frac, ic in info
    ]
    weights.sort(key=lambda w: (w.feature_set, -w.weight, w.feature))
    return weights


def structural_weight_map(weights: list[FeatureWeight]) -> dict[str, float]:
    """Collapse weights to ``{feature: weight}`` for structural-layer lookup.

    A feature appearing in more than one featureset (e.g. ``rDNA`` in both ``region`` and
    ``acrocentric``) keeps its **largest** weight — the most-distinctive interpretation —
    matching how the clusterer looks features up by name from the structural layer.
    """
    out: dict[str, float] = {}
    for w in weights:
        if w.feature not in out or w.weight > out[w.feature]:
            out[w.feature] = w.weight
    return out


def load_structural_weights(path: str | Path) -> dict[str, float]:
    """Read a ``genome-weights`` TSV back into ``{feature: weight}`` (collapsed by name).

    Skips the header and ``#`` comments; later/larger weights win for a repeated feature
    (consistent with :func:`structural_weight_map`).
    """
    out: dict[str, float] = {}
    for raw in Path(path).read_text().splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if fields[0] == "feature_set":  # header
            continue
        feature, weight = fields[1], float(fields[5])
        if feature not in out or weight > out[feature]:
            out[feature] = weight
    return out
