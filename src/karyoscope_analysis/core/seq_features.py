"""Per-sequence feature metrics for ``build-feature-matrix``.

Pure functions over one featureset's per-``seq_id`` annotation intervals (the C4
partition from :mod:`karyoscope_analysis.core.io.bed`), computing the columns
documented in ``docs/audit/feature_matrix_metrics.md``:

* coverage — ``bp`` (feature bp), ``frac`` (bp / total_bp), ``total_bp``;
* local sliding-window density (1-bp step) — ``dmax``/``dmin``/``dmedian``/``dfirst``/
  ``dlast``/``dterminal``/``dterminal_min``;
* longest contiguous block ``max_block_bp`` (bridging gaps ≤ ``gap_tol``);
* per-sequence ``interspersion`` (category transitions per kb).

Every magic constant (window size, block-gap tolerance, threshold factor/bounds) is
a parameter with a documented default (decision F4 — the defaults' rationale is
still TBD and tracked in ``OPEN_QUESTIONS.md``). Feature categories for interspersion
are derived from the database hierarchy (no hard-coded vocab, D4.2/D6).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import Interval

DEFAULT_WINDOW_SIZE = 1000
DEFAULT_BLOCK_GAP_TOL = 100
DEFAULT_THRESHOLD_FACTOR = 3.0
DEFAULT_THRESHOLD_MIN = 0.001
DEFAULT_THRESHOLD_MAX = 0.05


# --------------------------------------------------------------------- coverage
def total_bp(intervals: Sequence[Interval]) -> int:
    """Total bp spanned by the intervals (= the sequence length under the C4 partition)."""
    return sum(end - start for start, end, _ in intervals)


def feature_bp(intervals: Sequence[Interval]) -> dict[str, int]:
    """Per-feature total bp on the sequence."""
    out: dict[str, int] = defaultdict(int)
    for start, end, feature in intervals:
        out[feature] += end - start
    return dict(out)


def feature_fraction(intervals: Sequence[Interval]) -> dict[str, float]:
    """Per-feature coverage fraction (feature bp / total bp)."""
    total = total_bp(intervals)
    if total <= 0:
        return {}
    return {feat: bp / total for feat, bp in feature_bp(intervals).items()}


# --------------------------------------------------------------------- density
@dataclass(frozen=True)
class DensityStats:
    """Per-feature local-density statistics (densities are fractions in [0, 1])."""

    dmax: float
    dmin: float
    dmedian: float
    dfirst: float
    dlast: float
    dterminal: float
    dterminal_min: float
    max_block_bp: int


def max_block_bp(coverage: np.ndarray, gap_tol: int = DEFAULT_BLOCK_GAP_TOL) -> int:
    """Longest contiguous run of 1s in ``coverage``, merging gaps of ``<= gap_tol`` bp."""
    diffs = np.diff(np.concatenate([[0], coverage, [0]]))
    run_starts = np.where(diffs == 1)[0]
    run_ends = np.where(diffs == -1)[0]
    if len(run_starts) == 0:
        return 0
    merged_starts = [int(run_starts[0])]
    merged_ends = [int(run_ends[0])]
    for start, end in zip(run_starts[1:], run_ends[1:], strict=True):
        if start - merged_ends[-1] <= gap_tol:
            merged_ends[-1] = int(end)
        else:
            merged_starts.append(int(start))
            merged_ends.append(int(end))
    return max(end - start for start, end in zip(merged_starts, merged_ends, strict=True))


def window_densities(
    intervals: Sequence[Interval],
    *,
    window_size: int = DEFAULT_WINDOW_SIZE,
    gap_tol: int = DEFAULT_BLOCK_GAP_TOL,
) -> dict[str, DensityStats]:
    """Per-feature sliding-window density stats over the sequence's own coordinate frame.

    A window's density is the fraction of its bp covered by the feature; the window
    steps by 1 bp (fully overlapping). Sequences shorter than ``window_size`` get their
    whole-sequence coverage fraction broadcast to every density field.
    """
    if not intervals:
        return {}
    seq_start = min(start for start, _, _ in intervals)
    seq_end = max(end for _, end, _ in intervals)
    span = seq_end - seq_start
    if span <= 0:
        return {}

    by_feature: dict[str, list[Interval]] = defaultdict(list)
    for start, end, feature in intervals:
        by_feature[feature].append((start, end, feature))

    result: dict[str, DensityStats] = {}
    for feature, ivals in by_feature.items():
        coverage = np.zeros(span, dtype=np.int64)
        for start, end, _ in ivals:
            coverage[start - seq_start : end - seq_start] = 1
        block = max_block_bp(coverage, gap_tol)

        if span < window_size:
            frac = float(coverage.sum()) / span
            result[feature] = DensityStats(frac, frac, frac, frac, frac, frac, frac, block)
            continue

        cumsum = np.concatenate([[0], np.cumsum(coverage)])
        window_sums = cumsum[window_size:] - cumsum[:-window_size]
        first = float(coverage[:window_size].sum()) / window_size
        last = float(coverage[-window_size:].sum()) / window_size
        result[feature] = DensityStats(
            dmax=float(window_sums.max()) / window_size,
            dmin=float(window_sums.min()) / window_size,
            dmedian=float(np.median(window_sums)) / window_size,
            dfirst=first,
            dlast=last,
            dterminal=max(first, last),
            dterminal_min=min(first, last),
            max_block_bp=block,
        )
    return result


# --------------------------------------------------------------- interspersion
def classify_feature(feature: str, hierarchy: FeatureHierarchy) -> str:
    """Classify a feature (single- or two-layer ``a:b``) into an interspersion category.

    Categories (priority order): satellite > layer-2 telomere types > layer-1 telomere
    types > ct > arm > other. Two-layer labels arise from ``overlay-annotations``
    composite output; single-layer labels use ``layer1`` only.
    """
    parts = feature.split(":", 1)
    layer1 = parts[0]
    layer2 = parts[1] if len(parts) > 1 else ""

    if layer1 in hierarchy.satellite_features:
        return "satellite"
    if layer2 in hierarchy.canonical_telomere:
        return "canonical"
    if layer2 in hierarchy.noncanonical_telomere:
        return "noncanonical"
    if layer2 in hierarchy.its_tar1:
        return "ITS_TAR1"
    if layer1 in hierarchy.canonical_telomere:
        return "canonical"
    if layer1 in hierarchy.noncanonical_telomere:
        return "noncanonical"
    if layer1 in hierarchy.its_tar1:
        return "ITS_TAR1"
    if layer1 in hierarchy.ct_features:
        return "ct"
    if layer1 in hierarchy.arm_features:
        return "arm"
    return "other"


def interspersion(intervals: Sequence[Interval], hierarchy: FeatureHierarchy) -> dict[str, float]:
    """Per-sequence interspersion: category transitions per kb.

    Returns ``{'total', 'can_ncan', 'tel_sat', 'arm_tel'}``. ``intervals`` are assumed
    sorted by start (C4). ``total`` counts changes between adjacent categories; the
    typed counts ignore ``other``-category intervals.
    """
    zero = {"total": 0.0, "can_ncan": 0.0, "tel_sat": 0.0, "arm_tel": 0.0}
    if not intervals:
        return dict(zero)
    span_kb = (max(end for _, end, _ in intervals) - min(start for start, _, _ in intervals)) / 1000
    if span_kb <= 0:
        return dict(zero)

    categories = [classify_feature(feat, hierarchy) for _, _, feat in intervals]
    total = sum(1 for i in range(1, len(categories)) if categories[i] != categories[i - 1])

    filtered = [c for c in categories if c != "other"]
    can_ncan = tel_sat = arm_tel = 0
    for i in range(1, len(filtered)):
        if filtered[i] == filtered[i - 1]:
            continue
        pair = frozenset({filtered[i - 1], filtered[i]})
        if pair == frozenset({"canonical", "noncanonical"}):
            can_ncan += 1
        if "satellite" in pair and pair & {"canonical", "noncanonical"}:
            tel_sat += 1
        if "arm" in pair and pair & {"canonical", "noncanonical", "ITS_TAR1"}:
            arm_tel += 1

    return {
        "total": round(total / span_kb, 2),
        "can_ncan": round(can_ncan / span_kb, 2),
        "tel_sat": round(tel_sat / span_kb, 2),
        "arm_tel": round(arm_tel / span_kb, 2),
    }


# ------------------------------------------------------------- adaptive thresholds
def adaptive_thresholds(
    fractions_per_seq: Iterable[Mapping[str, float]],
    *,
    factor: float = DEFAULT_THRESHOLD_FACTOR,
    min_thresh: float = DEFAULT_THRESHOLD_MIN,
    max_thresh: float = DEFAULT_THRESHOLD_MAX,
) -> dict[str, float]:
    """Per-feature presence threshold = ``clamp(median(nonzero fracs) / factor, min, max)``.

    Consumed downstream by ``cluster-annotate`` to compute ``readpct`` columns (F5).
    Features never observed nonzero get ``min_thresh``.
    """
    all_features: set[str] = set()
    nonzero: dict[str, list[float]] = defaultdict(list)
    for row in fractions_per_seq:
        for feature, frac in row.items():
            all_features.add(feature)
            if frac > 0:
                nonzero[feature].append(frac)
    out: dict[str, float] = {}
    for feature in all_features:
        values = nonzero.get(feature)
        if values:
            out[feature] = max(min_thresh, min(max_thresh, float(np.median(values)) / factor))
        else:
            out[feature] = min_thresh  # never observed nonzero
    return out
