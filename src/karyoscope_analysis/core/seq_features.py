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


def max_block_bp_segments(
    segments: Sequence[tuple[int, int]], gap_tol: int = DEFAULT_BLOCK_GAP_TOL
) -> int:
    """Longest covered block over disjoint ``[start, end)`` segments, bridging gaps ``<= gap_tol``.

    The interval-native equivalent of :func:`max_block_bp` — same result, but ``O(n
    segments)`` instead of ``O(span)`` (no dense coverage array). Used by
    :func:`window_densities` so a feature's contiguity cost scales with its intervals,
    not the sequence length.
    """
    if not segments:
        return 0
    ordered = sorted(segments)
    best = 0
    run_start, run_end = ordered[0]
    for start, end in ordered[1:]:
        if start - run_end <= gap_tol:
            run_end = max(run_end, end)
        else:
            best = max(best, run_end - run_start)
            run_start, run_end = start, end
    return max(best, run_end - run_start)


def _window_sum_stats(
    segs: Sequence[tuple[int, int]], span: int, window_size: int
) -> tuple[int, int, float, int, int]:
    """Exact ``(max, min, median, first, last)`` of the 1-bp-step window coverage counts.

    ``segs`` are a feature's sorted, disjoint ``[start, end)`` intervals in ``[0, span)``;
    ``window_size <= span``. Let ``S(i)`` be the feature bp in window ``[i, i+window_size)``
    for ``i`` in ``[0, span-window_size]``. This returns the max/min/median of ``S`` plus
    ``S(0)`` and ``S(span-window_size)`` — **byte-identical** to ``np.max/min/median`` over
    the dense window-sum array — in ``O(len(segs) + window_size)`` instead of ``O(span)``.

    ``S`` is piecewise-linear (slope in ``{-1, 0, +1}``) with breakpoints only where a
    window edge crosses an interval boundary, so it is summarized run by run: max/min from
    run endpoints, and the value histogram (bounded by ``window_size``) via a difference
    array, from which the order-statistic median is read off.
    """
    w = window_size
    imax = span - w
    n = imax + 1  # number of windows

    # S(0): feature coverage in [0, w).
    s0 = 0
    for st, en in segs:
        if st >= w:
            break  # sorted -> no later interval overlaps [0, w)
        s0 += min(en, w) - max(st, 0)

    # S's slope (in {-1, 0, +1}) changes only where a window edge crosses an interval
    # boundary. Record each change as a delta at that i, then sweep — so the slope is
    # tracked incrementally (no per-run coverage probe). dS(i) = cov(i+w) - cov(i):
    # cov(i) toggles at i in {s (-1), e (+1)}; cov(i+w) at i in {s-w (+1), e-w (-1)}.
    deltas: dict[int, int] = defaultdict(int)
    for st, en in segs:
        deltas[st] -= 1
        deltas[en] += 1
        deltas[st - w] += 1
        deltas[en - w] -= 1
    positions = sorted(p for p in deltas if 0 < p <= imax)
    slope = sum(d for p, d in deltas.items() if p <= 0)  # slope on the run starting at i=0
    boundaries = [0, *positions, imax + 1]

    # Sparse value histogram as a difference map keyed by window-sum value (O(runs), not
    # O(window)): each run contributes a point mass (slope 0) or a unit-height value range.
    diff: dict[int, int] = defaultdict(int)
    cur = s0
    s_max = s_min = s_last = s0
    for k in range(len(boundaries) - 1):
        a, b = boundaries[k], boundaries[k + 1]
        count = b - a
        first_val = cur
        if slope == 0:
            if first_val > s_max:
                s_max = first_val
            elif first_val < s_min:
                s_min = first_val
            diff[first_val] += count
            diff[first_val + 1] -= count
            last_val = first_val
        else:
            last_val = cur + slope * (count - 1)
            hi, lo = (first_val, last_val) if slope < 0 else (last_val, first_val)
            if hi > s_max:
                s_max = hi
            if lo < s_min:
                s_min = lo
            diff[lo] += 1
            diff[hi + 1] -= 1
        s_last = last_val  # the final run ends at i = imax
        cur += slope * count
        if b <= imax:  # advance the slope by the event at the next boundary
            slope += deltas[b]

    # Median = mean of the two central order statistics (matches numpy for even n).
    # Walk the sparse histogram as run-length segments of constant count.
    lo_k, hi_k = (n - 1) // 2, n // 2
    keys = sorted(diff)
    v_lo = v_hi = 0
    lo_done = False
    level = cum = 0
    for j, key in enumerate(keys):
        level += diff[key]
        if level <= 0:
            continue
        width = keys[j + 1] - key  # next key always exists while level > 0
        seg_windows = level * width
        if not lo_done and lo_k < cum + seg_windows:
            v_lo = key + (lo_k - cum) // level
            lo_done = True
        if hi_k < cum + seg_windows:
            v_hi = key + (hi_k - cum) // level
            break
        cum += seg_windows
    return s_max, s_min, (v_lo + v_hi) / 2, s0, s_last


def window_densities(
    intervals: Sequence[Interval],
    *,
    window_size: int = DEFAULT_WINDOW_SIZE,
    gap_tol: int = DEFAULT_BLOCK_GAP_TOL,
) -> dict[str, DensityStats]:
    """Per-feature sliding-window density stats over the sequence's own coordinate frame.

    A window's density is the fraction of its bp covered by the feature; the window
    steps by 1 bp (fully overlapping). Sequences shorter than ``window_size`` get their
    whole-sequence coverage fraction broadcast to every density field. Computed straight
    from the intervals (``O(intervals + window_size)`` per feature) — no dense array.
    """
    if not intervals:
        return {}
    seq_start = min(start for start, _, _ in intervals)
    seq_end = max(end for _, end, _ in intervals)
    span = seq_end - seq_start
    if span <= 0:
        return {}

    by_feature: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for start, end, feature in intervals:
        by_feature[feature].append((start - seq_start, end - seq_start))  # local coords

    result: dict[str, DensityStats] = {}
    for feature, segs in by_feature.items():
        segs.sort()
        block = max_block_bp_segments(segs, gap_tol)

        if span < window_size:
            frac = sum(end - start for start, end in segs) / span
            result[feature] = DensityStats(frac, frac, frac, frac, frac, frac, frac, block)
            continue

        s_max, s_min, s_med, s_first, s_last = _window_sum_stats(segs, span, window_size)
        dfirst = s_first / window_size
        dlast = s_last / window_size
        result[feature] = DensityStats(
            dmax=s_max / window_size,
            dmin=s_min / window_size,
            dmedian=s_med / window_size,
            dfirst=dfirst,
            dlast=dlast,
            dterminal=max(dfirst, dlast),
            dterminal_min=min(dfirst, dlast),
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
