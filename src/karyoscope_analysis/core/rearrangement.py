"""Recurrent-rearrangement detection — Engine A differential test (see
``docs/audit/rearrangement_detection.md``).

Builds on the alignment-free colocalization measurement
(:mod:`karyoscope_analysis.core.colocalization`). The pipeline is:

1. **Aggregate** each sample's per-read pair gaps into per-(length-bucket, pair) support
   counts (a read *supports* pair ``{A, B}`` at window ``W`` if their min bp gap ≤ ``W``).
2. **Differential test** per ``(pair, W)``: a 2x2 (support x experiment/control) table per
   length bucket, combined across buckets with the **Cochran-Mantel-Haenszel** stratified
   test (so a long-read rate is never compared against a short-read denominator), then
   **BH-FDR** across all tested ``(pair, W)``. A call must also clear a minimum recurrence,
   a minimum effect size, and the **artifact floor** — the rate seen in the normal
   reference (annotated CHM13 reads), so a call has to beat normal-genome noise.

This is a **v1 for coauthor review** (the Group A statistics). Two caveats are deliberate:

* **Read independence is assumed, not enforced.** Recurrence counts reads as independent
  molecules; if the upstream pipeline does not remove duplicates, support counts can be
  inflated. Surface this with the data, not silently. (Open item in the design doc.)
* CMH uses the asymptotic (continuity-corrected) chi-square; with the recurrence gate the
  counts are not tiny, but an exact stratified test is a possible refinement.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from scipy.stats import chi2, false_discovery_control

from karyoscope_analysis.core.colocalization import FeaturePair, iter_read_gaps

#: Default colocalization windows (bp) at which to report rates (decision: report a few;
#: 0 = adjacency). Tight windows probe novel juxtaposition, generous ones probe abundance.
DEFAULT_WINDOWS: tuple[int, ...] = (0, 1_000, 10_000, 50_000)


def length_bucket(span: int, boundaries: Sequence[int]) -> int:
    """Index of the read-length bucket for ``span`` given sorted ascending ``boundaries``.

    ``boundaries=()`` -> a single bucket (0); ``boundaries=(25_000,)`` -> bucket 0 for
    spans < 25 kb, bucket 1 otherwise.
    """
    return bisect_right(boundaries, span)


@dataclass
class SampleColocalization:
    """Aggregated per-(length-bucket, pair) colocalization gaps for one sample."""

    boundaries: tuple[int, ...]
    bucket_totals: dict[int, int]  # bucket -> number of reads
    pair_gaps: dict[int, dict[FeaturePair, list[int]]]  # bucket -> pair -> sorted min-gaps

    def buckets(self) -> list[int]:
        return sorted(self.bucket_totals)

    def total(self, bucket: int) -> int:
        return self.bucket_totals.get(bucket, 0)

    def support(self, bucket: int, pair: FeaturePair, window: int) -> int:
        """Reads in ``bucket`` whose min gap for ``pair`` is ≤ ``window``."""
        gaps = self.pair_gaps.get(bucket, {}).get(pair)
        return bisect_right(gaps, window) if gaps else 0

    def pairs(self) -> set[FeaturePair]:
        out: set[FeaturePair] = set()
        for by_pair in self.pair_gaps.values():
            out.update(by_pair)
        return out

    def pooled_rate(self, pair: FeaturePair, window: int) -> float:
        """Support / total pooled across all buckets (the simple overall rate)."""
        support = total = 0
        for bucket in self.bucket_totals:
            support += self.support(bucket, pair, window)
            total += self.total(bucket)
        return support / total if total else 0.0


def aggregate(
    read_gaps: Iterable[tuple[str, int, dict[FeaturePair, int]]],
    boundaries: Sequence[int] = (),
) -> SampleColocalization:
    """Aggregate ``(seq_id, span, {pair: gap})`` records into a :class:`SampleColocalization`."""
    bounds = tuple(sorted(boundaries))
    totals: dict[int, int] = defaultdict(int)
    pair_gaps: dict[int, dict[FeaturePair, list[int]]] = defaultdict(lambda: defaultdict(list))
    for _seq_id, span, gaps in read_gaps:
        bucket = length_bucket(span, bounds)
        totals[bucket] += 1
        for pair, gap in gaps.items():
            pair_gaps[bucket][pair].append(gap)
    for by_pair in pair_gaps.values():
        for gaps in by_pair.values():
            gaps.sort()
    return SampleColocalization(bounds, dict(totals), {b: dict(d) for b, d in pair_gaps.items()})


def aggregate_bed(
    path: str, *, boundaries: Sequence[int] = (), min_occurrence_bp: int = 0
) -> SampleColocalization:
    """Aggregate an overlay-annotation BED into a :class:`SampleColocalization`."""
    return aggregate(iter_read_gaps(path, min_occurrence_bp=min_occurrence_bp), boundaries)


def _cmh(tables: Sequence[tuple[int, int, int, int]]) -> tuple[float, float]:
    """Cochran-Mantel-Haenszel test over strata ``(a, b, c, d)`` 2x2 tables.

    ``a``/``c`` = experiment/control support, ``b``/``d`` = their non-support. Returns
    ``(p_value, mantel_haenszel_odds_ratio)``. Strata with no information contribute nothing;
    ``p=1.0`` when there is no information at all.
    """
    sum_a = sum_e = sum_v = 0.0
    or_num = or_den = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2:
            continue
        r1, r2, c1, c2 = a + b, c + d, a + c, b + d
        sum_a += a
        sum_e += r1 * c1 / n
        sum_v += r1 * r2 * c1 * c2 / (n * n * (n - 1))
        or_num += a * d / n
        or_den += b * c / n
    if sum_v <= 0:
        p_value = 1.0
    else:
        stat = (abs(sum_a - sum_e) - 0.5) ** 2 / sum_v
        p_value = float(chi2.sf(stat, 1))
    if or_den > 0:
        odds_ratio = or_num / or_den
    elif or_num > 0:
        odds_ratio = math.inf
    else:
        odds_ratio = math.nan
    return p_value, odds_ratio


def _log2_ratio(exp_support: int, exp_total: int, ctrl_support: int, ctrl_total: int) -> float:
    """log2 rate ratio with a Haldane-Anscombe (+0.5 / +1) correction to avoid 0/inf."""
    exp_rate = (exp_support + 0.5) / (exp_total + 1)
    ctrl_rate = (ctrl_support + 0.5) / (ctrl_total + 1)
    return math.log2(exp_rate / ctrl_rate)


@dataclass(frozen=True)
class RearrangementCall:
    """One ``(pair, window)`` differential-colocalization result."""

    pair: FeaturePair
    window: int
    exp_support: int
    exp_total: int
    exp_rate: float
    ctrl_support: int
    ctrl_total: int
    ctrl_rate: float
    ref_rate: float  # normal (CHM13) rate = artifact floor
    log2_ratio: float
    odds_ratio: float
    p_value: float
    q_value: float
    direction: str  # "enriched" | "depleted"
    reference_abnormal: bool  # pair rarely colocalizes in the reference at this window
    passes: bool


def detect_rearrangements(
    experiment: SampleColocalization,
    control: SampleColocalization,
    reference: SampleColocalization | None = None,
    *,
    windows: Sequence[int] = DEFAULT_WINDOWS,
    min_support: int = 3,
    min_log2_ratio: float = 1.0,
    fdr_alpha: float = 0.05,
    floor_eps: float = 0.0,
    candidate_pairs: set[FeaturePair] | None = None,
) -> list[RearrangementCall]:
    """Differentially test colocalization rates between experiment and control.

    Tests every pair observed in the experiment (or ``candidate_pairs`` if given) at each
    window. A pair passes when it is FDR-significant, recurrent, has a large enough effect,
    and the higher condition clears the reference (artifact) floor. ``reference`` is the
    normal baseline (annotated CHM13 reads); without it the floor is 0.

    Returns all tested ``(pair, window)`` calls (so nothing is silently dropped), sorted with
    passing calls first by q-value.
    """
    pairs = candidate_pairs if candidate_pairs is not None else experiment.pairs()
    all_buckets = sorted(set(experiment.buckets()) | set(control.buckets()))

    records: list[dict] = []
    for pair in pairs:
        for window in windows:
            tables: list[tuple[int, int, int, int]] = []
            es = et = cs = ct = 0
            for bucket in all_buckets:
                a = experiment.support(bucket, pair, window)
                a_tot = experiment.total(bucket)
                c = control.support(bucket, pair, window)
                c_tot = control.total(bucket)
                tables.append((a, a_tot - a, c, c_tot - c))
                es += a
                et += a_tot
                cs += c
                ct += c_tot
            if et == 0 and ct == 0:
                continue
            p_value, odds_ratio = _cmh(tables)
            records.append(
                {
                    "pair": pair,
                    "window": window,
                    "es": es,
                    "et": et,
                    "cs": cs,
                    "ct": ct,
                    "ref_rate": reference.pooled_rate(pair, window) if reference else 0.0,
                    "p": p_value,
                    "or": odds_ratio,
                }
            )

    q_values = (
        list(false_discovery_control([r["p"] for r in records], method="bh")) if records else []
    )

    calls: list[RearrangementCall] = []
    for r, q in zip(records, q_values, strict=True):
        exp_rate = r["es"] / r["et"] if r["et"] else 0.0
        ctrl_rate = r["cs"] / r["ct"] if r["ct"] else 0.0
        log2_ratio = _log2_ratio(r["es"], r["et"], r["cs"], r["ct"])
        direction = "enriched" if exp_rate >= ctrl_rate else "depleted"
        higher_support = max(r["es"], r["cs"])
        higher_rate = max(exp_rate, ctrl_rate)
        passes = (
            q <= fdr_alpha
            and higher_support >= min_support
            and abs(log2_ratio) >= min_log2_ratio
            and higher_rate > r["ref_rate"] + floor_eps
        )
        calls.append(
            RearrangementCall(
                pair=r["pair"],
                window=r["window"],
                exp_support=r["es"],
                exp_total=r["et"],
                exp_rate=exp_rate,
                ctrl_support=r["cs"],
                ctrl_total=r["ct"],
                ctrl_rate=ctrl_rate,
                ref_rate=r["ref_rate"],
                log2_ratio=log2_ratio,
                odds_ratio=r["or"],
                p_value=r["p"],
                q_value=float(q),
                direction=direction,
                reference_abnormal=r["ref_rate"] <= floor_eps,
                passes=passes,
            )
        )

    calls.sort(key=lambda c: (not c.passes, c.q_value, c.pair, c.window))
    return calls
