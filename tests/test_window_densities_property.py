"""Property test: the analytic window_densities == the old dense computation.

build-feature-matrix's density metrics were rewritten from O(span) dense arrays to
O(intervals + window) interval math. This pins that rewrite to the exact behaviour of
the previous dense implementation (kept here as the reference) over thousands of random
partitions, covering the short-sequence broadcast, span == window, sparse and
whole-span features, touching same-feature intervals, and non-zero sequence offsets.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from karyoscope_analysis.core import seq_features as sf


def _dense_reference(intervals, window_size, gap_tol):
    """The previous dense implementation of window_densities (the oracle)."""
    if not intervals:
        return {}
    seq_start = min(s for s, _, _ in intervals)
    seq_end = max(e for _, e, _ in intervals)
    span = seq_end - seq_start
    if span <= 0:
        return {}
    by_feature: dict[str, list[tuple[int, int]]] = {}
    for s, e, f in intervals:
        by_feature.setdefault(f, []).append((s, e))

    out = {}
    for feature, segs in by_feature.items():
        cov = np.zeros(span, dtype=np.int64)
        for s, e in segs:
            cov[s - seq_start : e - seq_start] = 1
        block = sf.max_block_bp(cov, gap_tol)
        if span < window_size:
            frac = float(cov.sum()) / span
            out[feature] = sf.DensityStats(frac, frac, frac, frac, frac, frac, frac, block)
            continue
        cumsum = np.concatenate([[0], np.cumsum(cov)])
        ws = cumsum[window_size:] - cumsum[:-window_size]
        first = float(cov[:window_size].sum()) / window_size
        last = float(cov[-window_size:].sum()) / window_size
        out[feature] = sf.DensityStats(
            dmax=float(ws.max()) / window_size,
            dmin=float(ws.min()) / window_size,
            dmedian=float(np.median(ws)) / window_size,
            dfirst=first,
            dlast=last,
            dterminal=max(first, last),
            dterminal_min=min(first, last),
            max_block_bp=block,
        )
    return out


def _assert_stats_equal(got: dict, exp: dict) -> None:
    assert got.keys() == exp.keys()
    for feature in exp:
        g, e = got[feature], exp[feature]
        for field in (
            "dmax",
            "dmin",
            "dmedian",
            "dfirst",
            "dlast",
            "dterminal",
            "dterminal_min",
        ):
            assert getattr(g, field) == pytest.approx(getattr(e, field), abs=1e-12), (
                feature,
                field,
                g,
                e,
            )
        assert g.max_block_bp == e.max_block_bp, (feature, "max_block_bp", g, e)


def _random_partition(rng: random.Random) -> tuple[list, int]:
    """A random C4 partition of [offset, offset+span) into feature intervals."""
    offset = rng.randrange(0, 50)
    span = rng.randrange(2, 400)
    features = ["a", "b", "c", "d"]
    n_cuts = min(rng.randrange(0, 12), span - 1)
    cuts = sorted(rng.sample(range(1, span), n_cuts)) if n_cuts else []
    bounds = [0, *cuts, span]
    intervals = [
        (offset + bounds[i], offset + bounds[i + 1], rng.choice(features))
        for i in range(len(bounds) - 1)
    ]
    return intervals, span


def test_analytic_matches_dense_reference():
    rng = random.Random(20240604)
    for _ in range(4000):
        intervals, span = _random_partition(rng)
        # window sizes spanning < span, == span, and > span (broadcast branch)
        window_size = rng.randrange(1, span + 5)
        gap_tol = rng.choice([0, 1, 5, 100])
        got = sf.window_densities(intervals, window_size=window_size, gap_tol=gap_tol)
        exp = _dense_reference(intervals, window_size, gap_tol)
        _assert_stats_equal(got, exp)
