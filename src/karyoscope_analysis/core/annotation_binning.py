"""Hierarchical mode-filter ("binning") of a per-featureset annotation BED.

A denoise step run **before** ``overlay-annotations``: it replaces each base's feature
with the locally dominant feature in a centered rolling window, collapsing the tiny
fragmented segments that otherwise fragment the feature sequences (hurting Engine B
clustering and the read plots). The sequence length is unchanged — the output is still a
C4-valid gapless partition — so everything downstream is unaffected.

The vote is **hierarchy-aware** rather than a flat mode (decided with the maintainer). The
per-feature base counts in the window propagate up the database tree; the reported feature
is found by descending from the root, stepping into the dominant child while its subtree
holds a majority, and stopping at the deepest node that still does. This means related
siblings (e.g. ``aSat``/``bSat`` under ``centromeric``) reinforce each other instead of
splitting their vote and losing to an unrelated minority — and when the subtypes themselves
split, the honest call is their *ancestor* (``centromeric``), which the descent returns.

Two knobs:

* ``majority_fraction`` (τ, default 0.5) — the majority bar for descending.
* ``threshold_scope`` — the denominator τ is applied to:
    - ``"node"`` (default): τ · (bp at the current node). The *conditional* majority —
      "within what we've already committed to, does one sub-branch dominate?" More
      specific; can descend to a leaf holding < 50 % of the whole window when the
      top-level split is near-even.
    - ``"window"``: τ · (whole-window bp). Conservative; never reports a label covering
      less than τ of the window, so it climbs to internal nodes more readily.

If the descent can't even leave the root (no top-level group has a majority — e.g. a clean
50/50 boundary between two unrelated features), it falls back to flat plurality (most bp,
ties broken toward the deeper/more-specific label). That keeps boundaries sharp instead of
smearing the generic root label across them.
"""

from __future__ import annotations

import bisect
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from karyoscope_analysis.core.feature_vocab import NOVEL, FeatureHierarchy
from karyoscope_analysis.core.io.bed import Interval

#: Default rolling-window size (bp). Odd so the window is symmetric about each base.
DEFAULT_WINDOW = 101
#: Default majority fraction (τ) for descending into a child. 0 = always descend to a leaf.
DEFAULT_MAJORITY = 0.5
#: Default minimum window fraction for ``novel`` to win a window (see :func:`descend`).
DEFAULT_NOVEL_MIN = 0.5
#: Valid values for ``threshold_scope``.
THRESHOLD_SCOPES = ("node", "window")

_INF = float("inf")
_EPS = 1e-9


@dataclass(frozen=True)
class BinTree:
    """A single featureset's feature tree, indexed for the binning descent.

    ``children`` maps a node to its (sorted) children; ``parent`` is the inverse;
    ``depth`` is the root-to-node distance (root = 0). Labels outside the tree (only
    ``novel`` survives the C2 check upstream) are handled by the descent as top-level
    leaves, so they are not present here.
    """

    parent: Mapping[str, str]
    children: Mapping[str, tuple[str, ...]]
    depth: Mapping[str, int]
    root: str

    @classmethod
    def from_hierarchy(cls, hierarchy: FeatureHierarchy, feature_set: str) -> BinTree:
        """Build the tree for ``feature_set`` from a parsed database hierarchy."""
        nodes = set(hierarchy.features(feature_set))
        if not nodes:
            raise ValueError(
                f"feature set {feature_set!r} is not in the hierarchy "
                f"(have: {', '.join(sorted(hierarchy.feature_sets()))})"
            )
        parent: dict[str, str] = {}
        children: dict[str, tuple[str, ...]] = {}
        for p in nodes:
            kids = tuple(sorted(hierarchy.children(feature_set, p)))
            if kids:
                children[p] = kids
                for c in kids:
                    parent[c] = p
        roots = sorted(n for n in nodes if n not in parent)
        if len(roots) == 1:
            root = roots[0]
        else:  # a forest — splice in a synthetic super-root (not seen for the v2 sets)
            root = "__root__"
            children[root] = tuple(roots)
            for r in roots:
                parent[r] = root
        depth: dict[str, int] = {root: 0}
        stack = [root]
        while stack:
            n = stack.pop()
            for c in children.get(n, ()):
                depth[c] = depth[n] + 1
                stack.append(c)
        return cls(parent=parent, children=children, depth=depth, root=root)


def _subtree_weights(weights: Mapping[str, float], tree: BinTree) -> dict[str, float]:
    """Total window bp at each node = its own label's bp + all descendants' bp.

    A label outside the tree (``novel``) contributes to itself and to the root, as if it
    were a direct child of the root.
    """
    sub: dict[str, float] = defaultdict(float)
    root = tree.root
    for label, w in weights.items():
        if w <= 0:
            continue
        sub[label] += w
        node = label
        while node != root:
            par = tree.parent.get(node)
            if par is None:  # out-of-tree label: attach under the root
                sub[root] += w
                break
            node = par
            sub[node] += w
    return sub


def descend(
    weights: Mapping[str, float],
    tree: BinTree,
    *,
    majority_fraction: float = DEFAULT_MAJORITY,
    scope: str = "node",
    novel_min_fraction: float = DEFAULT_NOVEL_MIN,
) -> str | None:
    """The hierarchical majority-rule call for one window's per-feature bp ``weights``.

    Returns the deepest node whose dominant child clears the majority bar at every step,
    or the flat-plurality label if the descent can't leave the root. ``None`` for an empty
    window.

    ``novel`` (k-mer-not-in-index) is gated on an **absolute** fraction: it is only reported when
    it covers at least ``novel_min_fraction`` of the window; otherwise it is dropped from the vote
    and the dominant *non-novel* feature is reported. Because novel positions are identical across
    a database's featuresets (it is an index property, not a featureset call), this absolute gate
    makes the binned-novel extent featureset-independent — so overlaying the binned featuresets
    yields ``novel:novel`` rather than spurious ``chrN:novel`` / ``novel:feature`` mixes that a
    relative plurality (whose competing feature differs per featureset) would produce.
    """
    total = sum(w for w in weights.values() if w > 0)
    if total <= 0:
        return None
    novel_w = weights.get(NOVEL, 0.0)
    if novel_w > 0:
        if novel_w >= novel_min_fraction * total:
            return NOVEL
        weights = {k: v for k, v in weights.items() if k != NOVEL}
        total = sum(w for w in weights.values() if w > 0)
        if total <= 0:
            return NOVEL
    sub = _subtree_weights(weights, tree)
    root = tree.root
    node = root
    while True:
        kids = list(tree.children.get(node, ()))
        if node == root:  # out-of-tree present labels (novel) vote as top-level leaves
            kids += [lbl for lbl in weights if lbl not in tree.depth and weights[lbl] > 0]
        kids = [c for c in kids if sub.get(c, 0.0) > 0.0]
        if not kids:
            break
        best = max(kids, key=lambda c: (sub.get(c, 0.0), tree.depth.get(c, 1), c))
        denom = total if scope == "window" else sub.get(node, total)
        if sub.get(best, 0.0) > majority_fraction * denom:
            node = best
        else:
            break
    if node == root:  # no top-level majority -> flat plurality (keeps boundaries sharp)
        node = max(
            (lbl for lbl in weights if weights[lbl] > 0),
            key=lambda lbl: (weights[lbl], tree.depth.get(lbl, 1), lbl),
        )
    return node


def _merge(runs: list[Interval], start: int, end: int, feature: str) -> None:
    """Append ``[start, end) -> feature``, coalescing with the previous run if equal."""
    if runs and runs[-1][2] == feature and runs[-1][1] == start:
        runs[-1] = (runs[-1][0], end, feature)
    else:
        runs.append((start, end, feature))


def _cap_gt(b0: float, bs: float, r0: float, rs: float) -> float:
    """Max integer ``dt >= 0`` with ``b0 + bs*dt > r0 + rs*dt`` (given it holds at dt=0)."""
    diff0, dcoeff = b0 - r0, bs - rs  # diff0 > 0
    if dcoeff >= 0:
        return _INF
    return max(0, math.ceil(diff0 / (-dcoeff) - _EPS) - 1)


def _cap_le(c0: float, cs: float, r0: float, rs: float) -> float:
    """Max integer ``dt >= 0`` with ``c0 + cs*dt <= r0 + rs*dt`` (given it holds at dt=0)."""
    diff0, dcoeff = c0 - r0, cs - rs  # diff0 <= 0
    if dcoeff <= 0:
        return _INF
    return max(0, math.floor((-diff0) / dcoeff + _EPS))


def _cap_zero(c0: float, cs: float) -> float:
    """Max integer ``dt >= 0`` keeping ``c0 + cs*dt > 0`` (``c0 > 0``, ``cs < 0``)."""
    return max(0, math.ceil(c0 / (-cs) - _EPS) - 1)


def _descent_run(
    counts: Mapping[str, float],
    tree: BinTree,
    *,
    anc_e: frozenset[str],
    anc_x: frozenset[str],
    majority_fraction: float,
    scope: str,
    novel_min_fraction: float = DEFAULT_NOVEL_MIN,
) -> tuple[str, float]:
    """The descent feature plus a safe lower bound on how long it stays constant.

    Within a segment each window count is linear in the step offset ``dt`` with slope
    ``[node under entering] - [node under leaving]`` (``anc_e``/``anc_x`` are the
    ancestor-or-self sets of the entering/leaving features, empty if that edge is inactive).
    Returns ``(feature, max_dt)`` where the call is provably constant for ``dt in [0, max_dt]``
    (``max_dt`` may be ``inf``). Bounds are conservative — recompute-and-merge keeps the
    result exact regardless — and mirror :func:`descend` at ``dt = 0``.
    """
    # The ``novel`` absolute-fraction gate (see :func:`descend`) is a non-hierarchical threshold;
    # rather than thread its crossing through the analytical bounds, recompute per base whenever
    # novel is present in, entering, or leaving the window. Non-novel windows (the bulk) keep the
    # window-independent fast path below.
    if counts.get(NOVEL, 0.0) > 0 or NOVEL in anc_e or NOVEL in anc_x:
        feat = descend(
            counts, tree, majority_fraction=majority_fraction, scope=scope,
            novel_min_fraction=novel_min_fraction,
        )
        assert feat is not None
        return feat, 0.0

    sub0 = _subtree_weights(counts, tree)
    root = tree.root

    def slope(n: str) -> int:
        return (1 if n in anc_e else 0) - (1 if n in anc_x else 0)

    total0 = float(sum(counts.values()))
    total_slope = slope(root)
    run: float = _INF
    node = root
    while True:
        kids = list(tree.children.get(node, ()))
        if node == root:
            # out-of-tree features (novel) vote as top-level leaves — include any already
            # present, plus the entering one, which is about to appear (so its arrival is capped)
            extra = {f for f in counts if f not in tree.depth and counts[f] > 0}
            extra |= {f for f in anc_e if f not in tree.depth and f != root}
            kids += sorted(extra)
        for c in kids:  # an absent child gaining the entering feature would change the choice
            if sub0.get(c, 0.0) <= 0.0 and slope(c) > 0:
                run = 0
        present = [c for c in kids if sub0.get(c, 0.0) > 0.0]
        if not present:
            break
        best = max(present, key=lambda c: (sub0[c], tree.depth.get(c, 1), c))
        if scope == "window":
            denom0, denom_slope = total0, float(total_slope)
        else:
            denom0, denom_slope = sub0.get(node, total0), float(slope(node))
        b0, bs = sub0[best], float(slope(best))
        if b0 > majority_fraction * denom0:  # descend into best
            run = min(run, _cap_gt(b0, bs, majority_fraction * denom0, majority_fraction * denom_slope))
            for c in present:
                if c == best:
                    continue
                run = min(run, _cap_arg(b0, bs, sub0[c], float(slope(c))))
                if slope(c) < 0:
                    run = min(run, _cap_zero(sub0[c], float(slope(c))))
            node = best
        else:  # stop at this node
            for c in present:
                run = min(
                    run, _cap_le(sub0[c], float(slope(c)), majority_fraction * denom0, majority_fraction * denom_slope)
                )
            break
    if node == root:  # plurality fallback (rare boundary positions): recompute each step
        node = max((f for f in counts if counts[f] > 0), key=lambda f: (counts[f], tree.depth.get(f, 1), f))
        run = 0
    return node, run


def _cap_arg(b0: float, bs: float, c0: float, cs: float) -> float:
    """Max integer ``dt >= 0`` keeping ``best`` ahead of child ``c`` (best wins at dt=0)."""
    if b0 > c0:
        return _cap_gt(b0, bs, c0, cs)  # stays strictly ahead
    # tie at dt=0 (best won the depth/name tiebreak): holds unless c pulls ahead
    return _INF if (bs - cs) >= 0 else 0


def bin_intervals(
    intervals: Sequence[Interval],
    tree: BinTree,
    *,
    window: int = DEFAULT_WINDOW,
    majority_fraction: float = DEFAULT_MAJORITY,
    scope: str = "node",
    novel_min_fraction: float = DEFAULT_NOVEL_MIN,
) -> list[Interval]:
    """Mode-filter a single sequence's intervals; return the smoothed partition of ``[0, L)``.

    ``intervals`` must be a gapless partition of ``[0, L)`` (C4), sorted by start. For each
    base ``i`` the window is ``[i - a, i + b]`` clipped to the sequence (``a = (W-1)//2``,
    ``b = W//2``); the output feature is :func:`descend` over that window's per-feature bp.

    **O(intervals), window-independent.** As ``i`` sweeps the sequence the window composition
    changes linearly between O(intervals) breakpoints (where the entering/leaving base crosses
    an interval boundary); within each such segment :func:`_descent_run` returns the call plus
    a safe bound on how far it stays constant, so the descent is evaluated O(1) times per
    segment instead of per base — the cost no longer grows with the window size.
    """
    if not intervals:
        return []
    a, b = (window - 1) // 2, window // 2
    length = intervals[-1][1]
    starts = [iv[0] for iv in intervals]
    ends = [iv[1] for iv in intervals]

    def label_at(base: int) -> str:
        return intervals[bisect.bisect_right(starts, base) - 1][2]

    def interval_of(base: int) -> int:
        return bisect.bisect_right(starts, base) - 1

    anc_cache: dict[str, frozenset[str]] = {}

    def ancestors(feature: str) -> frozenset[str]:
        cached = anc_cache.get(feature)
        if cached is None:
            seen = {feature}
            node = feature
            while node in tree.parent:
                node = tree.parent[node]
                seen.add(node)
            if feature not in tree.depth:  # out-of-tree (novel): attach to root
                seen.add(tree.root)
            cached = frozenset(seen)
            anc_cache[feature] = cached
        return cached

    # initial window(0) counts = features over [0, min(L-1, b)]
    counts: dict[str, float] = defaultdict(float)
    hi0 = min(length - 1, b)
    k = 0
    while k < len(intervals) and starts[k] <= hi0:
        portion = min(ends[k] - 1, hi0) - starts[k] + 1
        if portion > 0:
            counts[intervals[k][2]] += portion
        k += 1

    runs: list[Interval] = []
    i = 0
    while i <= length - 1:
        if i == length - 1:  # last position: one window, no forward step
            feat = descend(
                counts, tree, majority_fraction=majority_fraction, scope=scope,
                novel_min_fraction=novel_min_fraction,
            )
            assert feat is not None
            _merge(runs, i, i + 1, feat)
            break

        enter, leave = i + 1 + b, i - a
        e_act, l_act = enter <= length - 1, leave >= 0
        anc_e = ancestors(label_at(enter)) if e_act else frozenset()
        anc_x = ancestors(label_at(leave)) if l_act else frozenset()
        E = label_at(enter) if e_act else None
        X = label_at(leave) if l_act else None

        # next position where the step's entering/leaving label (or clip status) changes
        seg_end = length - 1
        if e_act:
            seg_end = min(seg_end, ends[interval_of(enter)] - 1 - b)
        # leaving label changes at el + a; if leaving is inactive it activates at i = a
        seg_end = min(seg_end, ends[interval_of(leave)] + a if l_act else a)

        feat, run = _descent_run(
            counts, tree, anc_e=anc_e, anc_x=anc_x, majority_fraction=majority_fraction,
            scope=scope, novel_min_fraction=novel_min_fraction,
        )
        d = min(run, seg_end - 1 - i, length - 1 - i)  # constant call, within segment + sequence
        d = max(0, int(d))
        _merge(runs, i, i + d + 1, feat)

        steps = d + 1
        if e_act and E is not None:
            counts[E] = counts.get(E, 0.0) + steps
        if l_act and X is not None:
            counts[X] = counts.get(X, 0.0) - steps
            if counts[X] <= 0:
                del counts[X]
        i += steps
    return runs


def bin_intervals_naive(
    intervals: Sequence[Interval],
    tree: BinTree,
    *,
    window: int = DEFAULT_WINDOW,
    majority_fraction: float = DEFAULT_MAJORITY,
    scope: str = "node",
    novel_min_fraction: float = DEFAULT_NOVEL_MIN,
) -> list[Interval]:
    """Reference implementation: evaluate :func:`descend` at every base (for tests)."""
    if not intervals:
        return []
    a, b = (window - 1) // 2, window // 2
    length = intervals[-1][1]
    labels = [""] * length
    for s, e, f in intervals:
        for x in range(s, e):
            labels[x] = f
    runs: list[Interval] = []
    for i in range(length):
        lo, hi = max(0, i - a), min(length - 1, i + b)
        counts: dict[str, float] = defaultdict(float)
        for x in range(lo, hi + 1):
            counts[labels[x]] += 1
        feature = descend(
            counts, tree, majority_fraction=majority_fraction, scope=scope,
            novel_min_fraction=novel_min_fraction,
        )
        assert feature is not None
        _merge(runs, i, i + 1, feature)
    return runs


def bin_sequence(
    intervals: Sequence[Interval],
    tree: BinTree,
    *,
    window: int = DEFAULT_WINDOW,
    majority_fraction: float = DEFAULT_MAJORITY,
    scope: str = "node",
    novel_min_fraction: float = DEFAULT_NOVEL_MIN,
) -> list[Interval]:
    """:func:`bin_intervals` for a sequence whose first interval may not start at 0."""
    if not intervals:
        return []
    off = intervals[0][0]
    shifted = [(s - off, e - off, f) for s, e, f in intervals]
    binned = bin_intervals(
        shifted, tree, window=window, majority_fraction=majority_fraction, scope=scope,
        novel_min_fraction=novel_min_fraction,
    )
    return [(s + off, e + off, f) for s, e, f in binned]
