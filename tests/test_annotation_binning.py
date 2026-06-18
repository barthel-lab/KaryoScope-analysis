"""Tests for the hierarchical mode-filter (bin-annotations) core + CLI.

The descent tests use the real ``tests/data/hierarchy.tsv`` region tree, where ``aSat`` and
``bSat`` are siblings under ``centromeric`` and ``arm`` is a separate child of the root —
the exact shape the design discussion turned on.
"""

from __future__ import annotations

import itertools
import random
from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import annotation_binning as binning
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


@pytest.fixture(scope="module")
def region_tree() -> binning.BinTree:
    hierarchy = FeatureHierarchy.from_tsv(HIERARCHY_TSV)
    return binning.BinTree.from_hierarchy(hierarchy, "region")


# ------------------------------------------------------------------ tree construction
def test_bin_tree_shape(region_tree: binning.BinTree):
    assert region_tree.root == "categorized"
    assert region_tree.parent["aSat"] == "centromeric"
    assert region_tree.parent["bSat"] == "centromeric"
    assert region_tree.parent["arm"] == "categorized"
    assert region_tree.depth["categorized"] == 0
    assert region_tree.depth["centromeric"] == 1
    assert region_tree.depth["aSat"] == 2


def test_bin_tree_unknown_feature_set():
    hierarchy = FeatureHierarchy.from_tsv(HIERARCHY_TSV)
    with pytest.raises(ValueError, match="not in the hierarchy"):
        binning.BinTree.from_hierarchy(hierarchy, "nonexistent")


# ------------------------------------------------------------------ descent semantics
def test_descend_aggregates_siblings(region_tree: binning.BinTree):
    # aSat 20 + bSat 40 + arm 40: siblings sum to 60 (majority) -> descend into centromeric;
    # node-scope then descends to bSat (40 > 50% of 60); window-scope stops at centromeric.
    weights = {"aSat": 20.0, "bSat": 40.0, "arm": 40.0}
    assert binning.descend(weights, region_tree, scope="node") == "bSat"
    assert binning.descend(weights, region_tree, scope="window") == "centromeric"


def test_descend_node_scope_can_over_descend(region_tree: binning.BinTree):
    # arm 49 / aSat 40 / bSat 11: near-even top split. node-scope reports aSat (below the
    # window plurality, arm); window-scope conservatively reports centromeric.
    weights = {"arm": 49.0, "aSat": 40.0, "bSat": 11.0}
    assert binning.descend(weights, region_tree, scope="node") == "aSat"
    assert binning.descend(weights, region_tree, scope="window") == "centromeric"


def test_descend_root_falls_back_to_plurality(region_tree: binning.BinTree):
    # 50/50 unrelated split: no top-level majority -> flat plurality, tie broken to the
    # deeper/more-specific label (aSat at depth 2 over arm at depth 1).
    weights = {"arm": 50.0, "aSat": 50.0}
    assert binning.descend(weights, region_tree, scope="node") == "aSat"
    assert binning.descend(weights, region_tree, scope="window") == "aSat"


def test_descend_single_label(region_tree: binning.BinTree):
    assert binning.descend({"HSat3": 101.0}, region_tree) == "HSat3"


def test_descend_novel_as_top_level_leaf(region_tree: binning.BinTree):
    # novel wins only when it covers >= novel_min_fraction (default 0.5) of the window.
    assert binning.descend({"novel": 60.0, "aSat": 40.0}, region_tree, scope="node") == "novel"


def test_descend_novel_gated_below_threshold(region_tree: binning.BinTree):
    # 40% novel < 50%: novel is dropped from the vote; the dominant non-novel feature wins
    # (60% aSat -> aSat), NOT novel (which a plurality vote would have picked at tau=0).
    assert binning.descend({"novel": 40.0, "aSat": 60.0}, region_tree,
                           scope="node", majority_fraction=0.0) == "aSat"
    # exactly at the threshold (50%) novel still wins.
    assert binning.descend({"novel": 50.0, "aSat": 50.0}, region_tree, scope="node") == "novel"
    # the gate is configurable.
    assert binning.descend({"novel": 40.0, "aSat": 60.0}, region_tree,
                           scope="node", novel_min_fraction=0.3) == "novel"


def test_bin_intervals_novel_gate_fast_matches_naive(region_tree: binning.BinTree):
    # A window mostly non-novel with scattered novel: novel must not win, and the fast path's
    # per-base novel recompute must match the reference implementation.
    intervals = [(0, 30, "aSat"), (30, 45, "novel"), (45, 100, "bSat")]
    fast = binning.bin_intervals(intervals, region_tree, window=51, majority_fraction=0.0)
    naive = binning.bin_intervals_naive(intervals, region_tree, window=51, majority_fraction=0.0)
    assert fast == naive


def test_descend_majority_fraction_knob(region_tree: binning.BinTree):
    # aSat 40 / bSat 35 (under centromeric) / arm 25: centromeric=75 clears the bar at the
    # root for both tau below. Within centromeric (75): aSat 40 > 0.5*75=37.5 -> aSat, but
    # 40 < 0.6*75=45 -> stops at centromeric. So a higher bar is more conservative.
    weights = {"aSat": 40.0, "bSat": 35.0, "arm": 25.0}
    assert binning.descend(weights, region_tree, scope="node", majority_fraction=0.5) == "aSat"
    assert (
        binning.descend(weights, region_tree, scope="node", majority_fraction=0.6) == "centromeric"
    )


# ------------------------------------------------------------------ interval smoothing
def test_bin_smooths_tiny_fragment(region_tree: binning.BinTree):
    # A 5 bp bSat fragment inside a long aSat region is absorbed by a 51 bp window.
    intervals = [(0, 100, "aSat"), (100, 105, "bSat"), (105, 200, "aSat")]
    out = binning.bin_intervals(intervals, region_tree, window=51)
    assert out == [(0, 200, "aSat")]


def test_bin_preserves_length(region_tree: binning.BinTree):
    intervals = [(0, 100, "aSat"), (100, 105, "bSat"), (105, 200, "arm")]
    out = binning.bin_intervals(intervals, region_tree, window=51)
    assert out[0][0] == 0 and out[-1][1] == 200
    # gapless partition preserved
    for prev, cur in itertools.pairwise(out):
        assert prev[1] == cur[0]


def test_bin_sequence_offset(region_tree: binning.BinTree):
    # A sequence not starting at 0 keeps its absolute coordinates.
    intervals = [(1000, 1100, "aSat"), (1100, 1105, "bSat"), (1105, 1200, "aSat")]
    out = binning.bin_sequence(intervals, region_tree, window=51)
    assert out == [(1000, 1200, "aSat")]


def test_bin_empty(region_tree: binning.BinTree):
    assert binning.bin_intervals([], region_tree) == []


# ------------------------------------------------------------------ fast == naive (property)
_FEATURES = ["aSat", "bSat", "HSat3", "cenSat", "arm", "ct", "rDNA", "novel"]


def _random_partition(rng: random.Random, length: int) -> list:
    intervals = []
    pos = 0
    while pos < length:
        seg = min(rng.randint(1, 12), length - pos)
        intervals.append((pos, pos + seg, rng.choice(_FEATURES)))
        pos += seg
    return intervals


@pytest.mark.parametrize("scope", ["node", "window"])
def test_fast_matches_naive(region_tree: binning.BinTree, scope: str):
    rng = random.Random(20240608)
    for _ in range(600):
        length = rng.randint(1, 80)
        intervals = _random_partition(rng, length)
        window = rng.choice([1, 2, 3, 5, 11, 21, 31, 50, 101])
        tau = rng.choice([0.0, 0.3, 0.5, 0.7, 1.0])
        fast = binning.bin_intervals(
            intervals, region_tree, window=window, majority_fraction=tau, scope=scope
        )
        naive = binning.bin_intervals_naive(
            intervals, region_tree, window=window, majority_fraction=tau, scope=scope
        )
        assert fast == naive, (intervals, window, tau, scope)


# ------------------------------------------------------------------ strided engine (step > 1)
@pytest.mark.parametrize("scope", ["node", "window"])
def test_strided_step1_matches_naive(region_tree: binning.BinTree, scope: str):
    # At step=1 the strided engine samples the window centered on every base -- it must
    # reproduce the per-base reference exactly, validating the incremental window arithmetic.
    rng = random.Random(20240609)
    for _ in range(400):
        length = rng.randint(1, 80)
        intervals = _random_partition(rng, length)
        window = rng.choice([1, 3, 11, 21, 51])
        tau = rng.choice([0.0, 0.5, 1.0])
        strided = binning.bin_intervals_strided(
            intervals, region_tree, window=window, step=1, majority_fraction=tau, scope=scope
        )
        naive = binning.bin_intervals_naive(
            intervals, region_tree, window=window, majority_fraction=tau, scope=scope
        )
        assert strided == naive, (intervals, window, tau, scope)


def test_strided_preserves_gapless_partition(region_tree: binning.BinTree):
    intervals = [(0, 100, "aSat"), (100, 105, "bSat"), (105, 263, "arm")]
    out = binning.bin_intervals_strided(intervals, region_tree, window=51, step=10)
    assert out[0][0] == 0 and out[-1][1] == 263  # exact [0, L) coverage, odd length vs step
    for prev, cur in itertools.pairwise(out):
        assert prev[1] == cur[0]


def test_strided_boundaries_snap_to_step_grid(region_tree: binning.BinTree):
    # A clean aSat/arm boundary at 1000; with step=100 the called boundary lands on a grid
    # multiple (within one step of the true boundary), unlike the exact step=1 engine.
    intervals = [(0, 1000, "aSat"), (1000, 2000, "arm")]
    out = binning.bin_intervals_strided(intervals, region_tree, window=51, step=100)
    boundaries = [e for _, e, _ in out[:-1]]
    assert all(boundary % 100 == 0 for boundary in boundaries)
    assert [f for _, _, f in out] == ["aSat", "arm"]


def test_strided_rejects_bad_step(region_tree: binning.BinTree):
    with pytest.raises(ValueError, match="step must be >= 1"):
        binning.bin_intervals_strided([(0, 10, "aSat")], region_tree, window=11, step=0)


def test_bin_sequence_dispatches_on_step(region_tree: binning.BinTree):
    # step>1 path is reachable through the public sequence wrapper and keeps absolute coords.
    intervals = [(500, 600, "aSat"), (600, 605, "bSat"), (605, 800, "aSat")]
    out = binning.bin_sequence(intervals, region_tree, window=51, step=10)
    assert out == [(500, 800, "aSat")]


def test_descend_tau_zero_always_reaches_a_leaf(region_tree: binning.BinTree):
    # tau=0 descends into the heaviest subtree at every level -> a specific leaf, never an
    # internal/ambiguous node. aSat 20 / bSat 40 / arm 40 -> centromeric(60) -> bSat(40).
    leaves = {"aSat", "bSat", "HSat3", "arm", "ct", "rDNA", "p_arm", "q_arm", "HSat1A"}
    assert binning.descend({"aSat": 20.0, "bSat": 40.0, "arm": 40.0}, region_tree,
                           majority_fraction=0.0) == "bSat"
    for weights in ({"aSat": 30.0, "bSat": 30.0, "arm": 40.0}, {"HSat3": 10.0, "p_arm": 9.0}):
        assert binning.descend(weights, region_tree, majority_fraction=0.0) in leaves


# ------------------------------------------------------------------ CLI end-to-end
def test_bin_annotations_cli(cli_runner, tmp_path: Path):
    bed = tmp_path / "region.bed"
    # two reads; each has a tiny fragment that should be smoothed away
    bed.write_text(
        "r1\t0\t100\taSat\nr1\t100\t105\tbSat\nr1\t105\t200\taSat\n"
        "r2\t0\t120\tarm\nr2\t120\t124\tct\nr2\t124\t240\tarm\n"
    )
    out = tmp_path / "region.binned.bed"
    res = cli_runner.invoke(
        main,
        [
            "bin-annotations",
            "--input",
            str(bed),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--feature-set",
            "region",
            "--window",
            "51",
            "-o",
            str(out),
        ],
    )
    assert res.exit_code == 0, res.output
    lines = [ln.split("\t") for ln in out.read_text().splitlines()]
    # both fragments absorbed -> one interval per read, length preserved
    assert lines == [["r1", "0", "200", "aSat"], ["r2", "0", "240", "arm"]]


def test_bin_annotations_cli_rejects_unknown_feature(cli_runner, tmp_path: Path):
    bed = tmp_path / "bad.bed"
    bed.write_text("r1\t0\t100\taSat\nr1\t100\t200\tnot_a_feature\n")
    out = tmp_path / "out.bed"
    res = cli_runner.invoke(
        main,
        [
            "bin-annotations",
            "--input",
            str(bed),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--feature-set",
            "region",
            "-o",
            str(out),
        ],
    )
    assert res.exit_code != 0
    assert "not_a_feature" in res.output


def test_bin_annotations_cli_rejects_unknown_feature_set(cli_runner, tmp_path: Path):
    bed = tmp_path / "r.bed"
    bed.write_text("r1\t0\t100\taSat\n")
    out = tmp_path / "out.bed"
    res = cli_runner.invoke(
        main,
        [
            "bin-annotations",
            "--input",
            str(bed),
            "--hierarchy",
            str(HIERARCHY_TSV),
            "--feature-set",
            "bogus",
            "-o",
            str(out),
        ],
    )
    assert res.exit_code != 0
    assert "not in the hierarchy" in res.output
