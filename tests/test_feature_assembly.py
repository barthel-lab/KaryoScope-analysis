"""Tests for the Engine B overlap graph + clustering."""

from __future__ import annotations

import pytest

from karyoscope_analysis.core import feature_assembly as asm
from karyoscope_analysis.core.feature_align import reverse_segments

EXACT = lambda x, y: 1.0 if x == y else -1.0  # noqa: E731
PARAMS = dict(
    sub_score=EXACT, gap_factor=0.01, match_score=1.0, min_overlap_bp=50, min_identity=0.9
)


# ----------------------------------------------------------------- overlap graph
def test_dovetail_makes_an_edge():
    reads = {
        "x": [("P", 30), ("Q", 100), ("R", 80)],
        "y": [("Q", 100), ("R", 80), ("S", 40)],
    }
    edges = asm.build_overlap_graph(reads, **PARAMS)
    assert len(edges) == 1
    assert edges[0].kind == "dovetail"
    assert not edges[0].flipped
    assert edges[0].identity == pytest.approx(1.0)


def test_unrelated_reads_make_no_edge():
    reads = {"x": [("X", 100)], "y": [("Y", 100)]}
    assert asm.build_overlap_graph(reads, **PARAMS) == []


def test_internal_only_match_is_rejected():
    # A and B share only an internal X block (a "repeat") -> not a proper overlap.
    reads = {
        "x": [("P", 100), ("X", 60), ("Q", 100)],
        "y": [("R", 100), ("X", 60), ("S", 100)],
    }
    assert asm.build_overlap_graph(reads, **PARAMS) == []


def test_short_overlap_is_rejected():
    # The shared block is only 10 bp < min_overlap_bp (50).
    reads = {"x": [("P", 100), ("Q", 10)], "y": [("Q", 10), ("S", 100)]}
    assert asm.build_overlap_graph(reads, **PARAMS) == []


def test_reverse_orientation_makes_a_flipped_edge():
    base = [("P", 30), ("Q", 100), ("R", 80)]
    reads = {"x": base, "y": reverse_segments([("Q", 100), ("R", 80), ("S", 40)])}
    edges = asm.build_overlap_graph(reads, **PARAMS)
    assert len(edges) == 1
    assert edges[0].flipped


# ----------------------------------------------------------------- parity union-find
def test_parity_dsu_orientation_propagation():
    dsu = asm._ParityDSU()
    for n in ("a", "b", "c"):
        dsu.add(n)
    dsu.union("a", "b", 0)  # same orientation
    dsu.union("b", "c", 1)  # c flipped relative to b
    ra, pa = dsu.find("a")
    rc, pc = dsu.find("c")
    assert ra == rc
    assert (pa ^ pc) == 1  # c is flipped relative to a


def test_parity_dsu_detects_conflict():
    dsu = asm._ParityDSU()
    for n in ("a", "b", "c"):
        dsu.add(n)
    dsu.union("a", "b", 0)
    dsu.union("b", "c", 0)
    assert dsu.union("a", "c", 1) is False  # implies a^c=1 but already a^c=0
    root, _ = dsu.find("a")
    assert root in dsu.conflicts


# ----------------------------------------------------------------- clustering
def test_transitive_chain_forms_one_cluster():
    # x--y (dovetail) and y--z (dovetail), but x and z don't overlap -> one cluster of 3.
    reads = {
        "x": [("A", 100), ("B", 100)],
        "y": [("B", 100), ("C", 100)],
        "z": [("C", 100), ("D", 100)],
    }
    clusters, edges = asm.assemble(reads, **PARAMS)
    assert {e.kind for e in edges} == {"dovetail"}
    assert len(clusters) == 1
    assert clusters[0].size == 3
    assert set(clusters[0].members) == {"x", "y", "z"}


def test_unrelated_read_is_its_own_singleton():
    reads = {
        "x": [("A", 100), ("B", 100)],
        "y": [("B", 100), ("C", 100)],
        "lone": [("Z", 100), ("W", 100)],
    }
    clusters, _ = asm.assemble(reads, **PARAMS)
    sizes = sorted(c.size for c in clusters)
    assert sizes == [1, 2]
    singleton = next(c for c in clusters if c.size == 1)
    assert singleton.members == ("lone",)
    assert singleton.seed == "lone"


def test_reversed_member_marked_relative_to_seed():
    base = [("A", 200), ("B", 100), ("C", 100)]  # longest -> seed
    reads = {"seed": base, "rev": reverse_segments(base)}
    clusters, _ = asm.assemble(reads, **PARAMS)
    assert len(clusters) == 1
    c = clusters[0]
    assert c.seed == "seed"
    assert c.reversed_relative_to_seed["seed"] is False
    assert c.reversed_relative_to_seed["rev"] is True


# ----------------------------------------------------------------- consensus
def _consensus(reads, cluster):
    return asm.cluster_consensus(reads, cluster, sub_score=EXACT, gap_factor=0.01)


def test_consensus_singleton_is_the_read_itself():
    reads = {"x": [("A", 100), ("B", 50)]}
    clusters, _ = asm.assemble(reads, **PARAMS)
    cons = _consensus(reads, clusters[0])
    assert cons.segments() == [("A", 100), ("B", 50)]
    assert [p.support for p in cons.positions] == [1, 1]
    assert [p.coverage for p in cons.positions] == [1, 1]


def test_consensus_counts_agreement_support():
    reads = {
        "seed": [("A", 100), ("B", 100), ("C", 100)],  # longest tie -> "seed" wins by id
        "m": [("B", 100), ("C", 100), ("D", 100)],  # overlaps B, C
    }
    clusters, _ = asm.assemble(reads, **PARAMS)
    cons = _consensus(reads, clusters[0])
    assert cons.seed == "seed"
    assert [p.feature for p in cons.positions] == ["A", "B", "C"]
    assert [p.support for p in cons.positions] == [1, 2, 2]  # A seed-only; B,C agreed
    assert [p.coverage for p in cons.positions] == [1, 2, 2]
    assert cons.segments() == [("A", 100), ("B", 100), ("C", 100)]


def test_consensus_orients_reversed_members():
    seed = [("A", 300), ("B", 100), ("C", 100)]  # 500 bp -> unambiguous seed
    reads = {"seed": seed, "rev": reverse_segments([("A", 300), ("B", 100)])}
    clusters, _ = asm.assemble(reads, **PARAMS)
    cons = _consensus(reads, clusters[0])
    assert cons.seed == "seed"
    # rev (oriented back) covers A and B but not C
    assert [p.support for p in cons.positions] == [2, 2, 1]
