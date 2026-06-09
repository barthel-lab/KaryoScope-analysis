"""Tests for the Engine B overlap graph + clustering."""

from __future__ import annotations

import pytest

from karyoscope_analysis.core import feature_assembly as asm
from karyoscope_analysis.core.feature_align import reverse_segments

EXACT = lambda x, y: 1.0 if x == y else -1.0  # noqa: E731
PARAMS = dict(
    sub_score=EXACT, gap_factor=0.01, match_score=1.0, min_overlap_bp=50, min_identity=0.9
)


# ----------------------------------------------------------------- weighting (structural layer)
def test_weight_looks_up_structural_layer():
    w = {"q_arm": 0.03, "canonical_telomere": 0.5}
    assert asm._weight(w, "chr1:q_arm") == 0.03  # composite -> structural lookup
    assert asm._weight(w, "q_arm") == 0.03  # bare label still works
    assert asm._weight(w, "chrX:canonical_telomere") == 0.5
    assert asm._weight(w, "chr2:p_arm") == 1.0  # unknown -> default 1.0
    assert asm._weight(None, "chr1:q_arm") == 1.0  # no weighting


def test_idf_weights_keyed_by_structural_layer():
    reads = {
        "a": [("chr1:q_arm", 100), ("chr1:aSat", 50)],
        "b": [("chr2:q_arm", 100), ("chr2:bSat", 50)],
    }
    w = asm.idf_weights(reads)
    assert "q_arm" in w and "chr1:q_arm" not in w  # structural keys, not composites
    assert w["q_arm"] < w["aSat"]  # q_arm in both reads (ubiquitous), aSat in one


# ----------------------------------------------------------------- distinctive overlap
def test_min_distinctive_bp_rejects_filler_only_overlap():
    # x's suffix == y's prefix == a big down-weighted "arm" block (a filler dovetail).
    reads = {
        "x": [("aSat", 400), ("arm", 5000)],
        "y": [("arm", 5000), ("bSat", 400)],
    }
    weight = {"arm": 0.03, "aSat": 0.5, "bSat": 0.5}  # arm = filler
    common = dict(
        sub_score=EXACT, gap_factor=0.01, min_overlap_bp=1, min_identity=0.5, weight=weight
    )
    # without the distinctive criterion the shared arm alone makes an edge
    assert asm.build_overlap_graph(reads, **common)
    # requiring 100 bp of matched distinctive content (weight >= 0.15) -> arm doesn't count -> no edge
    assert not asm.build_overlap_graph(reads, min_distinctive_bp=100, **common)


def test_min_distinctive_bp_with_filler_set():
    # filler-based distinctiveness: a shared telomere (genome-rare, so high weight) is filler and
    # can't make an edge, but a shared satellite can.
    telo = {
        "x": [("aSat", 400), ("canonical_telomere", 5000)],
        "y": [("canonical_telomere", 5000), ("bSat", 400)],
    }
    filler = frozenset({"canonical_telomere"})
    common = dict(sub_score=EXACT, gap_factor=0.01, min_overlap_bp=1, min_identity=0.5)
    assert not asm.build_overlap_graph(telo, min_distinctive_bp=100, filler_features=filler, **common)
    sat = {  # now the shared block is a satellite (not filler) -> edge
        "x": [("ct", 400), ("aSat", 5000)],
        "y": [("aSat", 5000), ("ct", 400)],
    }
    assert asm.build_overlap_graph(sat, min_distinctive_bp=100, filler_features=filler, **common)


# ----------------------------------------------------------------- blocking index
def test_block_min_bp_only_compares_reads_sharing_a_major_feature():
    # x & y share a big "Q" block (dovetail); z shares nothing major with either.
    reads = {
        "x": [("P", 2000), ("Q", 3000)],
        "y": [("Q", 3000), ("R", 2000)],
        "z": [("S", 4000), ("T", 4000)],
    }
    pairs = asm._candidate_pairs(["x", "y", "z"], reads, block_min_bp=2500)
    # only (x, y) share a >= 2500 bp feature (Q); z is isolated
    assert pairs == [(0, 1)]
    # the blocked graph still finds the x-y dovetail edge, and z makes none
    edges = asm.build_overlap_graph(
        reads, sub_score=EXACT, gap_factor=0.01, min_overlap_bp=1, min_identity=0.5,
        block_min_bp=2500,
    )
    assert {tuple(sorted((e.a, e.b))) for e in edges} == {("x", "y")}
    # same edges as all-vs-all here (z shares nothing alignable anyway)
    full = asm.build_overlap_graph(
        reads, sub_score=EXACT, gap_factor=0.01, min_overlap_bp=1, min_identity=0.5
    )
    assert {tuple(sorted((e.a, e.b))) for e in edges} == {tuple(sorted((e.a, e.b))) for e in full}


# ----------------------------------------------------------------- community detection
def _edge(a, b, w, flipped=False):
    return asm.OverlapEdge(a, b, w, 1.0, w, "dovetail", flipped)


def test_label_propagation_splits_bridged_cliques():
    # two triangles {a,b,c} and {x,y,z}, strongly internally connected, joined by ONE weak
    # bridge edge c-x (a noisy hub link). LP should keep two communities; CC would merge them.
    strong, weak = 10000.0, 1100.0
    edges = [
        _edge("a", "b", strong), _edge("a", "c", strong), _edge("b", "c", strong),
        _edge("x", "y", strong), _edge("x", "z", strong), _edge("y", "z", strong),
        _edge("c", "x", weak),  # the bridge
    ]
    label = asm._label_propagation(["a", "b", "c", "x", "y", "z"], edges)
    assert label["a"] == label["b"] == label["c"]
    assert label["x"] == label["y"] == label["z"]
    assert label["a"] != label["x"]  # two communities, not one


def test_cluster_reads_communities_vs_components():
    reads = {n: [("P", 1000)] for n in ("a", "b", "c", "x", "y", "z")}  # lengths only matter for seed
    strong, weak = 10000.0, 1100.0
    edges = [
        _edge("a", "b", strong), _edge("a", "c", strong), _edge("b", "c", strong),
        _edge("x", "y", strong), _edge("x", "z", strong), _edge("y", "z", strong),
        _edge("c", "x", weak),
    ]
    # connected components: one cluster of 6
    assert [c.size for c in asm.cluster_reads(reads, edges)] == [6]
    # communities: two clusters of 3
    assert sorted(c.size for c in asm.cluster_reads(reads, edges, communities=True)) == [3, 3]


# ----------------------------------------------------------------- parallel == serial
def test_workers_match_serial():
    reads = {
        "a": [("P", 30), ("Q", 100), ("R", 80)],
        "b": [("Q", 100), ("R", 80), ("S", 40)],
        "c": [("R", 80), ("S", 40), ("T", 50)],
        "d": [("X", 90), ("Y", 90)],
        "e": [("Y", 90), ("Z", 90)],
    }
    common = dict(sub_score=EXACT, gap_factor=0.01, min_overlap_bp=50, min_identity=0.9)
    serial = asm.build_overlap_graph(reads, **common, workers=1)
    parallel = asm.build_overlap_graph(reads, **common, workers=2)
    key = lambda es: sorted((e.a, e.b, round(e.score, 6), e.flipped) for e in es)  # noqa: E731
    assert key(serial) == key(parallel)
    # and with the blocking index on
    assert key(asm.build_overlap_graph(reads, **common, block_min_bp=50, workers=1)) == key(
        asm.build_overlap_graph(reads, **common, block_min_bp=50, workers=2)
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


# ----------------------------------------------------------------- feature weighting (v2)
def test_idf_weights():
    reads = {
        "r1": [("X", 100), ("Y", 100)],
        "r2": [("X", 100), ("Z", 100)],
        "r3": [("X", 100), ("W", 100)],
    }
    w = asm.idf_weights(reads, floor=0.1)
    assert w["X"] == pytest.approx(0.1)  # in all 3 reads -> floored
    assert w["Y"] == pytest.approx(2 / 3)  # in 1 of 3 -> 1 - 1/3
    assert w["Z"] == pytest.approx(2 / 3)


def test_weighting_suppresses_ubiquitous_overlap():
    # x and y dovetail only through a ubiquitous "REP" block.
    reads = {"x": [("DISTX", 2000), ("REP", 2000)], "y": [("REP", 2000), ("DISTY", 2000)]}
    base = dict(
        sub_score=EXACT, gap_factor=0.01, match_score=1.0, min_overlap_bp=500, min_identity=0.9
    )
    assert len(asm.build_overlap_graph(reads, **base)) == 1  # uniform: 2000 bp -> edge
    weight = {"REP": 0.05, "DISTX": 1.0, "DISTY": 1.0}
    # weighted overlap 0.05 * 2000 = 100 < 500 -> rejected
    assert asm.build_overlap_graph(reads, **base, weight=weight) == []


def test_weighting_keeps_distinctive_overlap():
    reads = {"x": [("A", 2000), ("DIST", 2000)], "y": [("DIST", 2000), ("B", 2000)]}
    base = dict(
        sub_score=EXACT, gap_factor=0.01, match_score=1.0, min_overlap_bp=500, min_identity=0.9
    )
    weight = {"DIST": 1.0, "A": 1.0, "B": 1.0}
    assert len(asm.build_overlap_graph(reads, **base, weight=weight)) == 1


# ----------------------------------------------------------------- consensus-coordinate layout
def _layout(reads, cluster):
    return asm.consensus_layout(reads, cluster, sub_score=EXACT, gap_factor=0.01)


def test_consensus_layout_singleton_is_the_read_itself():
    reads = {"x": [("A", 100), ("B", 50)]}
    clusters, _ = asm.assemble(reads, **PARAMS)
    lo = _layout(reads, clusters[0])
    assert lo.seed == "x" and lo.width == 150
    assert list(lo.placed[0].segments) == [(0, 100, "A"), (100, 150, "B")]
    assert [(p.start, p.end, p.feature, p.support, p.coverage) for p in lo.consensus] == [
        (0, 100, "A", 1, 1),
        (100, 150, "B", 1, 1),
    ]


def test_consensus_layout_stacks_features_and_spans_union():
    # member's B,C must land *under* the seed's B,C (stacking), and its D overhang must extend
    # the consensus span beyond the seed (union span), not be truncated.
    reads = {
        "seed": [("A", 100), ("B", 100), ("C", 100)],  # 0..300
        "m": [("B", 100), ("C", 100), ("D", 100)],  # B,C align to seed; D overhangs
    }
    clusters, _ = asm.assemble(reads, **PARAMS)
    lo = _layout(reads, clusters[0])
    assert lo.seed == "seed" and lo.width == 400
    placed = {r.read_id: r for r in lo.placed}
    assert list(placed["seed"].segments) == [(0, 100, "A"), (100, 200, "B"), (200, 300, "C")]
    assert list(placed["m"].segments) == [(100, 200, "B"), (200, 300, "C"), (300, 400, "D")]
    # consensus spans the union (incl. D); B,C have support 2 (both reads agree)
    assert [(p.feature, p.support, p.coverage) for p in lo.consensus] == [
        ("A", 1, 1), ("B", 2, 2), ("C", 2, 2), ("D", 1, 1)
    ]


def test_progressive_layout_places_transitive_members():
    # a chain seed -- m1 -- m2 where m2 overlaps m1 (on S) but NOT the seed.
    reads = {
        "seed": [("P", 100), ("Q", 100), ("R", 100)],
        "m1": [("Q", 100), ("R", 100), ("S", 100)],  # overlaps seed (Q,R)
        "m2": [("S", 100), ("T", 100), ("U", 100)],  # overlaps m1 (S), not the seed
    }
    clusters, edges = asm.assemble(
        reads, sub_score=EXACT, gap_factor=0.01, min_overlap_bp=50, min_identity=0.5
    )
    assert len(clusters) == 1  # one connected chain
    neighbors: dict[str, list[str]] = {}
    for e in edges:
        neighbors.setdefault(e.a, []).append(e.b)
        neighbors.setdefault(e.b, []).append(e.a)

    def sstart(layout, rid):
        seg = {r.read_id: r for r in layout.placed}[rid]
        return next(s for s, _e, f in seg.segments if f == "S")

    # progressive: m2 is placed via m1, so its S stacks on m1's S
    prog = asm.consensus_layout(
        reads, clusters[0], neighbors=neighbors, sub_score=EXACT, gap_factor=0.01
    )
    assert sstart(prog, "m1") == sstart(prog, "m2")
    # star (no neighbors): m2 can't align to the seed -> own coords -> NOT stacked
    star = asm.consensus_layout(reads, clusters[0], sub_score=EXACT, gap_factor=0.01)
    assert sstart(star, "m1") != sstart(star, "m2")


def test_consensus_layout_orients_reversed_members():
    seed = [("A", 300), ("B", 100), ("C", 100)]  # 500 bp -> unambiguous seed
    reads = {"seed": seed, "rev": reverse_segments([("A", 300), ("B", 100)])}
    clusters, _ = asm.assemble(reads, **PARAMS)
    lo = _layout(reads, clusters[0])
    assert lo.seed == "seed"
    placed = {r.read_id: r for r in lo.placed}
    assert placed["rev"].reversed is True
    # oriented back, rev's A,B stack under the seed's A,B
    assert list(placed["rev"].segments) == [(0, 300, "A"), (300, 400, "B")]
    # consensus: A,B agreed (support 2), C seed-only (support 1)
    assert [p.support for p in lo.consensus] == [2, 2, 1]


def test_require_transition_rejects_single_stretch_overlap():
    # two reads sharing only ONE uniform stretch (a long shared satellite) -> no junction -> no edge
    one_stretch = {
        "a": [("chr1:p_arm", 5000), ("chr1:bSat", 8000)],
        "b": [("chr1:bSat", 8000), ("chr2:p_arm", 5000)],
    }
    # their overlap is just bSat (a single stretch type); without require_transition it's an edge...
    assert asm.build_overlap_graph(one_stretch, **PARAMS)
    # ...but with require_transition it is rejected (overlap crosses no junction)
    assert not asm.build_overlap_graph(one_stretch, require_transition=True, **PARAMS)


def test_require_transition_keeps_overlap_spanning_a_junction():
    # reads sharing a junction (bSat -> ITS, a transition between stretch types) -> edge survives
    spanning = {
        "a": [("chr1:p_arm", 5000), ("chr1:bSat", 8000), ("chr1:ITS", 2000)],
        "b": [("chr1:bSat", 8000), ("chr1:ITS", 2000), ("chr1:q_arm", 5000)],
    }
    assert asm.build_overlap_graph(spanning, require_transition=True, **PARAMS)
