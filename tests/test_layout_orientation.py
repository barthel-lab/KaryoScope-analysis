"""Regression test for progressive-layout orientation on a real 3-chromosome cluster.

Captured from full-sample U2OS clustering (``cluster_6``): a chr11 - chr13 - chr19 structure
where chr13 is a large, ~uniform ``active_hor`` hub. Reads cover either the chr11-chr13 junction
or the chr13-chr19 junction (none span all three). The chr11 reads come in both orientations
(``[q_arm, TAR1, chr13:active_hor]`` and its exact reverse). The layout must orient every read so
the chromosomes read in one consistent left-to-right order and the shared junctions line up,
*without* a uniform hub flipping reads. See ``docs/audit/rearrangement_detection.md`` (Engine B).
"""

from __future__ import annotations

from itertools import pairwise
from pathlib import Path

from karyoscope_analysis.core import feature_assembly as asm
from karyoscope_analysis.core.feature_align import chromosome_aware_substitution
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"

# Exact reads from cluster_6 (coalesced (composite_feature, length) segments).
CLUSTER6: dict[str, list[tuple[str, int]]] = {
    "b1cdc7e9": [("chr13:active_hor", 6857), ("chr19:active_hor", 92612)],
    "ce8fc11a": [
        ("acrocentric:noncanonical_telomere", 120), ("chr13:noncanonical_telomere", 43),
        ("chr13:active_hor", 2233), ("chr13:dhor", 910), ("chr13:active_hor", 14083),
        ("acrocentric:active_hor", 527), ("chr13:active_hor", 6532),
    ],
    "fa83bad4": [
        ("chr13:active_hor", 5906), ("chr13:dhor", 7), ("chr13:active_hor", 15236),
        ("chr13:dhor", 527), ("chr13:active_hor", 20699), ("chr19:active_hor", 3662),
    ],
    "r144772": [
        ("chr11:q_arm", 2274), ("chr11:ct", 442), ("chr11:q_arm", 9786), ("chr11:gSat", 48),
        ("chr11:TAR1", 1177), ("chr13:TAR1", 195), ("chr13:active_hor", 6220),
    ],
    "r221775": [
        ("chr13:active_hor", 8449), ("chr13:TAR1", 6), ("chr11:TAR1", 1338),
        ("chr11:gSat", 48), ("chr11:q_arm", 1800),
    ],
    "r166594": [
        ("chr13:active_hor", 7042), ("chr13:TAR1", 6), ("chr11:TAR1", 1343),
        ("chr11:gSat", 48), ("chr11:q_arm", 7681),
    ],
    "r212371": [
        ("chr13:alpha_hor", 1), ("chr13:active_hor", 9834), ("chr13:alpha_hor", 19),
        ("chr13:active_hor", 5591), ("chr13:TAR1", 6), ("chr11:TAR1", 1340),
        ("chr11:gSat", 48), ("chr11:q_arm", 9808),
    ],
    "r218256": [("chr19:active_hor", 9763), ("chr13:active_hor", 6845)],
    "r220268": [("chr19:active_hor", 5601), ("chr13:active_hor", 7232)],
    "r234684": [
        ("chr13:dhor", 36), ("chr13:active_hor", 7988), ("chr13:TAR1", 72),
        ("chr11:TAR1", 1265), ("chr11:gSat", 48), ("chr11:q_arm", 8008),
    ],
    "r963850": [
        ("chr11:q_arm", 7638), ("chr11:gSat", 48), ("chr11:TAR1", 1338),
        ("chr13:TAR1", 6), ("chr13:active_hor", 6185),
    ],
}

# Exact-or-mismatch structural scorer, made chromosome-aware (cross-chromosome penalty), matching
# how the cluster command scores composite labels.
_STRUCT = lambda a, b: 1.0 if a == b else -1.0  # noqa: E731
SUB = chromosome_aware_substitution(_STRUCT, cross_chromosome_penalty=-2.0)

# Genome-frequency-style weights (as in the real run): the chr11 arm is down-weighted to near
# zero, so the big chr13 ``active_hor`` hub dominates the *score* — the flank that disambiguates
# orientation contributes little, which is exactly what made the layout flip reads.
WEIGHTS = {
    "q_arm": 0.06, "p_arm": 0.06, "ct": 0.12, "active_hor": 0.42, "dhor": 0.42,
    "alpha_hor": 0.42, "TAR1": 0.32, "gSat": 0.28, "noncanonical_telomere": 0.12,
}


def _chrom_order(segments: tuple) -> list[str]:
    """The specific chromosomes a laid-out read spans, left to right (consecutive dups collapsed)."""
    order: list[str] = []
    for _s, _e, feat in sorted(segments):
        chrom = feat.split(":", 1)[0]
        if chrom.startswith("chr") and (not order or order[-1] != chrom):
            order.append(chrom)
    return order


def _layout_cluster6() -> asm.ClusterLayout:
    reads = {rid: list(segs) for rid, segs in CLUSTER6.items()}
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    seed = members[0]  # longest read (chr13-chr19)
    cluster = asm.Cluster(
        members=members, seed=seed, reversed_relative_to_seed={}, size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}  # all-vs-all within the cluster
    return asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1, weight=WEIGHTS
    )


def test_cluster6_chromosome_order_is_consistent():
    """Every read must agree on the relative order of chromosomes (no inversions)."""
    layout = _layout_cluster6()
    left_of: set[tuple[str, str]] = set()
    for read in layout.placed:
        order = _chrom_order(read.segments)
        for i, x in enumerate(order):
            for y in order[i + 1 :]:
                left_of.add((x, y))
    conflicts = {(x, y) for (x, y) in left_of if (y, x) in left_of}
    assert not conflicts, f"chromosomes ordered inconsistently across reads: {sorted(conflicts)}"


def test_cluster6_shared_junctions_align():
    """The chr11-chr13 and chr13-chr19 junctions land at a consistent consensus coordinate."""
    layout = _layout_cluster6()
    pair_boundary: dict[frozenset, list[float]] = {}
    for read in layout.placed:
        blocks: list[tuple[str, int, int]] = []
        for s, e, feat in sorted(read.segments):
            chrom = feat.split(":", 1)[0]
            if not chrom.startswith("chr"):
                continue
            if blocks and blocks[-1][0] == chrom:
                blocks[-1] = (chrom, blocks[-1][1], e)
            else:
                blocks.append((chrom, s, e))
        for (ca, _sa, ea), (cb, sb, _eb) in pairwise(blocks):
            pair_boundary.setdefault(frozenset({ca, cb}), []).append((ea + sb) / 2)
    for pair, coords in pair_boundary.items():
        if len(coords) >= 2:
            spread = max(coords) - min(coords)
            assert spread <= 3000, f"{set(pair)} junction varies by {spread:.0f}bp: {coords}"


# Derived from full-sample cluster_2 (a single chromosome, chr18): a bSat - TAR1 distinctive
# backbone (ITS is a small feature between them, below the landmark threshold), captured in both
# orientations with varying arm. The backbone here is the *structural* features, not chromosomes.
CHR18 = {
    "a_fwd": [
        ("chr18:p_arm", 30000), ("chr18:bSat", 7000), ("chr18:ITS", 600),
        ("chr18:p_arm", 5000), ("chr18:TAR1", 2000), ("chr18:p_arm", 20000),
    ],
    "b_fwd": [
        ("chr18:p_arm", 18000), ("chr18:bSat", 7400), ("chr18:ITS", 590),
        ("chr18:p_arm", 5000), ("chr18:TAR1", 1990), ("chr18:p_arm", 32000),
    ],
    "c_rev": [  # same molecule, reversed: TAR1 ... bSat
        ("chr18:p_arm", 20000), ("chr18:TAR1", 2000), ("chr18:p_arm", 5000),
        ("chr18:ITS", 588), ("chr18:bSat", 7400),
    ],
    "d_rev": [
        ("chr18:p_arm", 19000), ("chr18:TAR1", 1880), ("chr18:p_arm", 4800),
        ("chr18:ITS", 593), ("chr18:bSat", 3860),
    ],
    # reads lacking TAR1 (the user's misaligned ones): only bSat + ITS. They share the bSat->ITS
    # junction with the others and must align there, not drift because they have no TAR1.
    "e_noTAR1": [("chr18:p_arm", 25000), ("chr18:bSat", 7200), ("chr18:ITS", 595), ("chr18:p_arm", 3000)],
    "f_noTAR1": [("chr18:p_arm", 22000), ("chr18:bSat", 3000), ("chr18:ITS", 590), ("chr18:p_arm", 8000)],
    # a read lacking the *proximal* bSat (only ITS + TAR1): no adjacent bSat->ITS junction, so it
    # must still pin its ITS to the ITS slot rather than drifting via the distal TAR1.
    "g_noBSAT": [("chr18:p_arm", 24000), ("chr18:ITS", 588), ("chr18:p_arm", 5000), ("chr18:TAR1", 1990)],
}


def _layout_chr18() -> asm.ClusterLayout:
    filler = FeatureHierarchy.from_tsv(HIERARCHY_TSV).filler_features
    reads = {rid: list(segs) for rid, segs in CHR18.items()}
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    cluster = asm.Cluster(
        members=members, seed=members[0], reversed_relative_to_seed={}, size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}
    return asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1, filler=filler
    )


def _first_start(read, structural):
    """Consensus start of a read's first segment whose structural layer is ``structural``."""
    for s, _e, feat in sorted(read.segments):
        if feat.split(":", 1)[1] == structural:
            return s
    return None


def test_single_chromosome_uses_structural_backbone():
    """A one-chromosome cluster is oriented + ordered by its distinctive features (bSat - TAR1).

    The reverse-orientation reads must be flipped so every read reads bSat before TAR1.
    """
    layout = _layout_chr18()
    # every read carrying both must place bSat to the left of TAR1 (consistent order, no inversions)
    for read in layout.placed:
        b, t = _first_start(read, "bSat"), _first_start(read, "TAR1")
        if b is not None and t is not None:
            assert b < t, f"{read.read_id}: bSat={b} TAR1={t}"
    # the two reverse-orientation reads were flipped
    flipped = {r.read_id for r in layout.placed if r.reversed}
    assert flipped == {"c_rev", "d_rev"}, flipped


def test_single_chromosome_its_aligns_without_the_distal_feature():
    """ITS lines up across reads even when some lack the distal TAR1 feature.

    The reads sharing only bSat + ITS (no TAR1) must anchor on the bSat->ITS junction like the
    rest, not drift to a different landmark — the bug the maintainer spotted on cluster_2.
    """
    layout = _layout_chr18()
    its_starts = [s for r in layout.placed if (s := _first_start(r, "ITS")) is not None]
    assert len(its_starts) == len(CHR18)  # every read has ITS placed
    assert max(its_starts) - min(its_starts) <= 1500, its_starts


def test_reads_sorted_by_consensus_start():
    """Placed reads are ordered top-to-bottom by where they start in the consensus."""
    for layout in (_layout_cluster6(), _layout_chr18()):
        starts = [min(s for s, _e, _f in r.segments) for r in layout.placed]
        assert starts == sorted(starts), starts


def test_acrocentric_chromosomes_collapse_in_backbone():
    """Reads sharing chr4 but assigned to *different* acrocentrics (chr15 vs chr21) lay out as one
    chr4 - acrocentric structure, so the recombining acrocentric short arms don't tangle the backbone."""
    acro = frozenset({"chr13", "chr14", "chr15", "chr21", "chr22"})
    reads = {
        "a": [("chr4:p_arm", 20000), ("chr4:bSat", 6000), ("chr15:TAR1", 2000)],
        "b": [("chr4:p_arm", 18000), ("chr4:bSat", 6000), ("chr21:TAR1", 2000)],  # chr21, not chr15
        "c": [("chr4:p_arm", 22000), ("chr4:bSat", 6000), ("chr22:TAR1", 2000)],
    }
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    cluster = asm.Cluster(
        members=members, seed=members[0], reversed_relative_to_seed={}, size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}
    layout = asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1,
        acrocentric_chromosomes=acro,
    )
    # the chr4->acrocentric junction (chr4:bSat end / acrocentric:TAR1 start) lands at one coordinate
    junctions = []
    for read in layout.placed:
        for s, _e, feat in sorted(read.segments):
            if feat.split(":", 1)[1] == "TAR1":
                junctions.append(s)
                break
    assert max(junctions) - min(junctions) <= 500, junctions
