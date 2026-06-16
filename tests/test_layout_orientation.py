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
        members=members, seed=seed, size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}  # all-vs-all within the cluster
    return asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1,
        structureless=frozenset({"arm", "p_arm", "q_arm", "ct"}), weight=WEIGHTS,  # arms/ct not landmarks
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
        members=members, seed=members[0], size=len(members),
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
        members=members, seed=members[0], size=len(members),
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


def test_translocation_chromosomes_distinguishes_clean_from_chimeric():
    """The chromosome backbone is opt-in for a *clean* translocation — a consistent chromosome
    junction (one big chromosome + one ≥1 kb partner, the same pair in most reads). A chimera, a
    shared subtelomere, or a sliver array falls back to the structural backbone."""
    acro = frozenset()
    bp = asm.MAJOR_CHROM_BP
    # clean: every read spans the chr11-chr13 junction, both chromosomes large
    clean = {f"r{i}": [("chr11:q_arm", 9000), ("chr13:active_hor", 8000)] for i in range(4)}
    assert asm._translocation_chromosomes(clean, list(clean), acro, bp) == {"chr11", "chr13"}
    # a real junction with a SMALL partner (one big + one >=1kb) is still a translocation
    small = {f"p{i}": [("chr13:active_hor", 11000), ("chr19:active_hor", 2000)] for i in range(3)}
    assert asm._translocation_chromosomes(small, list(small), acro, bp) == {"chr13", "chr19"}
    # subtelomere cluster: 6 chr16-only reads + 2 with a stray chr2 — the chr16-chr2 junction is in
    # only 2/8 reads (< half), so it is not a translocation -> structural backbone
    sub = {f"s{i}": [("chr16:q_arm", 9000), ("chr16:TAR1", 1500)] for i in range(6)}
    sub["t1"] = [("chr16:q_arm", 9000), ("chr2:q_arm", 8000)]
    sub["t2"] = [("chr16:q_arm", 9000), ("chr2:q_arm", 8000)]
    assert asm._translocation_chromosomes(sub, list(sub), acro, bp) == set()
    # a shared subtelomere (a TAR1 array on chr4 + an acrocentric): adjacent but both sides small
    shared = {f"u{i}": [("chr4:TAR1", 1000), ("chr15:TAR1", 2000)] for i in range(4)}
    assert asm._translocation_chromosomes(shared, list(shared), acro, bp) == set()
    # chr4 arm, then a telomere (chr2), then an acrocentric centromere: the telomere run sits between
    # the two big chromosomes, so they're not consecutive -> a chr4-end-then-acrocentric structure,
    # not a fusion -> structural backbone
    telo = {f"v{i}": [("chr4:q_arm", 7000), ("chr2:canonical_telomere", 700), ("chr14:cenSat", 6000)]
            for i in range(3)}
    assert asm._translocation_chromosomes(telo, list(telo), frozenset({"chr14"}), bp) == set()
    # satellite slivers (a bSat array mis-assigned across chromosomes) never reach JUNCTION_PARTNER_BP
    slivers = {f"x{i}": [("chr18:bSat", 5500), ("chr22:bSat", 400), ("chr4:bSat", 350)] for i in range(3)}
    assert asm._translocation_chromosomes(slivers, list(slivers), acro, bp) == set()
    # a "hub" chimera: every read carries the chr19->chr8 junction (a recurring pair) BUT also threads
    # many other substantial chromosomes (a noisy centromeric-satellite read) -> not a *simple*
    # translocation, so it falls back to the structural backbone (HeLa cluster_12).
    hub = {f"h{i}": [("chr5:bSat", 2000), ("chr6:bSat", 2000), ("chr12:bSat", 2000),
                     ("chr19:active_hor", 4000), ("chr8:q_arm", 20000)] for i in range(4)}
    assert asm._translocation_chromosomes(hub, list(hub), acro, bp) == set()
    # but one off-target recurring chromosome is allowed (a real 3-way fusion whose third junction
    # falls just under the half-of-reads bar): chr11-chr13 in all reads, chr13-chr19 in fewer
    threeway = {f"a{i}": [("chr11:q_arm", 9000), ("chr13:active_hor", 8000), ("chr19:active_hor", 6000)]
                for i in range(3)}
    threeway["a0"] = [("chr11:q_arm", 9000), ("chr13:active_hor", 8000)]  # this one lacks chr19
    assert asm._translocation_chromosomes(threeway, list(threeway), acro, bp) >= {"chr11", "chr13"}


def test_anchors_on_filler_breakpoint_for_a_variable_length_feature():
    """A feature whose only neighbor is filler (bSat → arm) anchors on that boundary, so a
    variable-length feature lines up at its breakpoint instead of drifting on its far edge."""
    reads = {
        "a": [("chr4:bSat", 4000), ("chr4:q_arm", 9000)],
        "b": [("chr4:bSat", 5500), ("chr4:q_arm", 9000)],  # longer bSat
        "c": [("chr4:bSat", 4500), ("chr4:q_arm", 9000)],
    }
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    cluster = asm.Cluster(
        members=members, seed=members[0], size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}
    layout = asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1,
        structureless=frozenset({"arm", "p_arm", "q_arm", "ct"}),
    )
    ends = []  # the bSat -> q_arm boundary (bSat end) should land at one coordinate
    for read in layout.placed:
        last = None
        for _s, e, feat in sorted(read.segments):
            if feat.split(":", 1)[1] == "bSat":
                last = e
        ends.append(last)
    assert max(ends) - min(ends) <= 1, ends


def _orient_only(reads, rank, seed, min_bp=150):
    """Helper: run _orient_reads with a chromosome-stripped landmark fn (telomere collapsed)."""
    def landmark_of(feature):
        s = feature.split(":", 1)[1] if ":" in feature else feature
        if s in {"p_arm", "q_arm", "arm", "ct"}:
            return None
        return "telomere" if s in {"canonical_telomere", "noncanonical_telomere"} else s
    return asm._orient_reads(reads, list(reads), seed, rank, landmark_of, min_bp, SUB, 0.1)


def test_orientation_ignores_rare_distal_tail():
    """A rare distal tail landmark must not decide a read's orientation (HeLa cluster_6).

    With a rank where the tail ``telomere`` outranks ``bSat``, comparing the first vs last of *all* a
    read's landmarks makes a ``bSat … TAR1 … telomere`` tail read look ascending (so it is *not*
    flipped) while a plain ``bSat … TAR1`` core read is flipped — leaving their shared bSat-ITS-TAR1
    backbone reading in opposite directions. Orienting on the conserved core landmarks ignores the
    rare telomere and flips them together."""
    rank = {"TAR1": 0, "ITS": 1, "bSat": 2, "telomere": 3}  # telomere outranks bSat (the failure)
    reads = {
        "core1": [("c:bSat", 1000), ("c:ITS", 1000), ("c:TAR1", 1000)],
        "core2": [("c:bSat", 1000), ("c:ITS", 1000), ("c:TAR1", 1000)],
        "core3": [("c:TAR1", 1000), ("c:ITS", 1000), ("c:bSat", 1000)],  # reverse-orientation core
        "tail1": [("c:bSat", 1000), ("c:ITS", 1000), ("c:TAR1", 1000), ("c:ITS", 1000),
                  ("c:TAR1", 1000), ("c:canonical_telomere", 1000)],  # telomere in <half the reads
        "tail2": [("c:bSat", 1000), ("c:ITS", 1000), ("c:TAR1", 1000), ("c:canonical_telomere", 1000)],
    }
    orient = _orient_only(reads, rank, "core1")
    assert orient["tail1"] == orient["core1"], orient
    assert orient["tail2"] == orient["core1"], orient
    assert orient["core2"] == orient["core1"], orient
    assert orient["core3"] != orient["core1"], orient


def test_orientation_by_conserved_junction_ignores_isolated_repeat():
    """A palindromic read (a feature on both ends) orients by the run holding the cluster's conserved
    junction, not the isolated distal copy (HeLa cluster_2).

    The bulk share a ``TAR1 → ITS`` contact (gap 0). The seed has that contact too *and* a second
    TAR1 169 kb away (an isolated copy with only a telomere beside it). Splitting on the gap puts the
    isolated TAR1 in its own run, so the conserved ``TAR1 → ITS`` run orients the seed like the bulk
    rather than the read's palindromic ``TAR1 … TAR1`` ends flipping it on the distal telomere."""
    rank = {"TAR1": 0, "ITS": 1, "telomere": 2}
    reads = {
        "bulk1": [("c:q_arm", 9000), ("c:TAR1", 1785), ("c:ITS", 297)],
        "bulk2": [("c:q_arm", 9000), ("c:TAR1", 1780), ("c:ITS", 300)],
        "bulk3": [("c:ITS", 297), ("c:TAR1", 1785), ("c:q_arm", 9000)],  # reverse
        # seed: TAR1->ITS (gap 0), then 169 kb of arm, then an isolated TAR1 beside a telomere
        "seed": [("c:q_arm", 80000), ("c:TAR1", 1785), ("c:ITS", 297), ("c:q_arm", 169000),
                 ("c:TAR1", 1652), ("c:noncanonical_telomere", 800)],
    }
    orient = _orient_only(reads, rank, "seed")
    # the seed's TAR1->ITS must read the same direction as the bulk's TAR1->ITS
    assert orient["seed"] == orient["bulk1"] == orient["bulk2"], orient
    assert orient["bulk3"] != orient["bulk1"], orient


def test_telomere_subtypes_collapse_to_one_backbone_landmark():
    """A subtelomeric ``TAR1 → telomere`` array aligns even when reads abut different telomere
    *subtypes* (canonical vs noncanonical). Keeping the subtypes as separate backbone landmarks
    split the one ``TAR1 → telomere`` junction into competing ``TAR1 → canonical`` / ``TAR1 →
    noncanonical`` junctions that anchored reads ~5.7 kb apart (HeLa cluster_1); collapsing telomere
    subtypes to one ``telomere`` token lines them up."""
    hier = FeatureHierarchy.from_tsv(HIERARCHY_TSV)
    filler = hier.filler_features
    structureless = filler - hier.canonical_telomere - hier.noncanonical_telomere  # keep telomeres
    # arm → TAR1 → telomere; reads vary in which telomere subtype directly abuts TAR1 and how much
    # telomere they captured. The TAR1 → telomere boundary must land at one coordinate.
    reads = {
        "a": [("chr1:p_arm", 20000), ("chr1:TAR1", 2200), ("chr1:noncanonical_telomere", 500),
              ("chr1:canonical_telomere", 900)],
        "b": [("chr1:p_arm", 17000), ("chr1:TAR1", 2100), ("chr1:canonical_telomere", 1200)],
        "c": [("chr1:p_arm", 24000), ("chr1:TAR1", 2300), ("chr1:noncanonical_telomere", 1100)],
        "d": [("chr1:p_arm", 15000), ("chr1:TAR1", 2150), ("chr1:noncanonical_telomere", 400),
              ("chr1:canonical_telomere", 700)],
    }
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    cluster = asm.Cluster(
        members=members, seed=members[0], size=len(members), orientation_conflict=False
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}
    layout = asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1,
        filler=filler, structureless=structureless,
    )
    ends = []  # the TAR1 -> telomere boundary (TAR1 end) should land at one coordinate
    for read in layout.placed:
        last = None
        for _s, e, feat in sorted(read.segments):
            if feat.split(":", 1)[1] == "TAR1":
                last = e
        if last is not None:
            ends.append(last)
    assert len(ends) == len(reads)
    assert max(ends) - min(ends) <= 1000, ends


def test_small_distinctive_feature_orients_reads():
    """A distinctive feature too small to *anchor* on (a ~400 bp ITS, below MIN_FEATURE_BLOCK_BP)
    still *orients* its read, so reads carrying one big TAR1 don't flip on a coin-toss (cluster_1).

    The arm is near-zero-weighted (as in the real run), so a best-fit-to-seed flip can't use it; only
    the small ITS, seen at the finer ORIENT_MIN_BLOCK_BP, disambiguates orientation.
    """
    assert asm.ORIENT_MIN_BLOCK_BP < 400 < asm.MIN_FEATURE_BLOCK_BP  # the ITS sits in this band
    weights = {**WEIGHTS, "ITS": 0.40}
    reads = {
        "a": [("chr19:q_arm", 9000), ("chr19:TAR1", 1700), ("chr19:ITS", 400)],
        "b": [("chr19:ITS", 400), ("chr19:TAR1", 1700), ("chr19:q_arm", 9000)],  # exact reverse
        "c": [("chr19:q_arm", 8000), ("chr19:TAR1", 1700), ("chr19:ITS", 400)],
        "d": [("chr19:ITS", 400), ("chr19:TAR1", 1700), ("chr19:q_arm", 7000)],  # exact reverse
    }
    members = tuple(sorted(reads, key=lambda r: -sum(length for _f, length in reads[r])))
    cluster = asm.Cluster(
        members=members, seed=members[0], size=len(members),
        orientation_conflict=False,
    )
    neighbors = {r: [o for o in reads if o != r] for r in reads}
    layout = asm.consensus_layout(
        reads, cluster, neighbors=neighbors, sub_score=SUB, gap_factor=0.1,
        structureless=frozenset({"q_arm", "ct"}), weight=weights,
    )
    sides = set()
    for read in layout.placed:
        its, tar1 = _first_start(read, "ITS"), _first_start(read, "TAR1")
        assert its is not None and tar1 is not None
        sides.add(its < tar1)
    assert len(sides) == 1, "ITS lands on inconsistent sides of TAR1 -> orientation not driven by it"
