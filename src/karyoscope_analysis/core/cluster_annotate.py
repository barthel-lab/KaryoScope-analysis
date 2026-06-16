"""Label each cluster from its consensus structure (hierarchy-derived feature classes).

Re-models the legacy ``auto_label_cluster`` decision tree for Engine B: instead of per-read
telomere-density columns (the Ward-era ``sequence_annotate`` format, gone), the label is read off
the cluster **consensus** — which feature classes it contains and which *end* the telomeres sit on.
The feature classes (telomere / ITS-TAR1 / satellite / chromosome) come from the database
``FeatureHierarchy`` (no hardcoded feature names); only the numeric thresholds are parameters,
defaulting to a documented human/CHM13 preset.

Labels (decision order, matching the legacy biology):

1. **ECTR** — telomere at *both* ends (extrachromosomal telomeric repeat / t-circle).
2. **subtelomere** — telomere at *one* end; **Type II ALT subtelomere** if it carries a contiguous
   canonical-telomere block of at least ``alt_block_bp``.
3. **interstitial telomere** — telomere only *internally* (not at an end).
4. **interstitial ITS/TAR1** — internal interstitial-telomeric repeats.
5. **satellite-dominant** — satellite features cover at least ``satellite_fraction`` of the span.
6. **(unlabeled)** — none of the above.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.representatives import consensus_signature

Segment = tuple[int, int, str]  # (start, end, feature) — feature is composite "chrom:structural"


@dataclass(frozen=True)
class LabelConfig:
    """Numeric thresholds for consensus labeling (defaults = human/CHM13 preset)."""

    end_fraction: float = 0.15  # fraction of the consensus span counted as an "end"
    satellite_fraction: float = 0.8  # satellite bp / span to call a cluster satellite-dominant
    alt_block_bp: int = 6000  # contiguous canonical-telomere bp for a Type II ALT subtelomere


@dataclass(frozen=True)
class ClusterAnnotation:
    """One cluster's structural label + catalog metadata."""

    cluster_id: str
    size: int
    width: int
    label: str
    chromosomes: str  # specific chromosomes the consensus spans, e.g. "chr2+chr13" ('' if none)
    signature: str


def _structural(feature: str) -> str:
    """Structural layer of a composite ``chrom:structural`` label (``chr20:gSat`` -> ``gSat``)."""
    return feature.split(":", 1)[1] if ":" in feature else feature


def _chrom(feature: str) -> str:
    """Chromosome layer of a composite label (``chr20:gSat`` -> ``chr20``; ``gSat`` -> ``''``)."""
    return feature.split(":", 1)[0] if ":" in feature else ""


def chromosomes_of(segments: Sequence[Segment], hierarchy: FeatureHierarchy) -> list[str]:
    """Specific chromosomes (from ``hierarchy.chromosomes``) the consensus spans, in natural order."""
    chroms = hierarchy.chromosomes
    present = {c for s in segments if (c := _chrom(s[2])) in chroms}

    def _key(c: str):
        suffix = c[3:]
        return (0, int(suffix)) if suffix.isdigit() else (1, suffix)

    return sorted(present, key=_key)


def label_cluster(
    segments: Sequence[Segment], hierarchy: FeatureHierarchy, cfg: LabelConfig
) -> str:
    """Structural label for a cluster's consensus segments (see module docstring)."""
    if not segments:
        return ""
    telomere = hierarchy.telomere_features
    satellites = hierarchy.satellite_features
    its_tar1 = hierarchy.its_tar1
    canonical = hierarchy.canonical_telomere

    ordered = sorted(segments, key=lambda s: s[0])
    span_start = ordered[0][0]
    span_end = max(e for _s, e, _f in ordered)
    span = max(1, span_end - span_start)
    end_len = cfg.end_fraction * span
    start_edge = span_start + end_len
    end_edge = span_end - end_len

    tel_start = tel_end = tel_internal = its_present = False
    satellite_bp = 0
    canonical_block = run = 0
    for s, e, feature in ordered:
        struct = _structural(feature)
        if struct in telomere:
            if s < start_edge:
                tel_start = True
            elif e > end_edge:
                tel_end = True
            else:
                tel_internal = True
        if struct in its_tar1:
            its_present = True
        if struct in satellites:
            satellite_bp += e - s
        if struct in canonical:
            run += e - s
            canonical_block = max(canonical_block, run)
        else:
            run = 0

    if tel_start and tel_end:
        return "ECTR"
    if tel_start or tel_end:
        return "Type II ALT subtelomere" if canonical_block >= cfg.alt_block_bp else "subtelomere"
    if tel_internal:
        return "interstitial telomere"
    if its_present:
        return "interstitial ITS/TAR1"
    if satellite_bp / span >= cfg.satellite_fraction:
        return "satellite-dominant"
    return ""


def annotate(
    cluster_sizes: Mapping[str, int],
    cluster_widths: Mapping[str, int],
    consensus_by_cluster: Mapping[str, Sequence[Segment]],
    hierarchy: FeatureHierarchy,
    *,
    cfg: LabelConfig | None = None,
    min_cluster_size: int = 2,
) -> list[ClusterAnnotation]:
    """Label every cluster (>= ``min_cluster_size`` reads); sorted by descending size."""
    cfg = cfg or LabelConfig()
    out: list[ClusterAnnotation] = []
    for cluster_id, size in cluster_sizes.items():
        if size < min_cluster_size:
            continue
        segs = list(consensus_by_cluster.get(cluster_id, []))
        out.append(
            ClusterAnnotation(
                cluster_id=cluster_id,
                size=size,
                width=cluster_widths.get(cluster_id, 0),
                label=label_cluster(segs, hierarchy, cfg),
                chromosomes="+".join(chromosomes_of(segs, hierarchy)),
                signature=consensus_signature(segs),
            )
        )
    out.sort(key=lambda a: (-a.size, a.cluster_id))
    return out


def annotation_tsv(rows: Sequence[ClusterAnnotation]) -> str:
    """Render cluster annotations as a TSV (one row per cluster)."""
    lines = ["cluster_id\tsize\twidth\tlabel\tchromosomes\tconsensus_signature"]
    for r in rows:
        lines.append(
            f"{r.cluster_id}\t{r.size}\t{r.width}\t{r.label}\t{r.chromosomes}\t{r.signature}"
        )
    return "\n".join(lines) + "\n"
