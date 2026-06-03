"""overlay-annotations: combine per-featureset annotation tracks into one resolved BED.

Given one annotation BED per featureset (already C4-validated by the reader), for
each ``seq_id`` this refines the tracks at the union of their boundaries
(:func:`karyoscope_analysis.core.intervals.refine`), resolves each segment's
feature tuple to a single feature via a :class:`ResolutionSpec`, and coalesces
adjacent same-feature segments. The featuresets are exactly the spec's
``precedence`` (and define the refinement/composite order).
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping

from karyoscope_analysis.core import intervals
from karyoscope_analysis.core.annotation_resolution import ResolutionSpec
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, Interval


def validate_features(
    beds: Mapping[str, Mapping[str, list[Interval]]], hierarchy: FeatureHierarchy
) -> None:
    """Check every feature value in each track is valid for its featureset (C2)."""
    for feature_set, by_seq in beds.items():
        seen = {feat for ivals in by_seq.values() for (_, _, feat) in ivals}
        for feat in sorted(seen):
            hierarchy.require_valid_feature(feat, feature_set)


def overlay_annotations(
    beds: Mapping[str, Mapping[str, list[Interval]]],
    spec: ResolutionSpec,
    hierarchy: FeatureHierarchy | None = None,
) -> Iterator[BedRow]:
    """Yield resolved ``(seq_id, start, end, feature)`` rows.

    Args:
        beds: ``{featureset: {seq_id: [(start, end, feature), ...]}}``; keys must be
            exactly ``spec.precedence``.
        spec: the resolution spec (precedence + rules).
        hierarchy: if given, every input feature is validated (C2) before resolving.

    Raises:
        ValueError: if the BED featuresets don't match ``spec.precedence``, the
            featuresets cover different ``seq_id`` sets, or (via ``refine``) a
            sequence's tracks don't share a gapless common span.
    """
    featuresets = spec.precedence
    if set(beds) != set(featuresets):
        raise ValueError(
            f"BED featuresets {sorted(beds)} must match the spec precedence {list(featuresets)}"
        )
    if hierarchy is not None:
        validate_features(beds, hierarchy)

    reference = featuresets[0]
    ref_seqs = set(beds[reference])
    for fs in featuresets[1:]:
        if set(beds[fs]) != ref_seqs:
            only_ref = sorted(ref_seqs - set(beds[fs]))[:3]
            only_fs = sorted(set(beds[fs]) - ref_seqs)[:3]
            raise ValueError(
                f"featureset {fs!r} covers a different seq_id set than {reference!r} "
                f"(only in {reference}: {only_ref}; only in {fs}: {only_fs})"
            )

    for seq_id in beds[reference]:
        tracks = [beds[fs][seq_id] for fs in featuresets]
        try:
            segments = intervals.refine(tracks)
        except ValueError as exc:
            raise ValueError(f"seq_id {seq_id!r}: {exc}") from exc
        resolved = [
            (start, end, spec.resolve(dict(zip(featuresets, values, strict=True))))
            for start, end, values in segments
        ]
        for start, end, feature in intervals.coalesce(resolved):
            yield seq_id, start, end, feature
