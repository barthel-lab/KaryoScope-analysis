"""overlay-annotations: combine per-featureset annotation tracks into one resolved BED.

This is a **single-pass, streaming k-way overlay**. Given one annotation BED per
featureset — each a C4 partition of the same per-``seq_id`` span, with sequences in the
same order — it reads all the tracks concurrently and, for each ``seq_id``, sweeps a line
across the union of the tracks' boundaries. At every segment it has exactly one feature
per track; it resolves that tuple to a single feature via a :class:`ResolutionSpec` and
coalesces adjacent same-feature segments. Only the *current* interval of each track is
held, so inputs of any size overlay in ``O(featuresets)`` memory and exactly one pass each
(the goal stated for this tool; cf. the legacy in-memory ``merge_beds``).

Two lessons from the earlier hand-written streaming merges are baked in:

* **Don't order ``seq_id``s by name.** Lexicographic comparison misorders natural names
  (``chr1`` < ``chr10`` < ``chr2``). Instead every track is advanced in lockstep and
  required to present sequences in the *same* order; a disagreement is an error, so there
  is no name comparison to get wrong.
* **A true k-way sweep, not chained pairwise merges.** Resolution needs every featureset's
  value at a position simultaneously, which one k-way sweep gives directly.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping

from karyoscope_analysis.core.annotation_resolution import ResolutionSpec
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, Interval


class _Cursor:
    """A one-row-lookahead view over a single track's annotation-row stream."""

    __slots__ = ("_rows", "current")

    def __init__(self, rows: Iterator[BedRow]) -> None:
        self._rows = rows
        self.current: BedRow | None = next(rows, None)

    def advance(self) -> None:
        self.current = next(self._rows, None)


def overlay_streams(
    named_streams: Mapping[str, Iterator[BedRow]],
    spec: ResolutionSpec,
    hierarchy: FeatureHierarchy | None = None,
) -> Iterator[BedRow]:
    """Overlay several annotation-row streams in one pass, yielding resolved rows.

    Args:
        named_streams: ``{featureset: iterator of (seq_id, start, end, feature)}``; keys
            must be exactly ``spec.precedence``. Each stream must be a C4 partition per
            sequence, and all streams must present sequences in the same order (e.g. the
            streaming reader :func:`karyoscope_analysis.core.io.bed.iter_annotation_rows`).
        spec: the resolution spec (precedence + rules).
        hierarchy: if given, every input feature is validated (C2) as it is first seen.

    Yields:
        Resolved ``(seq_id, start, end, feature)`` rows, coalesced within each sequence.

    Raises:
        ValueError: if the featuresets don't match ``spec.precedence``, the streams
            disagree on sequence order / coverage, or a sequence's tracks don't share a
            common span.
    """
    featuresets = spec.precedence
    if set(named_streams) != set(featuresets):
        raise ValueError(
            f"BED featuresets {sorted(named_streams)} must match the spec precedence "
            f"{list(featuresets)}"
        )

    cursors = [_Cursor(named_streams[fs]) for fs in featuresets]
    validated: set[tuple[str, str]] = set()

    def resolve_here(seq_id: str) -> str:
        """Resolve the feature tuple covering the current sweep position."""
        values: list[str] = []
        for fs, cur in zip(featuresets, cursors, strict=True):
            feat = cur.current[3]  # type: ignore[index]  # current is non-None here
            if hierarchy is not None and (fs, feat) not in validated:
                hierarchy.require_valid_feature(feat, fs)
                validated.add((fs, feat))
            values.append(feat)
        return spec.resolve(dict(zip(featuresets, values, strict=True)))

    while True:
        if all(c.current is None for c in cursors):
            return

        # A new sequence begins: every track must be live and on the same seq_id.
        seq_id: str | None = None
        for fs, cur in zip(featuresets, cursors, strict=True):
            if cur.current is None:
                raise ValueError(
                    f"featureset {fs!r} ran out of sequences before the others; "
                    f"featuresets must cover the same seq_ids in the same order"
                )
            if seq_id is None:
                seq_id = cur.current[0]
            elif cur.current[0] != seq_id:
                raise ValueError(
                    f"featureset {fs!r} is at seq_id {cur.current[0]!r} but expected "
                    f"{seq_id!r}; featuresets must list seq_ids in the same order"
                )
        assert seq_id is not None

        lo = cursors[0].current[1]  # type: ignore[index]
        for fs, cur in zip(featuresets, cursors, strict=True):
            if cur.current[1] != lo:  # type: ignore[index]
                raise ValueError(
                    f"seq_id {seq_id!r}: featureset {fs!r} starts at {cur.current[1]} "  # type: ignore[index]
                    f"but {featuresets[0]!r} starts at {lo} (tracks must share a span)"
                )

        # Sweep this sequence, coalescing the resolved run as we go.
        pos = lo
        run_start = run_end = pos
        run_feat: str | None = None
        while True:
            boundary = min(c.current[2] for c in cursors)  # type: ignore[index]
            feat = resolve_here(seq_id)
            if run_feat is not None and feat == run_feat and run_end == pos:
                run_end = boundary
            else:
                if run_feat is not None:
                    yield seq_id, run_start, run_end, run_feat
                run_start, run_end, run_feat = pos, boundary, feat
            pos = boundary

            for cur in cursors:
                if cur.current[2] == boundary:  # type: ignore[index]
                    cur.advance()

            on_seq = [c.current is not None and c.current[0] == seq_id for c in cursors]
            if all(on_seq):
                continue  # more of this sequence remains
            if any(on_seq):
                raise ValueError(
                    f"seq_id {seq_id!r}: featuresets cover different spans "
                    f"(some tracks end before others)"
                )
            # Sequence complete: flush its coalesced run and move on.
            if run_feat is not None:
                yield seq_id, run_start, run_end, run_feat
            break


def overlay_annotations(
    beds: Mapping[str, Mapping[str, list[Interval]]],
    spec: ResolutionSpec,
    hierarchy: FeatureHierarchy | None = None,
) -> Iterator[BedRow]:
    """Overlay in-memory tracks (convenience wrapper over :func:`overlay_streams`).

    Accepts already-read ``{featureset: {seq_id: [(start, end, feature), ...]}}`` and
    streams the rows through the same single-pass engine. Unlike the file path it is
    order-tolerant: all featuresets are emitted in the reference featureset's seq_id order
    (and each sequence's intervals are sorted), so dicts built in any order work.

    Raises:
        ValueError: if the featuresets don't match ``spec.precedence``, cover different
            ``seq_id`` sets, or a sequence's tracks don't share a common span.
    """
    featuresets = spec.precedence
    if set(beds) != set(featuresets):
        raise ValueError(
            f"BED featuresets {sorted(beds)} must match the spec precedence {list(featuresets)}"
        )

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

    def stream_for(fs: str) -> Iterator[BedRow]:
        per_seq = beds[fs]
        for seq_id in beds[reference]:  # reference order, for all featuresets
            for start, end, feature in sorted(per_seq[seq_id]):
                yield seq_id, start, end, feature

    yield from overlay_streams({fs: stream_for(fs) for fs in featuresets}, spec, hierarchy)
