"""Assemble the wide per-``seq_id`` feature matrix from per-featureset annotation BEDs.

Produces the matrix consumed by ``cluster-annotate``/``cluster-diagnostics`` and an
adaptive-threshold sidecar (F5). Column schema (decision F2, ``__`` is the sole
delimiter):

* ``{featureset}__{metric}__{feature}`` for ``frac``, ``bp``, the density metrics
  (``dmax``/``dmin``/``dmedian``/``dfirst``/``dlast``/``dterminal``/``dterminal_min``),
  and ``max_block_bp``;
* ``{featureset}__total_bp``;
* ``interspersion__{type}`` (only when an interspersion featureset is chosen);
* key column ``seq_id``.

Alignment-QC columns are intentionally *not* here — they move to ``cluster-diagnostics``
(decision F6). Input features are validated against the hierarchy (C2).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping

from karyoscope_analysis.core import seq_features as sf
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, Interval

#: Density metric columns, in a stable order.
_DENSITY_METRICS = ("dmax", "dmin", "dmedian", "dfirst", "dlast", "dterminal", "dterminal_min")
_INTERSPERSION_TYPES = ("total", "can_ncan", "tel_sat", "arm_tel")


class FeatureMatrix:
    """A wide per-``seq_id`` feature matrix plus the adaptive-threshold sidecar."""

    def __init__(
        self,
        columns: list[str],
        rows: Mapping[str, Mapping[str, float]],
        thresholds: list[tuple[str, str, float]],
    ) -> None:
        self.columns = columns  # sorted feature columns (excludes the seq_id key)
        self.rows = rows  # seq_id -> {column: value} (sparse; missing => 0)
        self.thresholds = thresholds  # (featureset, feature, threshold), sorted

    def seq_ids(self) -> list[str]:
        return sorted(self.rows)


#: One per-``seq_id`` group: ``{featureset: [(start, end, feature), ...]}``.
SequenceGroup = tuple[str, Mapping[str, list[Interval]]]


def build_matrix_from_groups(
    groups: Iterable[SequenceGroup],
    featureset_names: set[str],
    hierarchy: FeatureHierarchy,
    *,
    window_size: int = sf.DEFAULT_WINDOW_SIZE,
    gap_tol: int = sf.DEFAULT_BLOCK_GAP_TOL,
    threshold_factor: float = sf.DEFAULT_THRESHOLD_FACTOR,
    threshold_min: float = sf.DEFAULT_THRESHOLD_MIN,
    threshold_max: float = sf.DEFAULT_THRESHOLD_MAX,
    interspersion_featureset: str | None = None,
) -> FeatureMatrix:
    """Build the matrix from a stream of per-``seq_id`` groups (the shared engine).

    Each group carries one sequence's intervals for every featureset, so only one
    sequence is in memory at a time. Both the in-memory :func:`build_feature_matrix`
    and the file-streaming :func:`build_feature_matrix_streaming` feed this.

    Args:
        groups: iterable of ``(seq_id, {featureset: [intervals]})``.
        featureset_names: all featuresets present (for the interspersion-name check).
        interspersion_featureset: featureset to compute interspersion over, or ``None``.

    Raises:
        ValueError: on an out-of-taxonomy feature (C2) or an unknown
            ``interspersion_featureset``.
    """
    if interspersion_featureset is not None and interspersion_featureset not in featureset_names:
        raise ValueError(
            f"interspersion featureset {interspersion_featureset!r} is not among the provided "
            f"featuresets {sorted(featureset_names)}"
        )

    known_sets = hierarchy.feature_sets()
    rows: dict[str, dict[str, float]] = defaultdict(dict)
    columns: set[str] = set()
    frac_accum: dict[str, list[dict[str, float]]] = defaultdict(list)
    validated: set[tuple[str, str]] = set()

    def put(seq_id: str, column: str, value: float) -> None:
        rows[seq_id][column] = value
        columns.add(column)

    for seq_id, by_fs in groups:
        for feature_set, ivals in by_fs.items():
            # C2: validate features of real hierarchy featuresets (skip composite/derived
            # featuresets like an overlay output, whose labels may be ``a:b``). Validated
            # once per (featureset, feature) and cached, before any column is emitted.
            if feature_set in known_sets:
                for _, _, feat in ivals:
                    key = (feature_set, feat)
                    if key not in validated:
                        hierarchy.require_valid_feature(feat, feature_set)
                        validated.add(key)

            # One pass for per-feature bp; total_bp and fractions derive from it
            # (no recomputation), matching sf.feature_fraction's empty-on-zero behaviour.
            bps = sf.feature_bp(ivals)
            total = sum(bps.values())
            put(seq_id, f"{feature_set}__total_bp", total)
            fracs = {feat: bp / total for feat, bp in bps.items()} if total > 0 else {}
            for feat, frac in fracs.items():
                put(seq_id, f"{feature_set}__frac__{feat}", frac)
                put(seq_id, f"{feature_set}__bp__{feat}", bps[feat])
            for feat, dens in sf.window_densities(
                ivals, window_size=window_size, gap_tol=gap_tol
            ).items():
                for metric in _DENSITY_METRICS:
                    put(seq_id, f"{feature_set}__{metric}__{feat}", getattr(dens, metric))
                put(seq_id, f"{feature_set}__max_block_bp__{feat}", dens.max_block_bp)
            frac_accum[feature_set].append(fracs)

        if interspersion_featureset is not None and interspersion_featureset in by_fs:
            inter = sf.interspersion(by_fs[interspersion_featureset], hierarchy)
            for itype in _INTERSPERSION_TYPES:
                put(seq_id, f"interspersion__{itype}", inter[itype])

    thresholds: list[tuple[str, str, float]] = []
    for feature_set in sorted(frac_accum):
        feature_thresholds = sf.adaptive_thresholds(
            frac_accum[feature_set],
            factor=threshold_factor,
            min_thresh=threshold_min,
            max_thresh=threshold_max,
        )
        thresholds.extend(
            (feature_set, feat, thr) for feat, thr in sorted(feature_thresholds.items())
        )

    return FeatureMatrix(sorted(columns), dict(rows), thresholds)


def build_feature_matrix(
    beds: Mapping[str, Mapping[str, list[Interval]]],
    hierarchy: FeatureHierarchy,
    *,
    window_size: int = sf.DEFAULT_WINDOW_SIZE,
    gap_tol: int = sf.DEFAULT_BLOCK_GAP_TOL,
    threshold_factor: float = sf.DEFAULT_THRESHOLD_FACTOR,
    threshold_min: float = sf.DEFAULT_THRESHOLD_MIN,
    threshold_max: float = sf.DEFAULT_THRESHOLD_MAX,
    interspersion_featureset: str | None = None,
) -> FeatureMatrix:
    """Compute the wide feature matrix from in-memory tracks (convenience wrapper).

    Groups the already-read ``{featureset: {seq_id: [intervals]}}`` by ``seq_id`` and
    runs them through :func:`build_matrix_from_groups`. Order-tolerant and supports
    featuresets that cover different ``seq_id`` sets (a sequence simply contributes
    whichever featuresets contain it). For files, prefer
    :func:`build_feature_matrix_streaming`.

    Raises:
        ValueError: on an out-of-taxonomy feature (C2) or an unknown
            ``interspersion_featureset``.
    """
    ordered_seqs: list[str] = []
    seen: set[str] = set()
    for by_seq in beds.values():
        for seq_id in by_seq:
            if seq_id not in seen:
                seen.add(seq_id)
                ordered_seqs.append(seq_id)

    def groups() -> Iterable[SequenceGroup]:
        for seq_id in ordered_seqs:
            yield seq_id, {fs: beds[fs][seq_id] for fs in beds if seq_id in beds[fs]}

    return build_matrix_from_groups(
        groups(),
        set(beds),
        hierarchy,
        window_size=window_size,
        gap_tol=gap_tol,
        threshold_factor=threshold_factor,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        interspersion_featureset=interspersion_featureset,
    )


def build_feature_matrix_streaming(
    named_streams: Mapping[str, Iterator[BedRow]],
    hierarchy: FeatureHierarchy,
    *,
    window_size: int = sf.DEFAULT_WINDOW_SIZE,
    gap_tol: int = sf.DEFAULT_BLOCK_GAP_TOL,
    threshold_factor: float = sf.DEFAULT_THRESHOLD_FACTOR,
    threshold_min: float = sf.DEFAULT_THRESHOLD_MIN,
    threshold_max: float = sf.DEFAULT_THRESHOLD_MAX,
    interspersion_featureset: str | None = None,
) -> FeatureMatrix:
    """Build the matrix by streaming all featureset BEDs concurrently (one read at a time).

    ``named_streams`` maps each featureset to an annotation-row iterator (e.g.
    :func:`karyoscope_analysis.core.io.bed.iter_annotation_rows`). Streams are walked in
    lockstep via :func:`~karyoscope_analysis.core.io.bed.iter_aligned_groups`, so the input
    is never fully loaded — peak memory is one sequence's intervals plus the output. All
    featuresets must list the same sequences in the same order.

    Raises:
        ValueError: on an out-of-taxonomy feature (C2), an unknown
            ``interspersion_featureset``, or streams that disagree on sequence order/coverage.
    """
    from karyoscope_analysis.core.io.bed import iter_aligned_groups

    return build_matrix_from_groups(
        iter_aligned_groups(named_streams),
        set(named_streams),
        hierarchy,
        window_size=window_size,
        gap_tol=gap_tol,
        threshold_factor=threshold_factor,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        interspersion_featureset=interspersion_featureset,
    )
