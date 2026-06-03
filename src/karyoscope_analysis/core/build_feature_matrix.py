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
from collections.abc import Mapping

from karyoscope_analysis.core import seq_features as sf
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import Interval

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
    """Compute the wide feature matrix.

    Args:
        beds: ``{featureset: {seq_id: [intervals]}}`` (each already C4-validated).
        hierarchy: for feature validation (C2) + interspersion categories.
        interspersion_featureset: which featureset to compute interspersion over
            (typically an overlay composite). ``None`` omits interspersion columns.

    Raises:
        ValueError: on an out-of-taxonomy feature (C2) or an unknown
            ``interspersion_featureset``.
    """
    # C2: validate features for real hierarchy featuresets. Composite/derived
    # featuresets (e.g. an overlay-annotations output, whose name isn't a hierarchy
    # feature set and whose labels may be ``a:b`` composites) are not validated here.
    known_sets = hierarchy.feature_sets()
    for feature_set, by_seq in beds.items():
        if feature_set not in known_sets:
            continue
        seen = {feat for ivals in by_seq.values() for (_, _, feat) in ivals}
        for feat in sorted(seen):
            hierarchy.require_valid_feature(feat, feature_set)

    if interspersion_featureset is not None and interspersion_featureset not in beds:
        raise ValueError(
            f"interspersion featureset {interspersion_featureset!r} is not among the provided "
            f"featuresets {sorted(beds)}"
        )

    rows: dict[str, dict[str, float]] = defaultdict(dict)
    columns: set[str] = set()
    thresholds: list[tuple[str, str, float]] = []

    def put(seq_id: str, column: str, value: float) -> None:
        rows[seq_id][column] = value
        columns.add(column)

    for feature_set in sorted(beds):
        frac_rows: list[dict[str, float]] = []
        for seq_id, ivals in beds[feature_set].items():
            fracs = sf.feature_fraction(ivals)
            bps = sf.feature_bp(ivals)
            put(seq_id, f"{feature_set}__total_bp", sf.total_bp(ivals))
            for feat, frac in fracs.items():
                put(seq_id, f"{feature_set}__frac__{feat}", frac)
                put(seq_id, f"{feature_set}__bp__{feat}", bps[feat])
            for feat, dens in sf.window_densities(
                ivals, window_size=window_size, gap_tol=gap_tol
            ).items():
                for metric in _DENSITY_METRICS:
                    put(seq_id, f"{feature_set}__{metric}__{feat}", getattr(dens, metric))
                put(seq_id, f"{feature_set}__max_block_bp__{feat}", dens.max_block_bp)
            frac_rows.append(fracs)

        feature_thresholds = sf.adaptive_thresholds(
            frac_rows,
            factor=threshold_factor,
            min_thresh=threshold_min,
            max_thresh=threshold_max,
        )
        thresholds.extend(
            (feature_set, feat, thr) for feat, thr in sorted(feature_thresholds.items())
        )

    if interspersion_featureset is not None:
        for seq_id, ivals in beds[interspersion_featureset].items():
            inter = sf.interspersion(ivals, hierarchy)
            for itype in _INTERSPERSION_TYPES:
                put(seq_id, f"interspersion__{itype}", inter[itype])

    return FeatureMatrix(sorted(columns), dict(rows), thresholds)
