"""Shared feature-vocabulary constants for KaryoScope analysis scripts.

The KaryoScope feature databases evolved between versions:

- KS_human_CHM13 (v1) used `_specific` / `_multigroup1` suffixes and lowercase
  satellite names (e.g. `bsat_specific`, `hsat3_specific`, `gsat_specific`).
- KS_human_CHM13_v2 dropped both suffixes and switched satellites to mixed-case
  bare names (`bSat`, `HSat3`, `gSat`, `cenSat`, ...) with finer subdivisions
  (`active_hor`, `alpha_hor`, `mon`, `dhor`, etc.).

To keep the annotation/labeling/representative-selection logic correct against
either database without per-call special-casing, this module exposes:

- The full v2-canonical satellite vocabulary
- A v1→v2 alias table
- Helpers (`is_satellite`, `lookup_satellite_col`) that try v2 names first
  then fall back to v1 aliases for back-compat with existing annotation TSVs.

The sets below are deliberately permissive — they include both v1 and v2 names
so that classification works regardless of which database produced the BEDs.
"""

# v2-canonical satellite features (per KS_human_CHM13_v2.features.txt)
SATELLITE_V2 = frozenset({
    'aSat', 'bSat', 'gSat', 'cenSat', 'centromeric',
    'HSat', 'HSat1', 'HSat1A', 'HSat1B', 'HSat2', 'HSat3',
    'dhor', 'hor', 'active_hor', 'alpha_hor', 'mixedAlpha', 'mon',
})

# v1 satellite aliases (per KS_human_CHM13.features.txt). These appear in
# annotation TSVs produced by older runs; lookups should fall back to them.
SATELLITE_V1 = frozenset({
    'asat_multigroup1', 'bsat', 'gsat', 'censet', 'censat',
    'hsat_multigroup1', 'hsat1_multigroup1',
    'hsat1A', 'hsat1B', 'hsat2', 'hsat3',
    'hor_multigroup1', 'monomeric', 'active', 'inactive', 'divergent',
})

# Union — used wherever the question is "is this any kind of satellite?"
SATELLITE_FEATURES = SATELLITE_V2 | SATELLITE_V1

# v1 → v2 aliases for satellite columns. The semantic mapping for `active` and
# `monomeric` is approximate: v1's `active` ≈ v2's `active_hor` (active alpha
# higher-order repeat), v1's `monomeric` ≈ v2's `mon` (monomeric satellite).
# Callers that need an exact substitution should use this table; otherwise
# `lookup_satellite_col` does best-effort lookup.
SATELLITE_V1_TO_V2 = {
    'bsat': 'bSat',
    'gsat': 'gSat',
    'censat': 'cenSat',
    'hsat1A': 'HSat1A',
    'hsat1B': 'HSat1B',
    'hsat2': 'HSat2',
    'hsat3': 'HSat3',
    'asat_multigroup1': 'aSat',
    'hor_multigroup1': 'hor',
    'hsat_multigroup1': 'HSat',
    'hsat1_multigroup1': 'HSat1',
    'active': 'active_hor',
    'monomeric': 'mon',
}

# Inverse for v2→v1 fallback when reading older TSVs.
SATELLITE_V2_TO_V1 = {v: k for k, v in SATELLITE_V1_TO_V2.items()}

# Chromosome-arm region features. v2 dropped the `_multigroup1` suffix and
# uses bare `arm` alongside `p_arm` / `q_arm`.
ARM_FEATURES = frozenset({'p_arm', 'q_arm', 'arm', 'arm_multigroup1'})

# Centromere-transition / inter-region features
CT_FEATURES = frozenset({'ct', 'ct_specific'})

# Telomere-type features (subtelomeric featureset). These names didn't change
# between v1 (`canonical_telomere_specific`) and v2 (`canonical_telomere`)
# beyond `_specific` stripping, which other code already normalizes.
CANONICAL_TELOMERE = frozenset({'canonical_telomere'})
NONCANONICAL_TELOMERE = frozenset({'noncanonical_telomere'})
ITS_TAR1 = frozenset({'ITS', 'TAR1'})


def is_satellite(feature: str) -> bool:
    """Return True for any v1 or v2 satellite feature name."""
    return feature in SATELLITE_FEATURES


def lookup_satellite_col(row, pfx: str, col_kind: str, name: str, default=0):
    """Read a satellite column from an annotations row, trying v2-then-v1.

    Args:
        row: dict-like row from cluster_analysis / sequence_annotations TSV.
        pfx: featureset prefix, e.g. 'region_subtelomere_flat'.
        col_kind: column kind, e.g. 'dmax', 'readpct', 'bppct', 'max_block_bp'.
        name: canonical (v2) feature name, e.g. 'bSat'.
        default: returned when neither variant is present (default 0).

    Returns:
        The first non-zero value found at `{pfx}_{col_kind}__{name}` (v2) or
        `{pfx}_{col_kind}__{v1_alias}` (v1) — or `default` if both are 0/missing.
    """
    val = row.get(f'{pfx}_{col_kind}__{name}', None)
    if val is None or val == 0:
        v1_alias = SATELLITE_V2_TO_V1.get(name)
        if v1_alias is not None:
            v1_val = row.get(f'{pfx}_{col_kind}__{v1_alias}', None)
            if v1_val is not None and v1_val != 0:
                return v1_val
    if val is None:
        return default
    return val
