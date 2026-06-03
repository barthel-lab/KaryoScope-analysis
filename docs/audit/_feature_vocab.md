# Audit: _feature_vocab.py (109 lines)

## 1. Purpose
FACT: Module docstring (lines 1-21) states it provides "Shared feature-vocabulary
constants for KaryoScope analysis scripts" to "standardize feature vocabulary
across v1/v2 KaryoScope databases". It reconciles two feature-database
vocabularies:
- `KS_human_CHM13` (v1): `_specific`/`_multigroup1` suffixes, lowercase satellite
  names (`bsat_specific`, `hsat3_specific`, `gsat_specific`) (lines 4-6).
- `KS_human_CHM13_v2` (v2): suffixes dropped, mixed-case bare names (`bSat`,
  `HSat3`, `gSat`, `cenSat`) with finer subdivisions (`active_hor`, `alpha_hor`,
  `mon`, `dhor`) (lines 7-9).

ASSESSMENT: It is a pure constants-plus-helpers module — no I/O, no side effects,
no CLI. The actual scope is narrower than "feature vocabulary" in general: it
covers satellites (with full v1/v2 reconciliation), plus four small hard-coded
sets for arm / centromere-transition / telomere-type features.

## 2. CLI surface
N/A — this is a library module. It has no `argparse`, no `main()`, no
`if __name__ == "__main__"` block.

## 3. Inputs & outputs
FACT: No file or stdin/stdout I/O. The only runtime input is the `row` argument to
`lookup_satellite_col` (line 86), a dict-like mapping (described as a row from
"cluster_analysis / sequence_annotations TSV", lines 91-92) plus string args
`pfx`, `col_kind`, `name`. Output is a looked-up cell value or `default` (line 109).
All other exports are module-level constants evaluated at import time.

## 4. Pipeline / control flow / public API
FACT: Exported symbols (no `__all__` is defined, so "exported" = module-level
public names):

Constants:
- `SATELLITE_V2` (frozenset) — line 24. v2-canonical satellite features.
- `SATELLITE_V1` (frozenset) — line 32. v1 satellite aliases.
- `SATELLITE_FEATURES` (frozenset) — line 40. `SATELLITE_V2 | SATELLITE_V1`.
- `SATELLITE_V1_TO_V2` (dict) — line 47. v1->v2 alias table (13 entries).
- `SATELLITE_V2_TO_V1` (dict) — line 64. Inverse of the above, built by comprehension.
- `ARM_FEATURES` (frozenset) — line 68. `{'p_arm','q_arm','arm','arm_multigroup1'}`.
- `CT_FEATURES` (frozenset) — line 71. `{'ct','ct_specific'}`.
- `CANONICAL_TELOMERE` (frozenset) — line 76. `{'canonical_telomere'}`.
- `NONCANONICAL_TELOMERE` (frozenset) — line 77. `{'noncanonical_telomere'}`.
- `ITS_TAR1` (frozenset) — line 78. `{'ITS','TAR1'}`.

Functions:
- `is_satellite(feature: str) -> bool` — lines 81-83. Returns
  `feature in SATELLITE_FEATURES`. (See section 9: never imported anywhere.)
- `lookup_satellite_col(row, pfx, col_kind, name, default=0)` — lines 86-109.
  Control flow: read `row[f'{pfx}_{col_kind}__{name}']` (v2 key, line 100); if that
  is `None` or `0`, look up the v1 alias via `SATELLITE_V2_TO_V1` (line 102) and try
  `row[f'{pfx}_{col_kind}__{v1_alias}']` (line 104), returning it if non-`None` and
  non-zero (lines 105-106); if the v2 value was `None`, return `default` (lines
  107-108); else return the v2 value (line 109).

## 5. Key design decisions
FACT:
- Permissive union set: `SATELLITE_FEATURES = SATELLITE_V2 | SATELLITE_V1` (line 40)
  is "used wherever the question is 'is this any kind of satellite?'" (line 39).
  Docstring (lines 19-20) states the sets "deliberately" include both v1 and v2
  names "so that classification works regardless of which database produced the BEDs."
- v2-first-then-v1 lookup: `lookup_satellite_col` tries the v2 key first and only
  falls back to the v1 alias (docstring lines 16-17, 87; logic lines 100-106). WHY
  (stated): "back-compat with existing annotation TSVs" (line 17).
- Approximate semantic mappings called out explicitly (lines 42-46): v1 `active`
  ≈ v2 `active_hor`, v1 `monomeric` ≈ v2 `mon`. The comment warns callers needing
  exact substitution to use the table directly (lines 45-46).
- Inverse table built programmatically (line 64) rather than hand-maintained,
  reducing drift between the two directions.
- Telomere-type names "didn't change between v1 and v2 beyond `_specific`
  stripping, which other code already normalizes" (lines 74-75) — WHY only bare
  names are listed in `CANONICAL_TELOMERE`/`NONCANONICAL_TELOMERE`.

ASSESSMENT: The v2-first ordering is sound for a forward-looking codebase. The
"falsy means missing" treatment (`val == 0` triggers v1 fallback, line 101) is a
deliberate but lossy choice — see section 9.

## 6. Assumptions (checkable statements)
1. Annotation TSV columns follow the exact pattern `{pfx}_{col_kind}__{name}`
   (double underscore before the feature name) — lines 100, 104. Checkable against
   real headers produced by `KaryoScope_sequence_annotate.py`.
2. A column value of literal `0` is semantically equivalent to "feature absent",
   so it is safe to fall back to the v1 alias and to never return a genuine `0`
   from the v2 column (lines 101, 105). Checkable: is `0` ever a meaningful result?
3. The v1<->v2 alias table is complete for every satellite a caller will ask about
   by v2 name. Any v2 name absent from `SATELLITE_V2_TO_V1` (e.g. `aSat`, `cenSat`,
   `centromeric`, `HSat`, `dhor`, `mixedAlpha`) has NO v1 fallback (line 102
   returns `None`). Checkable against the database `features.txt` files.
4. The v1/v2 name lists faithfully match the database `features.txt` files
   referenced in comments (lines 23, 30) — these files are the source of truth.
5. `SATELLITE_V1` contains `'censet'` AND `'censat'` (line 34) — `'censet'` is a
   probable typo (see section 9); assumption is that no TSV actually uses `censet`.

## 7. Dependencies
FACT:
- External libs: NONE. No imports at all (pure stdlib data structures).
- karyoplot usage: NONE.
- External tools: NONE.
- Who imports this (within `scripts/`, the only place it is referenced):
  - `KaryoScope_sequence_annotate.py` line 69-70: imports `SATELLITE_FEATURES`,
    `ARM_FEATURES`, `CT_FEATURES`, `CANONICAL_TELOMERE`, `NONCANONICAL_TELOMERE`,
    `ITS_TAR1` (aliased to module-level `_*_LAYER1/_LAYER2_*` names at lines
    262-267, consumed in `classify_bed_feature`, lines ~270-305).
  - `KaryoScope_cluster_annotate.py` line 54: imports `lookup_satellite_col` only
    (used at lines 368, 822, 827).
  - Both importers use `sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))`
    immediately before the import (sequence_annotate line 68; cluster_annotate line
    53) — a path hack required because the module is a sibling script, not a package.
- NOT imported by `KaryoScope_cluster_analysis.py` (the prompt's example is
  inaccurate — grep finds zero references there).
- `KaryoScope_plot_reads.py` defines its OWN unrelated `SATELLITE_FEATURES`
  (line 2350) that does NOT import from this module and uses a different,
  all-v1-style name set (`active`, `inactive`, `divergent`, `monomeric`, `hsat1A`,
  ..., `bsat`, `gsat`, `censat`, `noncentromeric`). This is duplication, not a
  dependency (see section 9).

## 8. Proposed home in new layout
RECOMMENDATION: `karyoscope_analysis/core/feature_vocab.py` (drop the leading
underscore once it is a real package module).

Justification:
- The v1/v2 reconciliation knowledge is specific to KaryoScope feature-database
  vocabularies, which is the analysis domain. It is NOT generic plotting (so it
  does NOT belong in `karyoplot`, which is theme/draw/layout oriented and currently
  has zero dependency on it).
- It is also NOT part of the core `karyoscope` engine: that package OWNS the
  database `features.txt` files, but this module is a *back-compat shim* for
  reading old annotation TSVs in the analysis layer. Pushing it into `karyoscope`
  would couple the engine to legacy analysis-output column conventions.
- Both consumers (`sequence_annotate`, `cluster_annotate`) become subcommands of
  `karyoscope_analysis`, so a shared `core/` module is the natural home and
  eliminates the `sys.path.insert` hack.

CAVEAT: If `karyoscope` ever exposes a canonical v1->v2 feature-name map (it is
the authority on `features.txt`), this module's alias tables should consume that
rather than re-encode it. That is a future de-duplication, not a reason to move
the file now.

## 9. Smells / risks / dead code / duplication
- DEAD CODE: `is_satellite` (lines 81-83) is never imported or called anywhere in
  the repo (grep confirms only the definition). `cluster_annotate.py` line 288
  defines an unrelated *local* boolean `is_satellite` with a different test
  (`any(s in cluster_name for s in ['aSat','bSat','CenSat','HSat','GSat'])`) — note
  even that hard-coded list disagrees with this module's vocabulary (`CenSat`/`GSat`
  capitalization, no `cenSat`/`gSat`).
- DUPLICATION: `KaryoScope_plot_reads.py` line 2350 has a parallel
  `SATELLITE_FEATURES` set that should arguably share this module's vocabulary but
  does not — divergent maintenance risk.
- LIKELY TYPO: `SATELLITE_V1` line 34 lists both `'censet'` and `'censat'`;
  `'censet'` looks like a typo and is unreachable via the alias tables (only
  `'censat'` is in `SATELLITE_V1_TO_V2`).
- SEMANTIC RISK: `lookup_satellite_col` treats a real `0` the same as "absent"
  (line 101 `val == 0`), so a genuine zero in the v2 column silently triggers v1
  fallback, and the function can never *return* `0` from the v2 column (it falls
  through to `default`, which happens to also be `0` by default — masking the
  behavior). If a caller ever passes `default != 0`, a true v2 zero would be
  reported as that non-zero default. Checkable / worth a test.
- PARTIAL ALIAS COVERAGE: many v2 names have no v1 fallback (section 6 #3); callers
  asking for those by v2 name get no back-compat. Not necessarily a bug (those
  features may not exist in v1) but it is undocumented per-name.
- NO `__all__`: every name including the typo'd `'censet'` member and the two dict
  internals is part of the implicit public surface.
- The approximate `active`->`active_hor` / `monomeric`->`mon` mappings (lines
  42-44) are encoded as if exact in the dict; the only guard is a comment.

## 10. Testability notes
- Pure functions with no I/O — trivially unit-testable. `lookup_satellite_col`
  takes a plain dict, so tests need only construct `{f'{pfx}_{col_kind}__{name}': v}`
  rows.
- Recommended unit tests: (a) v2 key present and non-zero -> returned; (b) v2 key
  zero, v1 alias present -> v1 value returned; (c) v2 key zero, v1 alias absent ->
  `default`; (d) v2 key missing entirely -> `default`; (e) a v2 name with no v1
  alias (e.g. `aSat`) -> never crashes, returns v2 value or default; (f) the
  `default != 0` edge case from section 9.
- Constant-level guard tests: assert `SATELLITE_V2_TO_V1 == {v:k for k,v in
  SATELLITE_V1_TO_V2.items()}` (catches future drift); assert every value in
  `SATELLITE_V1_TO_V2` is in `SATELLITE_V2`; consider a golden test that all v1/v2
  names match the database `features.txt` files.
- No mocking, no fixtures, no external tools needed.

## 11. Open questions for the user
1. Is `'censet'` (line 34) a typo for `'censat'`, and is it safe to remove?
2. Should `is_satellite` (lines 81-83) be kept as public API, or deleted as dead
   code? If kept, should `cluster_annotate.py` line 288 use it instead of its
   divergent inline list?
3. Is a column value of `0` ever a meaningful (present-but-zero) result, or is
   "0 == absent" always correct (affects the fallback logic at line 101)?
4. Should `KaryoScope_plot_reads.py`'s separate `SATELLITE_FEATURES` (line 2350) be
   unified with this module during the reorg, or is it intentionally a different
   (orientation-fallback) vocabulary?
5. Does the core `karyoscope` engine already expose a canonical v1->v2 feature
   mapping we should source from instead of re-encoding the alias tables here?
6. Are the `active`->`active_hor` and `monomeric`->`mon` mappings (lines 42-44)
   considered exact enough for production lookups, or should they be flagged/opt-in?
