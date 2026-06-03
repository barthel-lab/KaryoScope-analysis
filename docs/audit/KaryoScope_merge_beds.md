# Audit: KaryoScope_merge_beds.py (1281 lines)

## 1. Purpose
FACT. Merges multiple per-read BED featuresets "by position overlay", combining the
feature labels of overlapping intervals into a single output BED (header comment
lines 1-16). The basic mode produces composite labels by string-joining feature
names with a separator (e.g. `region:chromosome`). Several specialized
*priority* merge modes instead pick a single winning feature per position using
hard-coded biology-specific priority rules (telomere/satellite/acrocentric).
The output is intended as input to `KaryoScope_cluster_analysis.py` (line 14-16).

ASSESSMENT. This is fundamentally an interval-algebra tool over per-read
coordinate systems (the BED "chrom" column is a *read name*, not a genome
chromosome). The "merge" word is overloaded: it is an overlay/intersection join,
not a bedtools-style interval union.

## 2. CLI surface
FACT. Top-level `argparse` parser built at module import time (lines 27-81); there
is no `main()` and no `if __name__ == "__main__"` guard — all logic runs at import.
Arguments:
- `--bed` (required, `nargs='+'`): 2+ BED files, merged in order (lines 31-32).
- `--output` / `-o` (required): output BED, `.gz` triggers gzip (lines 33-34).
- `--separator` / `-s` (default `:`): label join separator (lines 35-36).
- `--reduce-features N` (int, default None): keep top-N features, collapse rest to
  `other` (lines 37-39).
- `--feature-filter` (`nargs='*'`, `INDEX:KEEP:COLLAPSE`): per-BED prefix filter
  before merge (lines 40-47).
- `--telomere-satellite-merge` (flag): 2-BED priority mode (lines 48-51).
- `--priority-merge` (flag): 3-BED subtel>region>repeat mode (lines 52-60).
- `--chromosome-acrocentric-merge` (flag): 2-BED mode (lines 61-65).
- `--telomere-acrocentric-merge` (flag): 2-BED composite-label mode (lines 66-70).
- `--collapse-non-acrocentric` (flag): remap chrom labels before merge (lines 71-77).
- `--log PATH`: tee stdout to a log file (lines 78-79).

FACT. The five mode flags are mutually exclusive only by convention; they are
checked in sequence (lines 1042, 1075, 1124, 1169) and each `sys.exit(0)`s on
completion, so the *first* matching flag wins and later ones are silently ignored
if multiple are passed.

## 3. Inputs & outputs (BED/TSV formats, columns, gzip handling)
FACT. Input BED is 4-column, tab-separated: `read`, `start` (int), `end` (int),
`feature` (`load_bed_file`, lines 118-148). Lines with `<4` tab fields or
non-integer start/end are counted as malformed and skipped (lines 128-140), with a
warning. gzip input is detected purely by `.gz` suffix (lines 120-121, mode `rt`).
Columns 5+ (strand, score) are ignored.

FACT. Output is the same 4-column headerless BED. Two different write paths exist:
- Priority modes use `pandas.to_csv(..., header=False, compression='gzip' if .gz)`
  (e.g. lines 1066-1069, 1115-1118, 1160-1163, 1202-1205).
- The default basic-merge path writes row-by-row with an explicit
  `gzip.open`/`open` loop (lines 1273-1278).
ASSESSMENT. Two divergent serialization paths for the same logical output — a
duplication/consistency smell; the row-by-row loop (lines 1276-1278) is also slow
for large outputs vs vectorized `to_csv`.

FACT. gzip handling is suffix-based only; no bgzip/tabix indexing is produced or
required.

## 4. Pipeline / control flow (key functions + line numbers)
FACT. Module-level procedural flow (no `main()`):
1. Parse args + optional `TeeLogger` setup (lines 81-111); require >=2 BEDs (113-115).
2. Define loaders/helpers/merge functions (lines 118-1014).
3. Parse `--feature-filter` into `feature_filter_map` (lines 1020-1035).
4. Branch on mode flags, each terminating with `sys.exit(0)`:
   - telomere-satellite (lines 1042-1072) → `telomere_satellite_merge` (488-608).
   - 3-way priority (lines 1075-1121) → `priority_merge_three_way` (270-444).
   - chromosome-acrocentric (lines 1124-1166) → `chromosome_acrocentric_merge` (611-731).
   - telomere-acrocentric (lines 1169-1208) → `telomere_acrocentric_merge` (770-936).
5. Default iterative pairwise overlay (lines 1210-1281): load BED1, optional remap
   (1212-1213) and filter (1217-1219); fold remaining BEDs left-to-right via
   `merge_two_beds` (lines 1226-1244); optional `--reduce-features` (1253-1266);
   sort + write (1268-1278).

Core interval functions:
- `merge_two_beds` (939-945) → `_merge_pyranges` (948-982) / `_merge_pandas`
  (985-1014): per-read pairwise overlap, label = `f1{sep}f2`.
- `subtract_intervals` (447-485): single-interval minus sorted blockers — the pure
  workhorse for all pandas-fallback priority modes.
- `apply_conditional_region_repeat_rules` (229-267): region+repeat label resolution.
- `_resolve_telomere_acro_feature` (734-767): subtel+acro label resolution.
- `_merge_adjacent_intervals` (790-807): coalesce same-read same-feature touching
  intervals (used only by telomere-acrocentric mode).

## 5. Key design decisions
FACT — pyranges with pandas fallback. Every merge mode tries `import pyranges`
inside a `try/except ImportError` and falls back to a pure-pandas implementation
(lines 282-286, 502-506, 624-628, 783-787, 941-945). The pandas fallbacks print
"Install pyranges for faster processing" (e.g. line 395).
ASSESSMENT. This is duplicated five times — each mode maintains two parallel
implementations (~10 functions) that must stay behaviorally identical. High
maintenance burden and a correctness-divergence risk (see §9).

FACT — overlay semantics. Basic merge is an inner join on overlap: only positions
where BOTH BEDs have an interval survive (`_merge_pyranges` uses `pr1.join(pr2)`,
lines 966-982; pandas does an O(n*m) nested loop, lines 998-1012). Overlap coords
are `max(starts)`/`min(ends)` (lines 972-974, 1004-1005). Non-overlapping regions
are dropped.

FACT — priority strategy. Priority modes keep designated high-priority intervals
verbatim and SUBTRACT them from the lower-priority layer (pyranges `.subtract()`:
lines 367, 543, 669; pandas via `subtract_intervals`). 3-way mode additionally
JOINs region+repeat first, applies `apply_conditional_region_repeat_rules`, then
subtracts subtel priority (lines 320-380).

FACT — hard-coded biological vocabularies (lines 177-206): `TELOMERE_PRIORITY_FEATURES`
(178), `SUBTEL_PRIORITY_FEATURES` (181), `REGION_BACKGROUND_FEATURES` (184),
`ACROCENTRIC_PRIORITY_FEATURES` (187-190), `TELOMERE_ACRO_TOP_PRIORITY` (193),
`TELOMERE_ACRO_COMPOSITE_SUBTEL` (196), `TELOMERE_ACRO_COMPOSITE_ACRO` (199),
`STRICT_NON_ACROCENTRIC_CHROMOSOMES` (202-206). Conditional rules (229-267): e.g.
`ct + nonrepeat -> ct` else repeat; `noncentromeric + rRNA -> rRNA` else `rDNA`;
arm features are background → use repeat. Composite labels like `DJ_TAR1` are
formed when TAR1/ITS overlaps DJ/PHR/rDNA (lines 754-756).
ASSESSMENT. No rationale ("WHY") is documented for the specific thresholds/feature
membership; these encode lab-specific CHM13 biology. They are NOT sourced from the
shared `_feature_vocab` module (unlike sequence_annotate.py) — they are private,
v1-flavored literals (`arm_multigroup1`, `telomere_like_multigroup1`), so they may
silently miss v2 names (`arm`). See §6/§11.

## 6. Assumptions (checkable statements)
FACT/ASSESSMENT (each checkable):
- BED column 1 is a per-read identifier; intervals are compared only within the
  same `read` value (all merges intersect read-sets first, e.g. line 952).
- BED is 0-based half-open `[start, end)`; `end > start` is required for an overlap
  to count (lines 974, 1006). Zero-length intervals are dropped.
- Input need NOT be pre-sorted: pandas fallbacks sort priority coords locally
  (line 423) and outputs are sorted before writing; pyranges sorts internally.
  EXCEPTION: `_merge_adjacent_intervals` sorts internally (line 795), so OK.
- Within a single BED, intervals of the *same* read may overlap; the basic overlay
  does NOT pre-merge them, so duplicated/overlapping inputs inflate output rows.
- Feature priority sets assume v1-style CHM13 naming (`*_multigroup1`); v2 bare
  names (`arm`, satellites) are NOT in these literals → may be misclassified.
- `--feature-filter` matching is prefix-based via `str.startswith` (line 165).
- chrom naming: `STRICT_NON_ACROCENTRIC_CHROMOSOMES` assumes `chrN`/`chrX`/`chrY`
  literal labels appear as the *feature* value in a chromosome BED (lines 202-217).

## 7. Dependencies
FACT.
- External libs: `pandas` (top-level import, line 22); `pyranges` OPTIONAL
  (lazy `import pyranges` in each merge fn; fallback to pandas on ImportError).
  Listed in `requirements.txt` as `pyranges>=0.0.129  # optional`.
- Stdlib: `argparse`, `atexit`, `datetime`, `gzip`, `sys`.
- NO numpy (unlike sequence_annotate).
- Inter-script / shared deps: NONE. Does not import `_feature_vocab`, karyoscope,
  or karyoplot. Defines its own feature-vocabulary literals (lines 177-206) and its
  own `TeeLogger` (85-103) and `load_bed_file` (118-148) — both duplicated from
  sequence_annotate.py / cluster scripts.
- External tools: NONE (no bgzip/tabix/bedtools/subprocess). gzip via stdlib only.
- karyoplot relevance: NONE — pure data processing, no plotting.

## 8. Proposed home in new layout
ASSESSMENT.
- Subcommand: `karyoscope-analysis merge-beds`.
- Decomposition:
  - `commands/merge_beds.py`: thin click command; mode flags become a
    `--mode {overlay,telomere-satellite,priority,chromosome-acrocentric,
    telomere-acrocentric}` choice (mutually-exclusive enforcement) instead of five
    independent boolean flags.
  - `core/io/bed.py`: shared `read_bed`/`write_bed` (replaces both local
    `load_bed_file` copies and the two divergent write paths).
  - `core/intervals.py`: pure `subtract_intervals`, `merge_adjacent_intervals`,
    `overlay_two`, `priority_subtract` — pyranges-backed with a pure-pandas/python
    fallback chosen ONCE behind a single dispatcher, not duplicated per mode.
  - `core/merge_rules.py`: the biological priority sets + `apply_conditional_*` /
    `_resolve_telomere_acro_feature` resolvers (ideally sourced from
    `_feature_vocab` so v1/v2 names both work).
- The feature-vocabulary literals (lines 177-206) should be reconciled with /
  moved into the shared `_feature_vocab` module to fix the v1-only naming gap.
- This script does NOT overlap with karyoscope `core/annotate.py` (that is
  FASTA→BED k-mer annotation; this is BED→BED overlay). Keep in analysis package.

## 9. Smells / risks / dead code / duplication (line-cited)
FACT/ASSESSMENT.
- No `main()` / no `__main__` guard; all work at import time (lines 81+). Blocks
  unit testing of the orchestration and makes the module unimportable for reuse
  without side effects. (Contrast sequence_annotate.py, which has `main()`.)
- Five duplicated pyranges/pandas implementation pairs (§5) — large surface for
  drift. The pandas fallbacks are O(n*m) per read (nested loops, e.g. lines
  1002-1012, 430-442) vs pyranges interval trees.
- DEAD import: `import pyranges as pr` at line 283 inside `priority_merge_three_way`
  is immediately followed by `return _priority_merge_pyranges(...)`; the `pr` bound
  there is unused (the real import is re-done at line 291). Harmless but dead.
- Divergent output serialization: row-by-row loop (1273-1278) vs `to_csv`
  (priority modes) — inconsistent and the loop is slow.
- `_merge_adjacent_intervals` is applied ONLY in telomere-acrocentric mode (lines
  866, 922); other modes can emit adjacent same-feature fragments (e.g. after
  subtract), producing more output rows than necessary — inconsistent behavior.
- `apply_conditional_region_repeat_rules` (229-267) silently returns the region
  feature for any unrecognized region label (fall-through, line 267); a typo'd
  label passes through unflagged.
- Priority-set literals are v1-only (§5/§6) — a v2 database could silently route
  features through the wrong branch with no warning.
- `_merge_adjacent_intervals` uses `df.iloc[1:].iterrows()` (line 799) — slow on
  large frames.
- `TeeLogger.close()` (lines 102-103) closes the log but `write` guards against a
  closed file (line 93) while `flush`/`__exit__` are absent here (this copy has no
  context-manager, unlike the sequence_annotate copy) — minor inconsistency.

## 10. Testability notes
ASSESSMENT. These interval functions are the MOST unit-testable code in the repo —
pure, deterministic, small inputs. Specific targets:
- `subtract_intervals(interval, blockers)` (447-485): pure; table-driven tests
  (blocker before/after/inside/spanning/multiple). HIGH priority.
- `apply_conditional_region_repeat_rules(region, repeat)` (229-267): pure mapping;
  enumerate all documented rules + fall-through.
- `_resolve_telomere_acro_feature(subtel, acro)` (734-767): pure; test all 5
  priority branches incl. composite `DJ_TAR1`.
- `_merge_adjacent_intervals(df)` (790-807): pure DataFrame→DataFrame; test
  touching vs gapped vs different-feature.
- `apply_feature_filter`, `remap_non_acrocentric_chromosomes` (151-174, 209-226):
  pure DataFrame transforms.
- `_merge_pandas` / each `_*_pandas` (985-1014, 383-444, 559-608, 685-731,
  883-936): deterministic given small DataFrames; ideal golden-output tests, and
  can be cross-checked against their `_*_pyranges` twins (parity tests) to lock the
  two implementations together.
PRE-REQ: extract these out of the import-time script into an importable module
(see §8) before they can be imported by a test.

## 11. Open questions for the user
1. Are the priority-feature sets (lines 177-206) authoritative/frozen, and should
   they cover v2 names (`arm`, `bSat`, etc.)? They are currently v1-only and NOT
   from `_feature_vocab` — was that intentional?
2. Should the five mode flags become one mutually-exclusive `--mode` option, and is
   it OK that today multiple flags silently let the first win?
3. Is the inner-join overlay (dropping positions present in only one BED) the
   intended semantic, or should non-overlapping regions be retained?
4. Should adjacent same-feature intervals be coalesced in ALL modes (currently only
   telomere-acrocentric)? Does downstream `cluster_analysis` care about fragmentation?
5. Is the pure-pandas fallback still needed, or can we make pyranges a hard
   dependency and delete ~5 duplicated functions?
6. Should output be bgzip+tabix-indexed for downstream tooling, or is plain gzip
   sufficient (current behavior)?
