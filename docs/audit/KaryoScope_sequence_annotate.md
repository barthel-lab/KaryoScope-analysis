# Audit: KaryoScope_sequence_annotate.py (1065 lines)

## 1. Purpose
FACT. Computes per-read (sequence-level) feature annotations from one or more
per-read BED files and writes a WIDE-format TSV with one row per read (module
docstring, lines 2-54). For each `(featureset, feature)` pair it emits coverage
fraction, bp, window-density statistics (max/min/median/first/last/terminal/
terminal_min), and longest contiguous block (lines 8-18). It also emits per-read
"interspersion" transition rates (lines 23-24), optional alignment statistics from
telogator `readnames.txt`/`stats.tsv` (lines 26-29), and an adaptive-threshold
sidecar TSV. Logic is explicitly lifted from `KaryoScope_cluster_annotate.py` /
`KaryoScope_annotate_sequences.py` (section comments, lines 196, 230, 311, 392,
447, 464, 618) but produces PER-READ rows rather than per-cluster aggregates.

ASSESSMENT. This is a feature-matrix builder over already-annotated BEDs — a
downstream consumer of karyoscope's annotate output, not a (re)annotator.

## 2. CLI surface
FACT. Real `main()` with `if __name__ == "__main__"` guard (lines 713, 1064-1065).
`argparse` parser built inside `main` (lines 714-752). Arguments:
- `--bed` (required, `nargs='+'`): BED/BED.gz; featureset auto-detected from
  filename OR given as `label:path` (lines 719-721).
- `--output` / `-o` (required): TSV path, `.gz` → gzip (lines 722-723).
- `--featuresets` (csv, default None): subset of detected featuresets (724-725).
- `--window-size` (int, default 1000): density window bp (726-727).
- `--readnames-dir` (default None): enables alignment-stats join (728-729).
- `--samples` (csv, default None): override auto-extracted sample names (730-732).
- `--reference` (default `CHM13`): used to locate `stats.tsv` (733-734).
- `--colors-dir` + `--database` (both default None, must be paired): pull canonical
  feature ordering from `{db}.{featureset}.colors.txt` (735-741, validated 755-758).
- `--log-file` (default True; truthy-string parser): write `<output>.log` (742-744).

FACT. Argparse defaults are snapshotted into module global `_argparse_defaults`
(lines 73, 746-750) for the parameter-printing table (`_print_params_and_command`,
664-706).

## 3. Inputs & outputs (BED/TSV formats, columns, gzip handling)
FACT. Input BED: 4-column TSV `read,start,end,feature`; a `length = end - start`
column is derived at load (`load_bed_file`, lines 80-102). Lines with `<4` fields
are silently skipped (line 89) — NO malformed-line warning (differs from
merge_beds). gzip detected by `.gz` suffix (lines 82-83). Non-integer start/end
will raise (no try/except, unlike merge_beds line 138).

FACT. Optional inputs (when `--readnames-dir`):
- `{dir}/{sample}/telogator/{sample}.readnames.txt` — 2-col headerless
  `read, sequencing_approach` (`load_readnames`, 467-496).
- `{dir}/{sample}/telogator/aligned/{sample}.{reference}.stats.tsv` — header'd TSV
  with columns `readname/read, is_primary, is_not_supplementary, align_len,
  read_len, mapq, de, align_fraction, is_mapped` (`load_stats`, 499-614).

FACT. Outputs:
- Main wide TSV (`--output`); gzip if `.gz` via `to_csv(compression=...)`
  (lines 1038-1040). One row per read; key column `sequence` (line 844). Column
  naming: `{fs}_frac__{feat}`, `{fs}_bp__{feat}`, `{fs}_dmax__{feat}`, …
  `{fs}_max_block_bp__{feat}`, `{fs}_total_bp`, `interspersion_*`, plus alignment
  columns (lines 919-968, 982-988, 1025-1031).
- Sidecar `{output_prefix}.adaptive_thresholds.tsv` (NOT gzipped, lines 1046-1050).
- `{output_prefix}.log` when `--log-file` (lines 761-764). `output_prefix` strips
  `.tsv`/`.tsv.gz` (`_strip_tsv_extension`, 651-657).
ASSESSMENT. The threshold sidecar is always plain-text even when the main output is
gzipped — minor inconsistency.

## 4. Pipeline / control flow (key functions + line numbers)
FACT. `main()` (713-1061):
1. Build parser, snapshot defaults, parse, validate colors/database pairing
   (714-758).
2. Set up `TeeLogger` (621-644) on stdout if `--log-file` (761-764); print params
   (766-769).
3. Parse `--bed` into `beds_by_featureset` dict via `label:path` split or
   `detect_featureset` (771-796); the disambiguation order is: colon-with-
   nonexistent-path → label; colon-with-existing-path → try auto-detect, else split;
   no colon → auto-detect (lines 775-790).
4. Filter by `--featuresets` (798-808); sort featuresets for reproducibility (811).
5. Load + concat BEDs per featureset (819-829); compute union of read ids (835-839).
6. Per featureset, compute bp/fractions/thresholds, resolve canonical feature order
   (851-910), compute window densities in bulk (912-914), assemble all columns into
   a dict then build a DataFrame (916-970) — explicitly to avoid fragmentation
   (842, 972-974).
7. Interspersion from `telomere_region` else `region` featureset (976-990).
8. Optional alignment-stats join (992-1035).
9. Write main TSV + thresholds sidecar + summary (1037-1061).

Core computation functions:
- `compute_read_feature_fractions` (199-210); `compute_per_read_feature_bp` (451-460).
- `compute_adaptive_thresholds` (213-226).
- `_max_block_length` (236-253) — numpy run-length with gap merging.
- `compute_per_read_window_densities_bulk` (314-388) — numpy event-array coverage,
  sliding-window sums, terminal stats.
- `classify_bed_feature` (270-307) + `compute_per_read_interspersion_bulk`
  (395-444) — typed transition counting per kb.
- `load_canonical_features` (158-192); `detect_featureset` (120-140);
  `extract_sample_name` (143-155).

## 5. Key design decisions
FACT — adaptive thresholds. `compute_adaptive_thresholds` (213-226): per feature,
`threshold = clamp(median_nonzero / 3, 0.001, 0.05)`; features with no nonzero
observations get `min_thresh` (0.001). The `/3` factor and the 0.1%–5% clamp are
hard-coded (lines 213, 225).
ASSESSMENT. No WHY is documented for `/3` or the 0.1%–5% bounds. These thresholds
are emitted to a sidecar but are NOT themselves applied to filter/zero the output
columns in this script (they appear to be advisory for downstream `cluster_*`). The
caller should confirm the threshold is consumed elsewhere.

FACT — block gap tolerance. `BLOCK_GAP_TOL = 100` (line 233): contiguous-block
detection merges runs separated by `<=100 bp` (lines 248, 236-253). Hard-coded.

FACT — window-density branching. Reads with span `< window_size` get a single
whole-read coverage fraction broadcast to all density fields (lines 344-358); reads
`>= window_size` get true sliding-window stats (359-384). Coverage is computed via
`np.add.at` event arrays + `cumsum` (lines 346-371) — vectorized and efficient.

FACT — interspersion classification. `classify_bed_feature` (270-307) handles
single-layer and `a:b` two-layer feature labels with a documented priority
(satellite > layer-2 telomere types > layer-1 telomere types > ct > arm > other,
lines 277-307). Feature membership comes from the SHARED `_feature_vocab` module
(lines 67-70, 262-267), explicitly so both v1 and v2 database vocabularies classify
correctly (comment lines 256-260).
ASSESSMENT. Good — this is the correct pattern that merge_beds.py does NOT follow.

FACT — interspersion source priority: `telomere_region` preferred, else `region`,
else skipped with a warning (lines 976-990).

FACT — alignment-stats semantics (load_stats, 499-614): "total aligned bases" sums
`align_len` only over NON-secondary alignments (primary OR supplementary), with an
inline rationale comment (lines 543-547). Two code paths: full aggregate when
`is_primary/is_not_supplementary/align_len/read_len` all present (532-576), else a
fallback that filters to primary-non-supplementary and back-fills placeholder
columns (577-601). Raises on duplicate reads (605-612).

## 6. Assumptions (checkable statements)
FACT/ASSESSMENT (checkable):
- BED col 1 is a per-read id; all grouping is `groupby('read')`.
- `length = end - start >= 0`; coverage math assumes 0-based half-open and
  `end >= start` (negative lengths would corrupt sums; not validated).
- Window/terminal stats assume the read's own coordinate frame: `read_start =
  min(start)`, `read_end = max(end)` per read (lines 335-337); reads with
  `span <= 0` are skipped (line 339).
- Interspersion assumes intervals are sortable by `start` (sorted at line 411);
  does NOT assume pre-sorted input.
- Filename auto-detection assumes specific patterns:
  `{sample}.(telogator|rDNA_filtered).{n}.{db}.{featureset}.smoothed.*.bed[.gz]`
  or `{sample}.{featureset}.merged.bed` (regexes lines 111-117). Other names MUST
  use `label:path`.
- `extract_sample_name` assumes sample is the substring before `.telogator` /
  `.rDNA_filtered`, else first dot-field (lines 143-155).
- Colors file path assumed `{colors_dir}/{database}.{featureset}.colors.txt`,
  whitespace-delimited, optional `feature` header, `_specific` suffix stripped,
  `_multigroup1` kept (lines 158-192).
- readnames/stats directory layout is the telogator convention (lines 480, 518).
- stats.tsv `is_primary`/`is_not_supplementary` are boolean (compared with `== True`,
  lines 548, 560-562).

## 7. Dependencies
FACT.
- External libs: `pandas` (top-level, line 65); `numpy` lazily imported inside
  `_max_block_length` (238) and `compute_per_read_window_densities_bulk` (326).
- Stdlib: `argparse`, `gzip`, `os`, `re`, `sys`.
- NO pyranges (no interval-join here; all per-read numpy/pandas).
- Inter-script / shared deps: imports `_feature_vocab` (SATELLITE_FEATURES,
  ARM_FEATURES, CT_FEATURES, CANONICAL_TELOMERE, NONCANONICAL_TELOMERE, ITS_TAR1)
  via a `sys.path.insert` on the script dir (lines 67-70). This is the ONLY shared
  dep. It does NOT import colors/karyoplot — `load_canonical_features` reads colors
  files directly (158-192). Duplicates `TeeLogger` (621-644) and `load_bed_file`
  (80-102) from sibling scripts.
- External tools: NONE (no bgzip/tabix/bedtools/subprocess); gzip via stdlib.
- karyoplot relevance: NONE (data processing; no plotting). Note karyoplot DOES
  have `core/colors` and `core/sample_metadata` modules that overlap conceptually
  with the local colors-file parsing — a candidate for future consolidation.

## 8. Proposed home in new layout
ASSESSMENT.
- Subcommand: `karyoscope-analysis annotate-sequences` (or `seq-annotate`). NOTE the
  name collision with karyoscope's `annotate` — keep them clearly distinct.
- Decomposition:
  - `commands/annotate_sequences.py`: thin click command.
  - `core/seq_features.py`: pure metric functions (fractions, bp, thresholds,
    window densities, block length, interspersion, classify).
  - `core/io/bed.py`: shared BED reader (replaces local `load_bed_file`).
  - `core/io/telogator_stats.py`: `load_readnames` / `load_stats`.
  - `core/io/colors.py` or reuse `_feature_vocab`/karyoplot for canonical features.
- DOES IT DUPLICATE karyoscope `core/annotate.py`? NO — IMPORTANT distinction:
  karyoscope `core/annotate.py` annotates a FASTA against a k-mer database
  (`get_featureIDs` C++ helper → feature-id BED → hierarchy-aware smoothing →
  per-featureset `*.smoothed.bed.gz`; see its commands/annotate.py docstring). That
  is the PRODUCER of the BEDs. THIS script CONSUMES those smoothed BEDs and computes
  per-read summary statistics into a wide TSV. The shared word "annotate" is the
  only overlap; there is zero functional duplication. It should NOT be dropped or
  redirected to karyoscope — it stays in the analysis package. Recommend renaming to
  avoid confusion (e.g. `seq-features` / `read-features`) since it computes features
  *matrices*, not annotations.

## 9. Smells / risks / dead code / duplication (line-cited)
FACT/ASSESSMENT.
- Silent skip of malformed BED lines with NO count/warning (line 89) — diverges from
  merge_beds and hides data issues; non-int coords raise uncaught (line 91-93).
- `_argparse_defaults` global + manual default snapshotting (73, 746-750) is fragile
  (click would give this for free).
- Duplicated `TeeLogger` (621-644) and `load_bed_file` (80-102) across sibling
  scripts — consolidate into `core/io`.
- `TeeLogger.write`/`flush` (627-634) do not guard against a closed file (unlike the
  merge_beds copy at line 93); writing after close would raise.
- Adaptive thresholds (`/3`, 0.001–0.05) are computed and written but never applied
  in this script (see §5) — possible dead/advisory output; confirm a consumer.
- Hard-coded magic numbers without rationale: `/3` (225), `0.001`/`0.05` (213),
  `BLOCK_GAP_TOL=100` (233), `window_size` default 1000 (726), `span_kb` /1000 (412).
- `--bed` colon-disambiguation (775-790) is subtle: a real path containing `:` that
  also fails auto-detect falls back to `split(':',1)`, which can mis-split a path.
- `load_stats` boolean comparisons use `== True` (548, 560-562) — brittle if the
  column is read as strings ("True") or 0/1; would silently select nothing.
- `compute_per_read_*_bulk` re-filter `bed_df[bed_df['read'].isin(read_ids)]` each
  call (328, 404) — `read_ids` is always the full union, so the filter is a no-op
  copy (minor inefficiency).
- Per-feature density assembly builds 8 python dicts then `Series.map` per feature
  (939-968) — O(reads*features); fine for typical sizes but the hot loop.

## 10. Testability notes
ASSESSMENT. Highly unit-testable; functions are pure and numeric.
- `_max_block_length(coverage, gap_tol)` (236-253): pure numpy; test gap-merging at
  the 100-bp boundary, empty coverage, single run. HIGH priority.
- `compute_adaptive_thresholds(fractions)` (213-226): pure; test clamp bounds,
  all-zero feature, median/3.
- `classify_bed_feature(feature)` (270-307): pure string→category; enumerate every
  priority branch incl. two-layer `a:b` and v1/v2 names. HIGH priority.
- `compute_per_read_interspersion_bulk` (395-444): deterministic given a small BED
  DataFrame; test transition counting and the typed pairs (can_ncan/tel_sat/arm_tel).
- `compute_read_feature_fractions` / `compute_per_read_feature_bp` (199-210, 451-460):
  pure groupby math; golden small-frame tests.
- `compute_per_read_window_densities_bulk` (314-388): both branches (span<window vs
  >=window); golden vectors.
- `detect_featureset` / `extract_sample_name` / `_strip_tsv_extension`
  (120-155, 651-657): pure string parsing; table-driven incl. ValueError path.
- `load_canonical_features` (158-192): file-based but tiny; tmp-file fixture.
- `load_stats` (499-614): both the aggregate and fallback paths + duplicate-read
  raise; tmp stats.tsv fixtures.
PRE-REQ: functions are already top-level (importable once on path), so testing
mostly needs them extracted from the `_feature_vocab` `sys.path.insert` hack into a
proper package import.

## 11. Open questions for the user
1. The adaptive-threshold formula (`median_nonzero/3`, clamp 0.1%–5%) and
   `BLOCK_GAP_TOL=100` — what is the biological/empirical rationale, and where are
   the thresholds actually consumed (cluster_analysis?)?
2. Should this be renamed to avoid collision with karyoscope `annotate` (e.g.
   `seq-features`)? Confirmed it does NOT duplicate karyoscope core/annotate.
3. Should malformed BED lines be counted/warned (as merge_beds does) instead of
   silently skipped, and should non-integer coords be handled gracefully?
4. Is the telogator `readnames.txt` / `stats.tsv` directory layout (lines 480, 518)
   stable, and should `--reference` default stay `CHM13`?
5. Canonical-feature ordering currently parses colors files directly; should it
   instead use karyoplot `core.colors` / the shared `_feature_vocab` to stay
   consistent with the rest of the ecosystem?
6. Are the two-layer `feat1:feat2` labels (consumed by `classify_bed_feature`)
   exactly the output of merge_beds.py's overlay mode? If so, the two scripts share
   an implicit label contract worth documenting/testing jointly.
