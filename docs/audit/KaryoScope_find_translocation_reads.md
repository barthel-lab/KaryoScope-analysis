# Audit: KaryoScope_find_translocation_reads.py (327 lines)

## 1. Purpose

FACT: Auto-discovers per-sample **chromosome** translocation BED files produced
by the core `karyoscope` engine, parses each gzipped BED to compute, per read,
the total read length and per-target-chromosome base-pair coverage, and writes a
single combined TSV plus a console summary (module docstring, lines 2-14;
`main`, lines 152-323).

ASSESSMENT: This is the **find** (stage 1) step of the
find → cluster → visualize trio. Its TSV output is the catalog of candidate
translocation reads consumed by stage 3 (`visualize`). It is a pure
discover/parse/aggregate utility with no plotting.

## 2. CLI surface

FACT (`argparse`, lines 153-177):

- `--results-dir` (required) — base results dir with per-sample subdirs.
- `--output-dir` (required) — output dir for the TSV.
- `--output-prefix` (default `translocation_reads`).
- `--database` (default `KS_human_CHM13`).
- `--translocation-types` (nargs+, default `chr1_chr21 chr2_chr13`).
- `--target-chromosomes` (nargs+, default `chr1 chr2 chr13 chr21`).
- `--sample-prefix` (default None → all sample dirs).
- `--log-file` / `--no-log-file` (`BooleanOptionalAction`, default True;
  writes `{output_dir}/{output_prefix}.log`).

FACT: Uses `argparse`, not click — diverges from the gold-standard click CLI in
`KaryoScope/src/karyoscope/cli.py`.

## 3. Inputs & outputs (formats passed between the three tools)

FACT — Inputs:
- Glob `*.*.*.{database}.*.chromosome.presmoothed.translocations.bed.gz`
  (line 62), matched against regex
  `^(.+?)\.(.+?)\.(\d+)\.{database}\.(chr\d+_chr\d+)\.chromosome\.presmoothed\.translocations\.bed\.gz$`
  (lines 58-60). Groups: `sample`, `data_type`, `replicate`, `trans_type`.
- Each BED is gzipped, tab-separated, ≥4 columns; columns consumed
  (`parse_chromosome_bed`, lines 104-119): col0=`read_id`, col1=`start`,
  col2=`end`, col3=`feature`. NOTE: col0 is the **read id**, not a chromosome —
  this is a per-read feature BED keyed by read, not a genome-coordinate BED.

FACT — Output (lines 255-277): one TSV `{output_dir}/{output_prefix}.tsv` with
dynamic header `read_id, sample, data_type, replicate, read_length`, then per
target chromosome `{chrom}_bp` and `{chrom}_pct` (2 decimals), then
`translocation_type` (lines 258-261). Reads with `read_length == 0` are dropped
(`build_results`, lines 129-130).

INTER-TOOL FACT: This TSV is the exact input contract for `visualize`'s
`--input-tsv` mode (`load_translocation_reads`,
`KaryoScope_visualize_translocation_reads.py` lines 190-204), which reads
columns `read_id, sample, data_type, replicate, read_length,
translocation_type`. The `_bp`/`_pct` columns are NOT consumed downstream.

## 4. Pipeline / control flow (key functions + line numbers)

1. `discover_combos` (46-88) — iterate sorted `results_dir` subdirs, optional
   `sample_prefix` filter, `rglob` the BED glob, regex-match, filter to
   requested `trans_types`, collect combo dicts, sort.
2. `parse_chromosome_bed` (91-121) — per read accumulate `length = max(end)` and
   `chr_coverage[feature] += end - start` for features in `target_set`.
3. `build_results` (124-149) — emit one row per non-zero-length read with
   `_bp`/`_pct` per target chromosome.
4. `main` (152-323) — validate, banner, discover, per-combo parse+build+extend,
   write TSV, print summary table and grand totals.

## 5. Key design decisions (cite lines)

FACT: Read length = max `end` seen across all rows for that read
(lines 115-116), not a stored read length. ASSESSMENT: assumes BED features tile
read coordinates from 0 and that the rightmost `end` equals the true read length.

FACT: `_pct` = `100 * bp / read_length` (line 144). FACT: coverage is summed
without merging overlaps (line 119) → overlapping features can inflate
`_bp`/`_pct` beyond 100%.

FACT: "Translocation read" detection is delegated upstream — this script does
NOT itself detect translocations; it only catalogs reads already present in the
engine-emitted `*.translocations.bed.gz` files (no breakpoint logic anywhere).
WHY: not stated in code.

FACT: `translocation_types` is used both as a discovery filter and as the
summary table's column set (lines 76, 281-313).

## 6. Assumptions (checkable)

- Filename regex (58-60) is rigid: `replicate` must be all digits (`\d+`),
  `trans_type` must match `chr\d+_chr\d+` (so `chrX_chr21` would NOT match).
- BED is gzipped, tab-delimited, ≥4 cols, col0=read_id (lines 104-110).
- `target_chromosomes` values match the literal strings in BED col3 (`feature`)
  (line 118) — chrom naming must agree exactly (e.g. `chr1`).
- Only the `chromosome` featureset is discovered (`presmoothed` + `.chromosome.`
  in glob, line 62); other featuresets are ignored here.
- Every read of interest appears with at least one row having `end > 0`.

## 7. Dependencies

FACT: Stdlib only — `argparse, gzip, re, sys, collections.defaultdict,
pathlib.Path` (lines 16-21). No pandas, no karyoplot, no subprocess, no
samtools/bedtools/rsvg. ASSESSMENT: the lightest of the trio.

## 8. Proposed home in new layout

ASSESSMENT:
- Subcommand: `karyoscope-analysis find-translocation-reads`.
- Thin command wrapper in `commands/find_translocation_reads.py`; discovery +
  parsing + row-building logic in `core/translocation.py` (or
  `core/find_translocation.py`).
- `discover_combos`'s filename-pattern parsing should live in a shared
  `core/io/result_layout.py` (the same `sample.data_type.replicate.database...`
  pattern is re-implemented in `cluster` and `visualize`).
- karyoplot push-down: BED reading should route through
  `karyoplot.core.io.load_bed` / `iter_bed_records` instead of the local
  `gzip.open` loop (lines 104-119).
- `TeeLogger` (lines 27-43) → shared `core/logging.py` (duplicated in 9 scripts).

## 9. Smells / risks / dead code / duplication

- DUP: `TeeLogger` (27-43) is byte-identical across all three scripts (and 6
  others).
- DUP: filename regex / discovery loop (46-88) overlaps heavily with
  `cluster.discover_translocation_beds` (different featureset, same shape).
- RISK: overlapping features double-count coverage; `_pct` can exceed 100
  (line 119, 144).
- RISK: read length inferred from max `end` (115-116) is wrong if features don't
  start at 0 or don't reach the read end.
- SMELL: redundant existence re-check `if not bed_path.exists()` (line 239) —
  paths just came from `rglob`.
- DEAD-ish: `_bp`/`_pct` columns are written but unused by downstream tools
  (only `read_length` and key columns are consumed by `visualize`).
- SMELL: `TeeLogger` opens a file but is never closed; `sys.stdout` is
  monkey-patched globally (lines 31, 192).
- SMELL: early `return` on no combos (line 225) skips the "Done" banner — minor
  inconsistency.

## 10. Testability notes

ASSESSMENT: Highly testable — pure functions. `parse_chromosome_bed` and
`build_results` take in-memory/file inputs and return plain dicts/lists; easy
golden tests on a tiny gzipped BED. `discover_combos` testable with a tmp dir
tree. Blockers: `main` mixes IO, formatting, and stdout side effects; `TeeLogger`
global stdout patch complicates capture; no return value from `main` to assert.

## 11. Open questions for the user

1. Should coverage merge overlapping intervals before summing `_bp` (to cap
   `_pct` ≤ 100)? Is current double-counting intentional?
2. Are `_bp`/`_pct` columns consumed by any tool outside this trio (they appear
   unused by `visualize`)? Can they be dropped or are they a downstream contract?
3. Is the `chr\d+_chr\d+` restriction acceptable, or must sex chromosomes
   (`chrX_chr21`) be supported in `trans_type`?
4. Should read length come from an authoritative source (e.g. BAM/engine
   metadata) rather than max-`end` inference?
