# Audit: KaryoScope_cluster_translocation_reads.py (519 lines)

## 1. Purpose

FACT: An **orchestrator** that auto-discovers per-sample translocation BED files
across multiple featuresets (`region`, `subtelomeric`, `chromosome`), groups
samples by data-type group, merges `region`+`subtelomeric` featuresets per
sample, concatenates replicates, builds per-group sample-metadata TSVs, then
**shells out** to `KaryoScope_cluster_analysis.py` and
`KaryoScope_cluster_plot.py` once per `(group, trans_type, featureset-set)`
(module docstring lines 2-15; `main` lines 212-515).

ASSESSMENT: This is the **cluster** (stage 2) step. Unlike stages 1 and 3 it does
no parsing/plotting itself — it is glue around three sibling scripts
(`KaryoScope_merge_beds.py`, `KaryoScope_cluster_analysis.py`,
`KaryoScope_cluster_plot.py`).

## 2. CLI surface

FACT (`argparse`, lines 213-249):

- `--results-dir` (required, `Path`).
- `--output-dir` (required, `Path`).
- `--scripts-dir` (required) — dir containing the 3 sibling scripts it calls.
- `--colors-dir` (required) — passed through to cluster_plot.
- `--database` (default `KS_human_CHM13`).
- `--translocation-types` (default `chr1_chr21 chr2_chr13`).
- `--featuresets` (default `region subtelomeric chromosome`).
- `--sample-metadata` (required) — TSV with columns `sample, group, color`.
- `--group-config` (optional JSON: `data_type_to_group`, `group_params`).
- `--sample-prefix` (default None).
- `--output-prefix` (default `translocation`).
- `--log-file` / `--no-log-file` (`BooleanOptionalAction`, default True).

FACT: `argparse`, not click — diverges from gold standard.

## 3. Inputs & outputs (formats passed between the three tools)

FACT — Inputs discovered (`discover_translocation_beds`, 69-103):
glob `*.{database}.*.*.presmoothed.translocations.bed.gz`, regex
`^(.+?)\.(.+?)\.(\d+)\.{database}\.(chr\d+_chr\d+)\.(\w+)\.presmoothed\.translocations\.bed\.gz$`
(76-78). Note this captures a **5th group** `featureset` (`\w+`) that
`find`'s regex hard-codes to `chromosome` — so the two discovery regexes differ
only in the featureset slot.

FACT — Other inputs: `--sample-metadata` TSV (`sample, group, color`; validated
294-299); optional `--group-config` JSON (280-291).

FACT — Outputs (under `output_dir/{group}/`):
- `merged_beds/{trans}/{sample}.{trans}.region_subtelomeric.merged.bed.gz`
  (414-416).
- `plot_beds/{trans}/{sample}.{trans}.{region|subtelomeric}.bed.gz` (419-425).
- `chromosome_beds/{trans}/{sample}.{trans}.chromosome.bed.gz` (428-430).
- `samples.tsv` per group (443-444).
- Whatever `cluster_analysis` writes under `out_prefix`
  (`{output_prefix}_{trans}` / `..._chromosome`) and the cluster_plot SVG
  `{out_prefix}.cluster_plot.svg` (465-503).

INTER-TOOL FACT: This stage does NOT read the stage-1 TSV. It independently
re-discovers BED files from `results_dir`. So `find` and `cluster` are siblings
fed by the same engine outputs, not a strict linear pipe. The merged/concatenated
BEDs are consumed only by the cluster_analysis/cluster_plot subprocesses.

## 4. Pipeline / control flow (key functions + line numbers)

1. `discover_translocation_beds` (69-103) → `{(sample,dt,rep,trans): {featureset: Path}}`.
2. `main` Step 1 (328-358): discover, then re-group into
   `groups[group][trans][sample] = [{data_type, replicate, beds}]` (341-349).
3. Step 2 (360-430): in a `TemporaryDirectory`, per sample: `merge_region_subtelomeric`
   (106-114) calls `KaryoScope_merge_beds.py`; `concatenate_beds` (117-124)
   concatenates replicate parts into final region/subtel/merged/chromosome BEDs.
4. Step 3 (432-445): `create_samples_tsv` (127-151) filters metadata to found
   samples, fills missing as `group='unknown', color='#999999'`.
5. Step 4 (447-507): per `(group, trans)`, `run_cluster_analysis` (154-186) on
   merged region+subtel (needs ≥2 samples), then `run_cluster_plot` (189-209) if
   rc==0; same again for chromosome-only.

## 5. Key design decisions (cite lines)

FACT: Fixed `data_type → group` map (`DEFAULT_DATA_TYPE_TO_GROUP`, 53-59) and
per-group clustering params (`DEFAULT_GROUP_PARAMS`, 62-66) baked in; overridable
via `--group-config` JSON (280-291). WHY: long vs fragmented vs HiFi reads need
different length windows / `min_k`.

FACT: cluster_analysis is invoked with a hard-coded analysis recipe
(`run_cluster_analysis`, 157-175): `--comparison-mode two-group`,
`--control-group primary`, `--exclude-features canonical_telomere*,novel,unknown`,
`--k-selection composite-knee`, `--min-cluster-size 3`, `--reduce-dims 500`,
`--umap`, `--circular-dendrogram`, `--background both`. ASSESSMENT: none of these
are exposed as CLI flags here — the recipe is frozen in code.

FACT: cluster_plot invoked with `--n-per-cluster 5 --background both`
(`run_cluster_plot`, 192-203).

FACT: Analyses are skipped unless ≥2 samples have the relevant BED set
(`len(beds) >= 2`, lines 464, 488). WHY (implied): clustering/two-group
comparison needs ≥2 samples.

FACT: Data-type filtering happens at discovery — combos whose `data_type` is not
a key in `data_type_to_group` are dropped (lines 97-98).

## 6. Assumptions (checkable)

- The 3 sibling scripts exist under `--scripts-dir` with those exact names and
  CLI contracts (106-114, 154-186, 189-209). The `--control-group primary`
  assumes a group literally named `primary` exists in the metadata.
- BED filename regex (76-78): `replicate` all digits, `trans_type` =
  `chr\d+_chr\d+`, `featureset` = `\w+`.
- `--sample-metadata` TSV has `sample, group, color` (294-299).
- `region` and `subtelomeric` featuresets both present per entry to produce a
  merged BED (line 397); otherwise no merged output for that entry.
- `python` interpreter = `sys.executable` can run the sibling scripts (110, 159,
  193).
- `KaryoScope_merge_beds.py` accepts `--bed a b --output out` (108-113).

## 7. Dependencies

FACT: stdlib `argparse, gzip, json, re, subprocess, sys, tempfile,
collections.defaultdict, pathlib` (17-25) + **pandas** (line 30, used for
metadata read/filter/concat 135-151, 294). FACT: **subprocess** shells out to 3
sibling Python scripts (`merge_beds`, `cluster_analysis`, `cluster_plot`). FACT:
No karyoplot import. No samtools/bedtools/rsvg directly (those, if any, live in
the called scripts). ASSESSMENT: heaviest inter-script coupling of the trio.

## 8. Proposed home in new layout

ASSESSMENT:
- Subcommand: `karyoscope-analysis cluster-translocation-reads`.
- Thin wrapper `commands/cluster_translocation_reads.py`; orchestration logic in
  `core/translocation_cluster.py`.
- CRITICAL: replace the three `subprocess.run(sys.executable, sibling.py, ...)`
  calls (114, 183, 206) with **in-process function calls** once `cluster_analysis`,
  `cluster_plot`, and `merge_beds` are migrated into the same package
  (`core/`). This removes `--scripts-dir` entirely.
- Shared `core/io/result_layout.py` for the discovery regex (shared with `find`,
  `visualize`).
- `concatenate_beds` (117-124) and `create_samples_tsv` (127-151) → shared
  `core/io/` (metadata logic overlaps `karyoplot.core.sample_metadata`, which
  cluster_plot/cluster_analysis already delegate to per the Phase-13 notes).
- karyoplot push-down: sample-metadata load/filter/merge → consolidate on
  `karyoplot.core.sample_metadata`; BED concatenation could use
  `karyoplot.core.io.smart_open`.
- `TeeLogger` (33-49) → shared `core/logging.py`.

## 9. Smells / risks / dead code / duplication

- DUP: `TeeLogger` (33-49) identical across the trio.
- DUP: discovery regex/loop (69-103) near-identical to `find.discover_combos`.
- RISK: frozen cluster_analysis recipe (157-175) — no override path for any of
  ~12 hard-coded flags; tuning requires editing source.
- RISK: `subprocess.run(..., capture_output=True)` in `merge_region_subtelomeric`
  (114) with `check=True` swallows stdout/stderr — on failure the traceback hides
  the merge error message.
- INCONSISTENCY: `run_cluster_analysis` (183) does NOT capture output (streams to
  console) but `merge` (114) does — uneven error visibility.
- SMELL: `--control-group primary` is hard-coded (162) but `create_samples_tsv`
  invents `group='unknown'` for missing samples (145) — a metadata gap silently
  yields a group that two-group mode won't treat as control.
- SMELL: deeply nested `defaultdict(lambda: defaultdict(lambda: defaultdict(...)))`
  (341, 363-365) is hard to follow.
- SMELL: `import pandas` placed after module code at line 30 (post-docstring,
  post-`_original_command`) — unconventional ordering.

## 10. Testability notes

ASSESSMENT: Moderately testable. `discover_translocation_beds`, `concatenate_beds`,
`create_samples_tsv` are pure-ish (tmp dir / DataFrame in, files out) and unit-
testable. The orchestration in `main` is hard to test because it spawns
subprocesses; tests would need to mock `subprocess.run` or stub the sibling
scripts. Migrating siblings in-process (see §8) would make end-to-end testing
with tiny fixtures feasible. No return value from `main` to assert; relies on
console summary.

## 11. Open questions for the user

1. Should the frozen cluster_analysis recipe (lines 157-175) become configurable
   flags / part of `--group-config`, or stay fixed for reproducibility?
2. Is `--control-group primary` a hard requirement (must a group literally named
   `primary` always exist)? How should runs without it behave?
3. After migration, can the three sibling scripts be called in-process (dropping
   `--scripts-dir`), or must subprocess isolation be preserved (e.g. matplotlib
   backend, memory)?
4. Should this stage consume stage-1's TSV (to restrict to confirmed
   translocation reads) instead of independently re-globbing the engine outputs?
5. Is the `≥2 samples` gate (464, 488) the right threshold, or should single-
   sample clustering be allowed in some modes?
