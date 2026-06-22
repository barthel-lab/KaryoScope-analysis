# KaryoScope-analysis — Known Issues & v1→v2 Vocabulary Errors

**Status:** findings from a code audit (June 2026), prepared for lab review.
**Confidence note:** items below were found by *reading* the code against the definitive v2
database files, not by running every pipeline. Items marked **(verify)** should be confirmed with a
quick test run before treating as settled. Everything else is a direct reading of the source.

**Definitive references used:**
- v2 taxonomy / palette: `KaryoScope-databases/KS_human_CHM13_v2/{hierarchy,features,colors}.tsv`
- v1→v2 rename map: `karyoscope_feature_name_mapping.tsv`

---

## 0. Background: the two feature vocabularies

KaryoScope's feature names changed between database versions:

- **v1** (`KS_human_CHM13`) used suffixes: `_multigroup1` (internal taxonomy nodes) and `_specific`
  (leaf features), with lowercase satellite names (`bsat`, `hsat3`, `gsat`, …).
- **v2** (`KS_human_CHM13_v2`) dropped both suffixes and uses mixed-case bare names (`bSat`, `HSat3`,
  `gSat`, `cenSat`, …) with a finer taxonomy (`alpha_hor`, `active_hor`, `dhor`, `hor`, `mon`, …).

**Suffix nuance (matters for any v1 support / migration):** in v1, the **presmoothed** BED files carried
*both* `_specific` and `_multigroup1` suffixes, while the **smoothed** BED files had the `_specific`
suffix stripped. So a correct v1→v2 mapping must match each v1 name **both with and without** its
`_specific` suffix.

**Decision taken:** standardize the entire analysis codebase on **v2 only**; drop v1 from all analysis
logic. (A separate one-shot v1→v2 file migrator may be built later if old result files must be reread.)

---

## 1. Vocabulary errors

| ID | Severity | Location | Issue |
|----|----------|----------|-------|
| V1 | **High** | `_feature_vocab.py:54` | Wrong mapping: `hor_multigroup1 → hor`. Correct: `hor_multigroup1 → alpha_hor`. (v2 `hor` actually comes from v1 `inactive_specific`.) |
| V2 | **High** | `_feature_vocab.py:47-61` | Missing mappings: `inactive_specific → hor` and `divergent_specific → dhor` are absent (both v1 names are listed in `SATELLITE_V1` but never mapped). |
| V3 | Medium | `_feature_vocab.py:47-61` | Inconsistent v1 key forms: some keys drop `_specific` (`bsat`, `gsat`, `censat`, `hsat1A/B/2/3`, `active`, `monomeric`); others keep `_multigroup1`. Real v1 columns carry their suffix, so the bare-name keys may never match (esp. presmoothed files — see suffix nuance). |
| V4 | Low | `_feature_vocab.py:33` | Typo `'censet'` in `SATELLITE_V1` (should be `censat`). Harmless (never matches) but indicates the set was hand-edited. |
| V5 | Medium | `_feature_vocab.py` (whole module) | Incomplete scope: only **satellite/arm** names are mapped. The definitive map covers all six feature sets (~87 renames), e.g. `noncentromeric_specific → rDNA`, `repeat_multigroup1 → repeat`, `autosome_multigroup1 → autosome`, `chr13_specific → chr13`. None of these are handled. |
| V6 | **High** | `KaryoScope_plot_reads.py:2350-2355` | A *second, private* `SATELLITE_FEATURES` set that is **100% v1 names** (`bsat`, `gsat`, `censat`, `hsat1A/B/2/3`, `active`, `inactive`, `divergent`, `monomeric`, `noncentromeric`). On v2 input it matches **nothing**. Also `noncentromeric` is semantically rDNA in v2, not a satellite. |
| V7 | **High** | `KaryoScope_merge_beds.py:184,180,186-189` | Priority/background/acrocentric feature sets use **v1-only** literals (`arm_multigroup1`, `telomere_like_multigroup1`, `array_multigroup1`, `acrocentric_multigroup1`). On v2 input the corresponding priority-merge rules **silently never fire**. |
| V8 | Medium | `_feature_vocab.py:101` | `lookup_satellite_col` treats a real value of **`0` as "absent"** and falls through to the v1 alias; it can never return a legitimate v2 `0`. |

### Downstream effects (which scripts are affected)

A count of hard-coded v1-suffix literals per script (`grep -c "_specific\|_multigroup1"`):

| Script | v1-suffix refs | Effect of the v1 vocabulary errors |
|---|---:|---|
| `_feature_vocab.py` | 14 | Source of V1–V5, V8. |
| `KaryoScope_plot_reads.py` | 10 | **V6**: satellite-dense read-orientation fallback never triggers on v2 data → reads may be drawn in the wrong orientation. (verify) |
| `KaryoScope_merge_beds.py` | 8 | **V7**: priority/background merges don't apply to v2 features → **wrong merged feature labels**, which then propagate to clustering, annotation, and every plot downstream. (verify) |
| `KaryoScope_sequence_annotate.py` | 6 | Consumes `_feature_vocab` constant sets for classification; carries V3/V4 baggage. On v2-only data the v2 names work, but the module mixes v1+v2. |
| `KaryoScope_cluster_annotate.py` | (via import) | Consumes `lookup_satellite_col` → inherits **V1** (reads `hor` instead of `alpha_hor` on v1 fallback) and **V8** (0-as-absent) → wrong satellite columns → **wrong cluster auto-labels, representative selection, and feature-importance**. (verify) |
| `KaryoScope_cluster_plot.py` | 5 | Private `p_arm_specific/q_arm_specific → arm` label-collapsing logic; v1-suffix handling becomes inert on v2 (low impact) but is duplicated, drift-prone logic. |
| `KaryoScope_draw_legend.py` | 1 | Minor `_specific` suffix strip. |

**Net:** the two highest-impact items are **V7 (merge_beds)** — because wrong merges corrupt everything
downstream — and **V1/V8 (cluster_annotate via `_feature_vocab`)** — because they silently mislabel
satellite content used for cluster naming. **V6 (plot_reads)** affects figure orientation, not the
underlying data.

---

## 2. Database-file inconsistency

Cross-checking the mapping file against the three v2 db files (script in the appendix), **every v2
feature in `hierarchy.tsv`, `colors.tsv`, and `features.tsv` is present in the mapping** — with one
exception:

- The feature **`repeat`** appears in `colors.tsv` (`repeat	repeat	#B0C4DE`) and in the mapping
  (`repeat_multigroup1 → repeat`), but it is **not a node in `hierarchy.tsv`** and is **not used as a
  leaf in `features.tsv`**. In the hierarchy, the `repeat` feature set is rooted at
  `Interspersed_Repeat` / `Satellite` / `Noninterspersed` / `nonrepeat` (all children of `categorized`);
  there is no `repeat` node.

**Question for the db owner:** is `repeat` supposed to be a hierarchy node (e.g. the parent of the
repeat subtree), or is the `repeat` row in `colors.tsv`/mapping stale and should be removed? This is the
one place the files disagree — consistent with a post-mapping edit to the db files.

---

## 3. Code correctness bugs (independent of vocabulary)

| ID | Severity | Location | Issue | Effect |
|----|----------|----------|-------|--------|
| B1 | **High** | `KaryoScope_compare_clusterings.py:84` | Imports `adjusted_rand_index` (correct name: `adjusted_rand_score`); the `ImportError` is silently caught | **ARI has likely never been computed** in any comparison report (NMI is unaffected). (verify) |
| B2 | Medium | `KaryoScope_cluster_analysis.py:~355` | structure mode does `groupby(['read', …])` but the column is `sequence` | structure mode likely broken for multi-read chromosomes. (verify) |
| B3 | — | (same as **V7** above) | merge priority sets are v1-only | see §1 |
| B4 | **High** | `KaryoScope_visualize_translocation_reads.py` | If `rsvg-convert` is missing, the SVG is deleted in a `finally` block | the run produces **no output at all**, silently |
| B5 | Medium | `KaryoScope_cluster_analysis.py` | Three inconsistent significance rules across the pipeline (k-loop odds>1.5; raw p<0.05; FDR q<thresh); per-sample FDR applied only over per-cluster min-p | inconsistent / under-corrected enrichment calls |
| B6 | Medium | `KaryoScope_select_representatives.py` | `centroid_distance` is computed but never drives selection; the required `--cluster-analysis` arg is unused | the option/name misleads; it does nothing |

---

## 4. Systemic issues (summary)

- **Hard-coded single-experiment biology**: e.g. `cluster_diagnostics` hard-codes `'E6E7'`/`'primary'`
  enrichment groups; `cluster_annotate` has a sample-name-specific "Type I ALT" relabel. These won't
  generalize to other experiments. (Decision: parameterize/remove.)
- **Silent failures**: broad `except Exception: pass`; malformed BED lines skipped without warning;
  unseeded jitter (non-reproducible figures).
- **Dead code / no-op flags**: several documented CLI flags do nothing (e.g. `cluster_plot
  --enrichment-normalization total`, `create_panning_animation --uniform-zoom`). (Decision: remove.)
- **Three separate feature vocabularies** (`_feature_vocab`, `plot_reads`, `merge_beds`) instead of one.

---

## 5. Recommended remediation

1. **Single source of truth for vocabulary**, derived from the definitive `hierarchy.tsv` (+ palette
   from `colors.tsv`) rather than hand-copied Python sets — eliminates V1–V8 by construction and removes
   hard-coded biology.
2. **v2-only**: delete v1 names, the alias tables, and the v1-fallback branch. Migrate `plot_reads` (V6)
   and `merge_beds` (V7) off their private literals.
3. **Fix B1–B6** as each script is migrated, locking the corrected behavior with a test.
4. **Resolve the `repeat` db inconsistency** with the database owner.

---

## Appendix: mapping cross-check (reproducible)

```python
# Compares karyoscope_feature_name_mapping.tsv against the v2 db files.
import os
DB="…/KaryoScope-databases/KS_human_CHM13_v2"; MAP="…/karyoscope_feature_name_mapping.tsv"
def rows(p):
    return [l.rstrip("\n").split("\t") for l in open(p)
            if l.strip() and not l.lstrip().startswith("#")]
hier=set(); [hier.update((r[1],r[2])) for r in rows(f"{DB}/hierarchy.tsv")[1:]]
colors={r[1] for r in rows(f"{DB}/colors.tsv")[1:]}
mv2={r[1] for r in rows(MAP)[1:] if r[1]!="N/A"}
print("missing from mapping:", sorted((hier|colors)-mv2))   # -> []  (complete)
print("in mapping not in hierarchy:", sorted(mv2-hier))      # -> ['repeat']
print("in colors not in hierarchy:", sorted(colors-hier))    # -> ['repeat']
```

Result: the mapping is **complete** on the v2 side; the only discrepancy is `repeat` (present in
`colors.tsv` + mapping, absent from `hierarchy.tsv`).
