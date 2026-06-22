# `cluster` — methods audit (rigorous / skeptical)

A step-by-step, deliberately skeptical review of the clustering + enrichment pipeline in
`KaryoScope_cluster_analysis.py` (enrichment mode; structure mode is being dropped). Code refs are line
numbers in the current script. **The goal is to expose every hidden assumption and every statistical
weakness so we can decide what to keep, fix, or rethink.** Severity tags: 🔴 critical, 🟠 high, 🟡 medium.

> Bottom line up front: the most serious issue is **circularity** — the number of clusters `k` is chosen
> partly to *maximize group enrichment/purity*, and then the same clusters are tested for group
> enrichment and reported with p-values. That makes the headline enrichment results **optimistically
> biased by construction**. Details in §4.1.

---

## 1. Proposed decomposition (your "split into separate scripts" request)

The current script is ~3,100 lines of mostly top-level code doing four separable jobs. Proposed split into
**three subcommands + the existing diagnostics consumer**, each independently inspectable and testable:

| New subcommand | Input → Output | Responsibility |
|---|---|---|
| **`build-matrix`** | BEDs → matrix `.npz` (+ feature/edge vocab, seq order) | Adjacency-matrix construction: edges + abundance, transforms, optional SVD. Pure NumPy; no statistics. |
| **`cluster`** | matrix `.npz` → `sequence_assignments.tsv` + linkage `.npz` | Distance + Ward linkage + k-selection + cut. No enrichment. |
| **`test-enrichment`** | assignments + sample-metadata → `cluster_analysis.tsv` | Per-cluster group/sample enrichment (Fisher + FDR). Decoupled from clustering. |
| *(`cluster-diagnostics`, downstream)* | assignments + feature matrix → QC figures | already separate |

**Why this split matters beyond tidiness:** decoupling **`test-enrichment`** from **`cluster`** is the
structural fix for the circularity in §4.1 — k-selection can be made purely unsupervised, and the
enrichment test becomes a clearly-labeled, separate step on clusters that were *not* chosen to maximize
enrichment.

---

## 2. The pipeline, step by step

1. **Load BEDs** → per-`seq_id` ordered feature list (sorted by start). Sample label = first dotted token
   of filename (L501-508).
2. **Filters:** `--sequence-list` whitelist; **span filter `[10000, 50000]`** (L94-96; span computed
   *before* feature exclusion); **`--exclude-features` default `"novel,canonical_telomere*"`** (L100).
3. **Adjacency matrix** (`_build_matrix_from_features` L1042-1184): an **edge block** + an **abundance
   block** (§3).
4. **Optional TruncatedSVD** to `reduce_dims=500` components (L1313-1323, `random_state=42`).
5. **Distance + linkage:** `pdist(euclidean)` + **Ward** (L1401-1402), on the SVD scores if reduction is on.
6. **k-selection** (L1452-1832): per-k loop computing silhouette / cosine-silhouette / Calinski-Harabasz /
   Davies-Bouldin **plus enrichment counts**, a **composite score**, early stopping, then a **knee**;
   `selected_k` per `--k-selection` (default `composite-knee`).
7. **Cut** `fcluster(maxclust=selected_k)` (L1473/1839).
8. **Per-cluster enrichment** (two-group or per-sample Fisher) + centroid read (L1850-1964).
9. **FDR** across clusters; relabel `mixed` if `q ≥ threshold` (L2088-2110).
10. **Outputs:** assignments TSV, cluster_analysis TSV, `feature_matrix.npz`, sample_metadata, UMAP, PDFs.

---

## 3. Adjacency-matrix construction (the feature space being clustered)

**Edge block** (`get_edges` L511-536; columns L1064-1106): one column per *possible* feature pair.
- `symmetric` (default): all `C(n_features, 2)` unordered pairs; pair key alphabetically sorted.
- `directional`: all `n_features·(n_features−1)` ordered pairs.
- Cell = count (or presence, or length-weighted) of that adjacent transition along the sequence.
- 🟠 **Dimensionality explosion:** the edge block is **O(n_features²)** columns, overwhelmingly zero (only
  transitions that actually occur are nonzero). This is why SVD is "recommended" — but it means the raw
  feature space is enormous and sparse.
- 🟡 **Column-key fragility:** symmetric edge columns are named `f"{f1}->{f2}"` in `feature_list` order
  (L1066-1068) while `get_edges` returns alphabetically-sorted pairs. These match **only if `feature_list`
  is alphabetically sorted**. If it isn't, transitions silently fail their column lookup (`if pair_name in
  pair_to_idx`) and are **dropped without warning**. Must assert sorted vocab.
- 🟡 **`length_weighted` is asymmetric/odd:** directional weight = `from_len/read_len` (only the source
  feature's length, L1101); symmetric = `avg(from,to)/read_len` (L1036). Two different definitions.

**Abundance block** (L1112-1140): one column per feature = total bp of that feature on the sequence
(count modes), `/read_len` (length_weighted), or 0/1 (binary).

**Transforms** (matrix-type; default `count_log1p_zscore_blockweight`):
- `log1p` both blocks (L1103-1104, 1142-1143).
- `zscore`: per-column `StandardScaler` on each block (L1147-1151).
- `blockweight`: multiply the **abundance** block by `sqrt(n_edge/n_abund)` (L1160-1161) so the two blocks
  have **equal total variance**. Stated rationale: "equal variance contribution."

🔴/🟠 **Skeptical critique of the default matrix:**
1. 🟠 **Z-scoring sparse count columns is dangerous.** Most edge columns are ~all-zero with a few 1s/2s.
   Z-scoring divides by a tiny std, so a *single* occurrence of a rare transition becomes a huge value and
   ends up **dominating Euclidean distance**. The clustering is then driven by rare, possibly-noise
   transitions rather than by the bulk feature content. This is a well-known pitfall of standardizing
   sparse/binary features.
2. 🟡 **"Equal variance per block" is a heuristic, not a model.** Why should transition structure and
   abundance contribute *equally*? And the `sqrt(n_edge/n_abund)` correction assumes columns are
   independent; edges are highly correlated (co-occurring transitions), so the edge block's *effective*
   dimensionality ≪ `n_edge`, and the equalization over- or under-weights in practice.
3. 🟡 **Ward + Euclidean** assume roughly spherical, equal-variance clusters in the (z-scored, reweighted,
   SVD-reduced) space — there is no reason the data satisfy this; Ward will still return a tree and a cut,
   so "it ran" is not evidence the structure is real.
4. 🟡 **SVD to 500 dims** is arbitrary (no variance-explained criterion); `StandardScaler` densifies the
   matrix first, so SVD runs on a dense `n × O(n_features²)` matrix (memory-heavy).

---

## 4. Statistics — the skeptical core

### 4.1 🔴 Circularity / double-dipping (the headline problem)
The **composite** k-selection score (L1544-1546) is
`0.5·silhouette_norm + 0.1·enriched_ratio + 0.4·perfect_ratio`, where `enriched_ratio`/`perfect_ratio` are
the fraction of clusters that are **group-enriched / group-pure** (`fast_enrichment_check`, odds>1.5).
So **k is chosen partly to maximize group separation**, and then steps 8-9 **test those same clusters for
group enrichment and report p-values/q-values**. The model-selection objective and the hypothesis test are
the same quantity, evaluated on the same data, with no accounting for the selection. Consequences:
- Reported enrichment p/q-values are **optimistically biased** (selection-induced); they do **not** support
  a claim that groups differ.
- `perfect_ratio` uses `perfect_threshold = 1.0` and `min_cluster_size = 3`: a size-3 cluster that is
  coincidentally single-group counts as "perfect." As **k grows, clusters shrink and are more likely to be
  coincidentally pure**, so the composite has a **built-in upward bias in k** and a built-in tendency to
  "find" enrichment.

**Proposed fix:** decouple (the §1 split). Make k-selection **purely unsupervised by default**
(silhouette/CH/DB); keep the enrichment-aware composite only as an explicit, clearly-labeled *exploratory*
option; and run `test-enrichment` as a separate step. If group-difference testing is a goal, it needs a
design that doesn't pick k to maximize the tested effect.

### 4.2 🔴 Per-sample: min-p over samples, then FDR over clusters
`calculate_enrichment_per_sample` (L715-824) runs a one-sided Fisher test for **each sample vs rest**, then
sets the cluster's `p_value` to the **minimum** p across samples (L798-806). FDR (`false_discovery_control`,
L2094) is then applied **across clusters** on these min-p values.
- 🔴 The within-cluster multiple testing across **S samples is never corrected** — taking the min of S
  p-values is the classic maximum-statistic inflation. And the min-of-S statistic is **not a valid p-value**
  (its null is not Uniform), so applying BH/BY to it (L2094) is **statistically invalid**.
- Correct options: correct over the **full sample×cluster grid**, or use a proper per-cluster omnibus test
  (e.g. one Fisher/χ² across all samples), then FDR over clusters.

### 4.3 🟠 Three different definitions of "enriched" in one pipeline
1. **k-selection** (`fast_enrichment_check` L1444-1446): odds-ratio (0.5 pseudocount) **> 1.5**, **no
   p-value at all**.
2. **per-cluster** (`calculate_enrichment_*` L688, 798): raw Fisher **p < 0.05**.
3. **final labels** (L2103): **FDR q < threshold**.
A cluster called "enriched" while choosing k can be "mixed" in the output, and vice versa. The numbers in
the k-selection diagnostic plots and the final table are **not the same quantity**.

### 4.4 🟠 Size-factor "normalization" by rounding
For `telomeric`/`total` normalization (L770-795, and again L1929-1947), counts are divided by a per-sample
size factor and **`round()`ed**, then Fisher is run on the rounded pseudo-counts, with `max(0, …)` clamps.
- 🟠 Rounding **zeros out small samples** and the clamps distort the 2×2 table; Fisher on fabricated
  pseudo-counts is **not a valid test**. A principled compositional approach (proportion test / logistic
  regression with a library-size offset, or CLR) is warranted.
- 🟡 Clarification vs the earlier audit: `telomeric` and `total` are **not** identical — they share the
  size-factor *machinery* but use **different denominators** (`telomeric` = telomeric read counts → equalizes
  telomeric library size; `total` = total genomic reads). They do give different results; but `telomeric`
  normalizing telomeric counts by telomeric counts just forces every sample's scaled total to ≈ the median —
  worth questioning whether that's meaningful.

### 4.5 🟡 One-sided vs two-sided inconsistency
Per-sample Fisher is one-sided `alternative='greater'` (L766, 793) — only tests *enrichment*; two-group
Fisher is **two-sided** (L678). Direction handling differs between modes.

### 4.6 🟡 Bare `except:` hides failures
The per-group enrichment inside per-sample mode catches **all** exceptions and returns `odds, pval = 1.0,
1.0` (L1927, 1946) — any bug silently becomes "not significant." Violates C2 (fail loudly).

### 4.7 🟡 composite-knee instability (the default k-selector)
- Normalization uses the **observed** k-range (which depends on `min_k`, `max_k`, **and where early stopping
  fired**), so `selected_k` is **not robust** to those arbitrary settings. The code itself comments that the
  knee "may vary slightly with max-k" (L1705).
- 🟡 **Diagnostic ≠ decision:** `selected_k` uses the **smoothed** knee (L1816), but the diagnostic plot
  labels the **raw** knee as "Composite-knee k" (L1779) — the plotted k can differ from the chosen k.
- Smoothing window = `max(3, 0.2·k_range)` (L1719) — another arbitrary, range-dependent knob.

---

## 5. Hidden constants to expose as documented parameters (D6)

| Constant | Value | Where | Concern |
|---|---|---|---|
| span filter | `[10000, 50000]` bp | L94-96 | telomere-specific; should be explicit/assay-agnostic |
| `--exclude-features` default | `novel,canonical_telomere*` | L100 | 🟠 **excludes the canonical telomere signal by default** — big, undocumented biological choice |
| `--min-k` | 40 | L79 | 🟠 strong prior that ≥40 clusters exist |
| `--max-k` | 300 | L81 | interacts with knee normalization |
| composite weights | 0.5 / 0.1 / 0.4 | L1544 | arbitrary; and supervised terms cause §4.1 |
| `--perfect-threshold` / `--strong-threshold` | 1.0 / 0.80 | L143-145 | arbitrary purity cutoffs |
| fast-enrichment odds cutoff | 1.5 | L1446 | undocumented |
| `--reduce-dims` | 500 | L160 | no variance criterion |
| `--silhouette-sample-size` | 2000 | L149 | silhouette becomes a noisy subsample estimate |
| `--early-stopping` | 150 | L147 | large; truncates the k-curve, affecting the knee |
| `random_state` | 42 (everywhere) | many | reproducible, but hides run-to-run variance of SVD/UMAP/silhouette |
| `--fdr-threshold` / `--fdr-method` | 0.05 / bh | L173-175 | OK but see §4.2 (input p-values invalid in per-sample) |

---

## 6. Other smells (confirmed)
- ~2,200 lines of **top-level code**, no `main()` — not importable/testable; helpers call `sys.exit()`
  (L417, 1463) which a CLI/test harness can't catch. Must wrap.
- `cluster_to_enrichment` / `enrichment_colors` rebuilt 3-4× across plot blocks (L2279, 2596, 2722, …).
- `groupby().apply(lambda)` (L869, 1012) — slow/deprecated on new pandas.
- NPZ writes a duplicate `read_names` == `seq_names` (L2210-2211) for back-compat — droppable (D1).

---

## 7. Decisions to work through together

1. **Decomposition** — OK with `build-matrix` → `cluster` → `test-enrichment` (+ diagnostics)? (§1)
2. **Circularity (§4.1)** — make k-selection **unsupervised by default**, demote the enrichment-composite to
   an explicit exploratory option, and run enrichment strictly as a separate post-hoc step? This is the big
   methodological call.
3. **Per-sample testing (§4.2)** — replace min-p-over-samples + cluster-FDR with either a full sample×cluster
   FDR grid or a per-cluster omnibus test?
4. **Compositional normalization (§4.4)** — replace the round-and-Fisher hack with a proportion/offset model?
   Keep `telomeric` vs `total` as distinct, or drop `telomeric` (questionable)?
5. **Default k-selector (§4.7)** — switch default off `composite-knee` to silhouette (or Calinski) for
   reproducibility?
6. **Matrix default (§3)** — is z-scoring sparse edge counts acceptable, or should the default be a
   non-standardized / abundance-only / presence-based representation? Do we need the edge block at all, or is
   feature *abundance* (which is interpretable) enough?
7. **`exclude canonical_telomere` default (§5)** — intended? It removes the dominant telomere signal from the
   clustering.
8. **Length window (§5)** — make explicit (you said yes); confirm there's no default that silently filters.

(Each becomes a recorded decision once we work through it.)
