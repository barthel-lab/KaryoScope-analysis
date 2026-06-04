# `v2_subset` — end-to-end test fixtures

A tiny, committed subset of the real **`KS_human_CHM13_v2`** HeLa annotation BEDs,
used by the end-to-end tests for `overlay-annotations` and `build-feature-matrix`
(`tests/test_overlay_annotations_e2e.py`, `tests/test_build_feature_matrix_e2e.py`).

These fixtures let the e2e tests run **fast and by default**. The full per-sample
BEDs in `data/raw_bed/` are large and exercised only by the `@pytest.mark.integration`
tests (deselected by default; run with `pytest -m integration`).

## Files

Six annotation BEDs (one per featureset), gzip-compressed like the real data:

```
HeLa.v2.acrocentric.bed.gz
HeLa.v2.chromosome.bed.gz
HeLa.v2.gene.bed.gz
HeLa.v2.region.bed.gz
HeLa.v2.repeat.bed.gz
HeLa.v2.subtelomeric.bed.gz
```

Each is a 4-column annotation BED (`seq_id  start  end  feature`) satisfying the
C4 invariant (rows grouped by `seq_id`; each sequence a gapless, non-overlapping
partition). All six share the same four `seq_id`s and the same per-read span.

## Provenance

Carved from `data/raw_bed/HeLa.telogator.1.KS_human_CHM13_v2.<featureset>.smoothed.features.bed.gz`
by `make_subset.py`. Reads are kept **whole** (every row, original order).

**Selected reads** (the four smallest, by total interval count across the six
featuresets, among reads with 60–220 total intervals that carry a telomere feature
— small but structurally rich, so overlay resolution and interspersion are
exercised):

| seq_id | span (bp) | total intervals |
|---|---|---|
| `ee18d619-9d89-4858-b420-b6fd840107af` | 3391 | 63 |
| `e48de5f7-0a89-4242-aa39-f233414ec711` | 3136 | 65 |
| `m84132_240112_233844_s3/234947230/ccs` | 6966 | 72 |
| `m84132_240112_233844_s3/133367142/ccs` | 8597 | 73 |

(The two read-id styles — UUID and PacBio `…/ccs` with slashes — are deliberate:
they confirm `seq_id`s containing `/` round-trip through the tools.)

## Regenerating

Requires the full `data/raw_bed/` v2 BEDs locally (maintainer-only):

```
python tests/data/v2_subset/make_subset.py
```

The selection is deterministic, so regenerating reproduces these exact fixtures.
If `data/raw_bed/` changes, regenerate and review the diff before committing.
