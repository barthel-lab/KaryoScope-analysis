"""Regenerate the committed v2 end-to-end test fixtures from the full raw BEDs.

The end-to-end tests for ``overlay-annotations`` and ``build-feature-matrix`` run
against a *tiny* subset of the real ``KS_human_CHM13_v2`` HeLa annotation BEDs so
they stay fast and run by default. This script carves that subset out of the full
per-sample BEDs in ``data/raw_bed/`` (which are large and not part of the default
test run) and writes the small fixtures next to it.

It is **not** part of the test suite — it is provenance + a regeneration recipe.
The fixtures it writes (``HeLa.v2.<featureset>.bed.gz``) are committed; this script
only needs the full ``data/raw_bed/`` files, which the maintainer has locally.

Selection rule (deterministic, documented in ``README.md``): among reads whose total
interval count across the six featuresets is in ``[MIN_INTERVALS, MAX_INTERVALS]``
(small, but not the degenerate single-feature reads) and that contain at least one
telomere feature in the subtelomeric track (so interspersion is exercised), take the
``N_READS`` smallest by total interval count, ties broken by ``seq_id``. Reads are
kept **whole** (every row, original order) so each fixture still satisfies the C4
partition invariant and all featuresets share the same span per read. Run from
anywhere::

    python tests/data/v2_subset/make_subset.py
"""

from __future__ import annotations

import gzip
from collections import OrderedDict
from pathlib import Path

SAMPLE = "HeLa"
FEATURESETS = ("acrocentric", "chromosome", "gene", "region", "repeat", "subtelomeric")
N_READS = 4
MIN_INTERVALS = 60
MAX_INTERVALS = 220
TELOMERE_FEATURES = frozenset({"canonical_telomere", "noncanonical_telomere", "telomere_like"})

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]  # tests/data/v2_subset -> repo root
RAW = REPO / "data" / "raw_bed"


def _raw_path(featureset: str) -> Path:
    return RAW / f"{SAMPLE}.telogator.1.KS_human_CHM13_v2.{featureset}.smoothed.features.bed.gz"


def _read_groups(path: Path) -> OrderedDict[str, list[tuple[int, int, str]]]:
    groups: OrderedDict[str, list[tuple[int, int, str]]] = OrderedDict()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            seq_id, start, end, feature = line.split("\t")[:4]
            groups.setdefault(seq_id, []).append((int(start), int(end), feature))
    return groups


def main() -> None:
    missing = [str(_raw_path(fs)) for fs in FEATURESETS if not _raw_path(fs).is_file()]
    if missing:
        raise SystemExit(
            "Cannot regenerate fixtures — missing full raw BEDs:\n  " + "\n  ".join(missing)
        )

    all_groups = {fs: _read_groups(_raw_path(fs)) for fs in FEATURESETS}

    # Pick small-but-structured reads: total interval count in band + a telomere
    # present (so the interspersion metric has telomere<->X transitions to count).
    read_ids = set(all_groups["region"])
    sizes = {rid: sum(len(all_groups[fs][rid]) for fs in FEATURESETS) for rid in read_ids}

    def has_telomere(rid: str) -> bool:
        return any(feat in TELOMERE_FEATURES for _, _, feat in all_groups["subtelomeric"][rid])

    eligible = [
        rid
        for rid in read_ids
        if MIN_INTERVALS <= sizes[rid] <= MAX_INTERVALS and has_telomere(rid)
    ]
    chosen = sorted(sorted(eligible), key=lambda r: (sizes[r], r))[:N_READS]
    if len(chosen) < N_READS:
        raise SystemExit(
            f"Only {len(chosen)} eligible reads found (wanted {N_READS}); widen the band."
        )

    print(f"Selected {len(chosen)} reads (smallest by total interval count):")
    for rid in chosen:
        span = all_groups["region"][rid][-1][1]
        print(f"  {rid}  span={span}bp  total_intervals={sizes[rid]}")

    for fs in FEATURESETS:
        out = HERE / f"{SAMPLE}.v2.{fs}.bed.gz"
        groups = all_groups[fs]
        n_rows = 0
        with gzip.open(out, "wt", newline="") as fh:
            for rid in chosen:  # whole reads, in the chosen (sorted) order
                for start, end, feature in groups[rid]:
                    fh.write(f"{rid}\t{start}\t{end}\t{feature}\n")
                    n_rows += 1
        print(f"  wrote {out.name}: {n_rows} rows")


if __name__ == "__main__":
    main()
