#!/usr/bin/env bash
#
# compute_genome_weights.sh — recompute the genome-frequency feature weights for Engine B.
#
# Tallies how much of the annotated CHM13 reference each feature covers (one BED per featureset)
# and writes information-content weights in (0, 1] — ubiquitous features (chromosome arm) -> ~0,
# rare/distinctive ones (a satellite, telomere) -> 1. See docs/audit/rearrangement_detection.md §12
# and src/karyoscope_analysis/core/genome_weights.py.
#
# The committed data/chm13v2_feature_weights.tsv is the derived artifact; clustering uses it
# directly and never needs to recompute. Re-run this ONLY when the reference changes or to add a
# featureset (just extend --featuresets / add its reference BED below).
#
# ---------------------------------------------------------------------------------------
# WHAT YOU NEED
# ---------------------------------------------------------------------------------------
#   * The annotated CHM13 reference BEDs, one per featureset, named
#       <REF_PREFIX>.<featureset>.smoothed.bed.gz
#     e.g. data/raw_bed/chm13v2.0.KS_human_CHM13_v2.region.smoothed.bed.gz . These are gitignored
#     (~475 MB) and needed ONLY here — copy them to the machine to recompute.
#   * The KaryoScope database directory with hierarchy.tsv (validates feature names; C2).
#
# Each reference BED must be C4 (chrom start end feature, seq_id = chromosome) and cover the WHOLE
# genome for that featureset — the weight uses each feature's fraction of the featureset total, so a
# partial reference inflates every fraction.
#
# ---------------------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------------------
#   scripts/compute_genome_weights.sh \
#       --db /path/to/KS_human_CHM13_v2 \
#       [--ref-prefix data/raw_bed/chm13v2.0.KS_human_CHM13_v2] \
#       [--featuresets "region subtelomeric chromosome repeat gene acrocentric"] \
#       [-o data/chm13v2_feature_weights.tsv]
#
# --db is the only flag you usually must set. Run with -h for the option list.
#
set -euo pipefail

# ======================================================================================
# CONFIG — defaults; override with the matching --flag.
# ======================================================================================
DB=""                                                       # dir holding hierarchy.tsv
REF_PREFIX="data/raw_bed/chm13v2.0.KS_human_CHM13_v2"       # reference BED path prefix
FEATURESETS="region subtelomeric chromosome repeat gene acrocentric"  # one --bed per featureset
OUTPUT="data/chm13v2_feature_weights.tsv"                   # output weights TSV

# ======================================================================================
usage() { sed -n '2,38p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db)          DB="$2"; shift 2 ;;
    --ref-prefix)  REF_PREFIX="$2"; shift 2 ;;
    --featuresets) FEATURESETS="$2"; shift 2 ;;
    -o|--output)   OUTPUT="$2"; shift 2 ;;
    -h|--help)     usage 0 ;;
    *) echo "Unknown argument: $1" >&2; usage 1 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -n "$DB" ]] || die "--db is required (directory containing hierarchy.tsv)"
[[ -f "$DB/hierarchy.tsv" ]] || die "no hierarchy.tsv in DB dir: $DB"

# karyoscope-analysis entry point: prefer the repo venv, else whatever is on PATH.
if [[ -x ".venv/bin/karyoscope-analysis" ]]; then
  KS=".venv/bin/karyoscope-analysis"
elif command -v karyoscope-analysis >/dev/null 2>&1; then
  KS="karyoscope-analysis"
else
  die "karyoscope-analysis not found (activate the venv or 'pip install -e .')"
fi

# Build a --bed FEATURESET=PATH argument per featureset, checking each reference BED exists.
bed_args=()
for fs in $FEATURESETS; do
  bed="$REF_PREFIX.$fs.smoothed.bed.gz"
  [[ -f "$bed" ]] || die "missing reference BED for featureset '$fs': $bed"
  bed_args+=(--bed "$fs=$bed")
done

echo "=== genome-weights ============================================="
echo "  featuresets : $FEATURESETS"
echo "  ref prefix  : $REF_PREFIX"
echo "  db          : $DB"
echo "  output      : $OUTPUT"
echo "================================================================"

"$KS" genome-weights "${bed_args[@]}" --hierarchy "$DB/hierarchy.tsv" -o "$OUTPUT"

echo "DONE. Weights: $OUTPUT"
