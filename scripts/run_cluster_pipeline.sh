#!/usr/bin/env bash
#
# run_cluster_pipeline.sh — Engine B whole-sample structural clustering, end to end.
#
# Reproduces the figures in plots_preview/<SAMPLE>.final.all_clusters.svg starting from
# the per-featureset telogator BEDs (the same inputs we started with). It runs the four
# pipeline stages in order:
#
#   1. bin-annotations     denoise each featureset BED (rolling-window mode filter)
#   2. overlay-annotations resolve the three featuresets into one annotation BED
#   3. cluster             overlap-layout-consensus clustering (Engine B)
#   4. cluster-plot        render every multi-read cluster to one SVG
#
# ---------------------------------------------------------------------------------------
# WHAT YOU NEED ON DISK (per sample)
# ---------------------------------------------------------------------------------------
#   * Three telogator featureset BEDs, named  <PREFIX>.<featureset>.smoothed.features.bed.gz
#     for featureset in { region, subtelomeric, chromosome }. <PREFIX> is whatever comes
#     before ".<featureset>.smoothed.features.bed.gz" — e.g.
#       MySample.telogator.1.KS_human_CHM13_v2
#   * The KaryoScope database directory (the same DB used to make those BEDs), containing
#     hierarchy.tsv and (optionally) colors.tsv.
#   * The genome feature-weights TSV. The CHM13v2 weights are committed at
#     data/chm13v2_feature_weights.tsv — use those if you used the CHM13v2 reference.
#     Only recompute (genome-weights) if you used a different reference.
#
# ---------------------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------------------
#   scripts/run_cluster_pipeline.sh \
#       --sample   MySample \
#       --prefix   /path/to/MySample.telogator.1.KS_human_CHM13_v2 \
#       --db       /path/to/KS_human_CHM13_v2 \
#       [--weights data/chm13v2_feature_weights.tsv] \
#       [--outdir  results/MySample] \
#       [--workers 8]
#
# All flags have defaults (see CONFIG below); --sample and --prefix are the two you must
# almost always set. Run with -h/--help for the option list.
#
set -euo pipefail

# ======================================================================================
# CONFIG — defaults; override any of these with the matching --flag.
# ======================================================================================
SAMPLE=""                                        # sample name, e.g. MySample
PREFIX=""                                        # input path prefix (see above)
DB=""                                            # dir holding hierarchy.tsv (+ colors.tsv)
WEIGHTS="data/chm13v2_feature_weights.tsv"       # genome-freq feature weights TSV
OUTDIR=""                                        # output dir (default: results/<SAMPLE>)
WORKERS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"  # parallel alignment workers

# Clustering knobs (the values used for the published figures — change only if you know why)
WINDOW=1001                                      # bin-annotations rolling window (bp)
MAJORITY_FRACTION=0                              # tau=0 => always descend to a leaf
PRESET="chromosome-telomere-satellite"           # overlay-annotations resolution preset
MIN_OVERLAP_BP=1000                              # edge size + anti-chaining gate (distinctive bp)
BLOCK_MIN_BP=2000                                # blocking index: only align reads sharing >= this
MIN_CLUSTER_SIZE=2                               # plot only clusters with >= this many reads
MIN_SEGMENT_BP=0                                 # drop draw segments shorter than this (0 = keep all)

# ======================================================================================
usage() { sed -n '2,46p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)   SAMPLE="$2"; shift 2 ;;
    --prefix)   PREFIX="$2"; shift 2 ;;
    --db)       DB="$2"; shift 2 ;;
    --weights)  WEIGHTS="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2"; shift 2 ;;
    --workers)  WORKERS="$2"; shift 2 ;;
    -h|--help)  usage 0 ;;
    *) echo "Unknown argument: $1" >&2; usage 1 ;;
  esac
done

# --- validate ------------------------------------------------------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -n "$SAMPLE" ]] || die "--sample is required"
[[ -n "$PREFIX" ]] || die "--prefix is required"
[[ -n "$DB"     ]] || die "--db is required (directory containing hierarchy.tsv)"
[[ -f "$DB/hierarchy.tsv" ]] || die "no hierarchy.tsv in DB dir: $DB"
[[ -f "$WEIGHTS" ]] || die "genome weights file not found: $WEIGHTS"

OUTDIR="${OUTDIR:-results/$SAMPLE}"
mkdir -p "$OUTDIR"

# karyoscope-analysis entry point: prefer the repo venv, else whatever is on PATH.
if [[ -x ".venv/bin/karyoscope-analysis" ]]; then
  KS=".venv/bin/karyoscope-analysis"
elif command -v karyoscope-analysis >/dev/null 2>&1; then
  KS="karyoscope-analysis"
else
  die "karyoscope-analysis not found (activate the venv or 'pip install -e .')"
fi

# colors.tsv is optional — cluster-plot falls back to an auto-palette without it.
COLORS_ARG=()
if [[ -f "$DB/colors.tsv" ]]; then
  COLORS_ARG=(--colors "$DB/colors.tsv")
else
  echo "NOTE: $DB/colors.tsv not found — plot will use an auto-generated palette." >&2
fi

# Confirm the three input featureset BEDs exist before doing any work.
for fs in region subtelomeric chromosome; do
  f="$PREFIX.$fs.smoothed.features.bed.gz"
  [[ -f "$f" ]] || die "missing input featureset BED: $f"
done

echo "=== KaryoScope Engine B pipeline ==============================="
echo "  sample   : $SAMPLE"
echo "  inputs   : $PREFIX.{region,subtelomeric,chromosome}.smoothed.features.bed.gz"
echo "  db       : $DB"
echo "  weights  : $WEIGHTS"
echo "  outdir   : $OUTDIR"
echo "  workers  : $WORKERS"
echo "================================================================"

# ======================================================================================
# 1. bin each featureset (rolling-window mode filter; tau 0 = always descend to a leaf)
# ======================================================================================
for fs in region subtelomeric chromosome; do
  echo ">>> [1/4] bin-annotations: $fs"
  "$KS" bin-annotations \
    --input "$PREFIX.$fs.smoothed.features.bed.gz" \
    --hierarchy "$DB/hierarchy.tsv" \
    --feature-set "$fs" \
    --window "$WINDOW" \
    --majority-fraction "$MAJORITY_FRACTION" \
    -o "$OUTDIR/$SAMPLE.$fs.binned.bed.gz"
done

# ======================================================================================
# 2. overlay the three binned featuresets into one resolved annotation BED
# ======================================================================================
echo ">>> [2/4] overlay-annotations"
"$KS" overlay-annotations \
  --bed chromosome="$OUTDIR/$SAMPLE.chromosome.binned.bed.gz" \
  --bed region="$OUTDIR/$SAMPLE.region.binned.bed.gz" \
  --bed subtelomeric="$OUTDIR/$SAMPLE.subtelomeric.binned.bed.gz" \
  --hierarchy "$DB/hierarchy.tsv" \
  --preset "$PRESET" \
  -o "$OUTDIR/$SAMPLE.overlay.bed"

# ======================================================================================
# 3. cluster ALL reads (genome-freq weights + distinctive gate + blocking + parallel)
#    Writes <out>.clusters.tsv plus sidecars .clusters.consensus.bed and .clusters.layout.tsv
# ======================================================================================
echo ">>> [3/4] cluster"
"$KS" cluster \
  --input "$OUTDIR/$SAMPLE.overlay.bed" \
  --hierarchy "$DB/hierarchy.tsv" \
  --min-length 0 \
  --weight-method genome-freq \
  --genome-weights "$WEIGHTS" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --block-min-bp "$BLOCK_MIN_BP" \
  --workers "$WORKERS" \
  -o "$OUTDIR/$SAMPLE.clusters.tsv"

# ======================================================================================
# 4. plot every multi-read cluster into one SVG
#    --chromosome-track    : chromosome-colored track under each read (translocations show two)
#    --no-consensus-track  : omit the union-consensus top row (clearer read-to-read alignment)
#    These are the flags used for the published plots_preview/<sample>.final.all_clusters.svg.
# ======================================================================================
echo ">>> [4/4] cluster-plot"
"$KS" cluster-plot \
  --layout "$OUTDIR/$SAMPLE.clusters.layout.tsv" \
  --consensus "$OUTDIR/$SAMPLE.clusters.consensus.bed" \
  "${COLORS_ARG[@]}" \
  --min-cluster-size "$MIN_CLUSTER_SIZE" \
  --min-segment-bp "$MIN_SEGMENT_BP" \
  --chromosome-track \
  --no-consensus-track \
  -o "$OUTDIR/$SAMPLE.final.all_clusters.svg"

echo "================================================================"
echo "DONE. Figure: $OUTDIR/$SAMPLE.final.all_clusters.svg"
echo "      Clusters table: $OUTDIR/$SAMPLE.clusters.tsv"
echo "================================================================"
