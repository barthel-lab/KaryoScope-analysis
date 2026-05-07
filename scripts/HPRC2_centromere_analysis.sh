#!/usr/bin/env bash
# =============================================================================
# HPRC2_centromere_analysis.sh
#
# End-to-end driver for the HPRC Release 2 pangenome centromere clustering
# analysis (Figure 6 of the KaryoScope paper). Runs every step from
# KaryoScope's per-chromosome BED outputs through to the three publication SVGs.
#
# Inputs (treat KaryoScope as upstream data source):
#   - per-chromosome pre-smoothed region BEDs:
#       <PER_CHR_BED_DIR>/pangenome.chr{1..22,X,Y}.centromere.KS_human_CHM13.presmoothed.region.bed
#   - NucFlag QC summary TSV with per-(chrom, contig) flag counts:
#       <QC_TSV>
#   - KaryoScope satellite-class color palettes:
#       <COLORS_DIR>  (e.g. KaryoScope/resources/databases/KS_human_CHM13)
#
# Outputs (under <OUT_DIR>, default = $ANALYSIS_DIR/agent_results):
#   - pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed       (concat)
#   - pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed  (QC pass)
#   - allchr_structure.feature_matrix.npz, .sequence_assignments.tsv,
#       .k_selection.pdf
#   - allchr_dendrogram.svg                  (clean, intermediate)
#   - allchr_outlier_annotations.tsv         (intermediate)
#   - allchr_dendrogram_annotated.svg        FIGURE 6A (final)
#   - allchr_barplot.svg                     FIGURE 6C (final)
#   - allchr_allele_heatmap.svg, allchr_allele_heatmap_2x2.svg   FIGURE 6D
#
# Usage (local):
#   ./scripts/HPRC2_centromere_analysis.sh
#
# Usage (SLURM):
#   sbatch scripts/HPRC2_centromere_analysis.sh
#
# Environment overrides (any of):
#   KARYOSCOPE_DIR=/path/to/KaryoScope
#   ANALYSIS_DIR=/path/to/KaryoScope-analysis
#   PER_CHR_BED_DIR=...
#   QC_TSV=...
#   COLORS_DIR=...
#   OUT_DIR=...
# =============================================================================

# ---- SLURM directives (ignored on local execution) --------------------------
#SBATCH --job-name=HPRC2_centromere
#SBATCH --output=logs/HPRC2_centromere_%j.out
#SBATCH --error=logs/HPRC2_centromere_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
# Edit for your cluster:
##SBATCH --partition=defq
##SBATCH --account=<your_account>

set -euo pipefail

# ---- Path resolution --------------------------------------------------------
# Repo root: prefer SLURM_SUBMIT_DIR if running under sbatch from the repo,
# otherwise resolve from the location of this script.
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="${ANALYSIS_DIR:-$(cd "${SCRIPT_PATH}/.." && pwd)}"
SCRIPTS_DIR="${ANALYSIS_DIR}/scripts"

# Data inputs (adjust for your cluster filesystem)
KARYOSCOPE_DIR="${KARYOSCOPE_DIR:-/Users/ychen/Documents/GitHub/KaryoScope}"
PER_CHR_BED_DIR="${PER_CHR_BED_DIR:-${KARYOSCOPE_DIR}/local_data/centromere_region_beds}"
QC_TSV="${QC_TSV:-${KARYOSCOPE_DIR}/local_data/aggregate_qc_v4.manual_curation.tsv}"
COLORS_DIR="${COLORS_DIR:-${KARYOSCOPE_DIR}/resources/databases/KS_human_CHM13}"

# Outputs
OUT_DIR="${OUT_DIR:-${ANALYSIS_DIR}/agent_results}"
mkdir -p "${OUT_DIR}" "${ANALYSIS_DIR}/logs"

# ---- Environment ------------------------------------------------------------
# Adapt to your cluster: uncomment the right lines.
# module load R/4.3 python/3.10
# source ~/miniconda3/bin/activate karyoscope

# ---- Pre-flight checks ------------------------------------------------------
required_chroms=( $(seq 1 22) X Y )
for c in "${required_chroms[@]}"; do
    f="${PER_CHR_BED_DIR}/pangenome.chr${c}.centromere.KS_human_CHM13.presmoothed.region.bed"
    [[ -f "$f" ]] || { echo "ERROR: missing input bed: $f" >&2; exit 2; }
done
[[ -f "$QC_TSV"           ]] || { echo "ERROR: missing QC TSV: $QC_TSV"   >&2; exit 2; }
[[ -d "$COLORS_DIR"       ]] || { echo "ERROR: missing colors dir: $COLORS_DIR" >&2; exit 2; }
command -v Rscript  >/dev/null || { echo "ERROR: Rscript not on PATH"     >&2; exit 2; }
command -v python3  >/dev/null || { echo "ERROR: python3 not on PATH"     >&2; exit 2; }

# ---- Banner -----------------------------------------------------------------
echo "=========================================================================="
echo " KaryoScope HPRC centromere analysis"
echo " Started:           $(date '+%F %T')"
echo " Host:              $(hostname)"
echo " Job ID (SLURM):    ${SLURM_JOB_ID:-n/a}"
echo " ANALYSIS_DIR:      $ANALYSIS_DIR"
echo " PER_CHR_BED_DIR:   $PER_CHR_BED_DIR"
echo " QC_TSV:            $QC_TSV"
echo " COLORS_DIR:        $COLORS_DIR"
echo " OUT_DIR:           $OUT_DIR"
echo "=========================================================================="

# ---- Step 1: Concatenate per-chr BEDs -> ALLchr.bed -------------------------
echo
echo "[1/7] Concatenating per-chr BEDs into ALLchr.bed..."
# concatenate_centromere_regions.R reads from the hard-coded relative path
# "../centromere_region_beds/" and writes its output to cwd. To satisfy that,
# we cd into a directory that has a *sibling* called centromere_region_beds
# pointing at PER_CHR_BED_DIR. The simplest setup is a scratch dir with two
# subdirs: "data" (a symlink to PER_CHR_BED_DIR named centromere_region_beds)
# and "run" (where we cd and Rscript writes its output).
ALLCHR_BED_NAME="pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed"
ALLCHR_BED="${OUT_DIR}/${ALLCHR_BED_NAME}"
CONCAT_TMP="$(mktemp -d -t HPRC2_concat.XXXXXX)"
trap 'rm -rf "${CONCAT_TMP}"' EXIT
ln -sfn "$PER_CHR_BED_DIR" "${CONCAT_TMP}/centromere_region_beds"
mkdir -p "${CONCAT_TMP}/run"
(
    cd "${CONCAT_TMP}/run"
    Rscript "${SCRIPTS_DIR}/concatenate_centromere_regions.R"
    mv "${ALLCHR_BED_NAME}" "${ALLCHR_BED}"
)
trap - EXIT
rm -rf "${CONCAT_TMP}"
echo "    -> ${ALLCHR_BED}  ($(wc -l < "$ALLCHR_BED") rows)"

# ---- Step 2: NucFlag QC filter -> pass.bed ----------------------------------
echo
echo "[2/7] Applying NucFlag QC filter (Err=0, COLLAPSE=0, COLLAPSE_VAR=0)..."
PASS_BED="${OUT_DIR}/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed"
(
    # QCfilter_apply.R writes <basename input>.pass.bed and .pass.report.txt
    # to the cwd, so we run it inside OUT_DIR.
    cd "$OUT_DIR"
    Rscript "${SCRIPTS_DIR}/QCfilter_apply.R" \
        "${ALLCHR_BED}" \
        "${QC_TSV}" \
        centromere \
        Err=0 COLLAPSE=0 COLLAPSE_VAR=0
)
echo "    -> ${PASS_BED}  ($(wc -l < "$PASS_BED") rows)"

# ---- Step 3: Per-chromosome clustering --------------------------------------
echo
echo "[3/7] Per-chromosome clustering (Ward linkage, structure mode, k in [2,10])..."
python3 "${SCRIPTS_DIR}/KaryoScope_cluster_analysis.py" \
    --bed "${PASS_BED}" \
    --output-prefix "${OUT_DIR}/allchr_structure" \
    --analysis-mode structure \
    --edges directional \
    --no-abundance \
    --max-sequence-length 50000000 \
    --exclude-features "novel" \
    --matrix-type count_log1p_zscore_blockweight \
    --k-selection silhouette \
    --background white
echo "    -> ${OUT_DIR}/allchr_structure.{feature_matrix.npz, sequence_assignments.tsv, k_selection.pdf}"

ASSIGNMENTS="${OUT_DIR}/allchr_structure.sequence_assignments.tsv"

# ---- Step 4: Clean dendrogram (intermediate) --------------------------------
echo
echo "[4/7] Building clean (un-annotated) dendrogram..."
python3 "${SCRIPTS_DIR}/KS_allchr_dendrogram.py" \
    --assignments "${ASSIGNMENTS}" \
    --bed "${PASS_BED}" \
    --colors "${COLORS_DIR}" \
    --output "${OUT_DIR}/allchr_dendrogram.svg" \
    --matrix-type count_log1p_zscore \
    --sil-threshold 0.5 \
    --centroid-sd 5 \
    --row-height 12 --bar-height 10
echo "    -> ${OUT_DIR}/allchr_dendrogram.svg"

# ---- Step 5: Outlier annotation TSV -----------------------------------------
echo
echo "[5/7] Generating outlier annotations..."
python3 "${SCRIPTS_DIR}/KS_allchr_annotate.py" \
    --svg "${OUT_DIR}/allchr_dendrogram.svg" \
    --assignments "${ASSIGNMENTS}" \
    --bed "${PASS_BED}" \
    --output "${OUT_DIR}/allchr_outlier_annotations.tsv"
echo "    -> ${OUT_DIR}/allchr_outlier_annotations.tsv"

# ---- Step 6: Annotated dendrogram (FIGURE 6A) -------------------------------
echo
echo "[6/7] Building annotated dendrogram (Figure 6A)..."
python3 "${SCRIPTS_DIR}/KS_allchr_dendrogram.py" \
    --assignments "${ASSIGNMENTS}" \
    --bed "${PASS_BED}" \
    --colors "${COLORS_DIR}" \
    --output "${OUT_DIR}/allchr_dendrogram_annotated.svg" \
    --matrix-type count_log1p_zscore \
    --sil-threshold 0.5 \
    --centroid-sd 5 \
    --row-height 12 --bar-height 10 \
    --annotations "${OUT_DIR}/allchr_outlier_annotations.tsv"
echo "    -> ${OUT_DIR}/allchr_dendrogram_annotated.svg"

# ---- Step 7a: Stacked Major/Minor barplot (FIGURE 6C) -----------------------
echo
echo "[7a/7] Stacked Major/Minor barplot (Figure 6C)..."
python3 "${SCRIPTS_DIR}/KS_allchr_barplot.py" \
    --assignments "${ASSIGNMENTS}" \
    --bed "${PASS_BED}" \
    --output "${OUT_DIR}/allchr_barplot.svg" \
    --matrix-type count_log1p_zscore \
    --sil-threshold 0.5 \
    --centroid-sd 5
echo "    -> ${OUT_DIR}/allchr_barplot.svg"

# ---- Step 7b: Allele co-occurrence heatmap (FIGURE 6D) ----------------------
echo
echo "[7b/7] Allele co-occurrence heatmap, chr3/8/11/12 (Figure 6D)..."
python3 "${SCRIPTS_DIR}/KS_allchr_allele_heatmap.py" \
    --assignments "${ASSIGNMENTS}" \
    --bed "${PASS_BED}" \
    --output-dir "${OUT_DIR}" \
    --matrix-type count_log1p_zscore \
    --sil-threshold 0.5 \
    --centroid-sd 5
echo "    -> ${OUT_DIR}/allchr_allele_heatmap.svg"
echo "    -> ${OUT_DIR}/allchr_allele_heatmap_2x2.svg"

echo
echo "=========================================================================="
echo " Finished: $(date '+%F %T')"
echo " All outputs in: ${OUT_DIR}"
echo "=========================================================================="
