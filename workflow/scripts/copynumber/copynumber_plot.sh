#!/bin/bash
# Draw the two-haplotype copy-number plot from the per-haplotype segmentation.
#
# Required positional args:
#   $1  TUMOR           tumor sample name
#   $2  NORMAL          normal sample name
#   $3  OUTPUT_DIR      dir holding the per-haplotype tsv/cbs/ploidy/ref.table
#   $4  WORK_DIR        scratch dir
#   $5  CHM13_LENGTHS   chromosome-length table (copynumber_reference.sh)
#   $6  CENSAT_BED      cenSat BED in contig coords, or empty to build one
#   $7  HAP1_SATELLITE  hap1 satellite BED.gz (fallback annotation source)
#   $8  HAP2_SATELLITE  hap2 satellite BED.gz
#   $9  BIN_WIDTH       bin width
#   $10 HAP1_LABEL      hap1 plot label
#   $11 HAP2_LABEL      hap2 plot label
#   $12 PLOT_SEX_CHROM  true | false
#   $13 CHM13_CENSAT    CHM13 cenSat BED(.gz), or empty
#   $14 SCRIPT_DIR      absolute path to workflow/scripts

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
OUTPUT_DIR=$3
WORK_DIR=$4
CHM13_LENGTHS=$5
CENSAT_BED=${6:-}
HAP1_SATELLITE=$7
HAP2_SATELLITE=$8
BIN_WIDTH=${9:-0.05}
HAP1_LABEL=${10:-Haplotype1}
HAP2_LABEL=${11:-Haplotype2}
PLOT_SEX_CHROM=${12:-true}
CHM13_CENSAT=${13:-}
SCRIPT_DIR=${14}

mkdir -p "${WORK_DIR}"

PLOIDY_HAP1=$(< "${OUTPUT_DIR}/${TUMOR}.hap1.ploidy")
PLOIDY_HAP2=$(< "${OUTPUT_DIR}/${TUMOR}.hap2.ploidy")

# Centromere/satellite annotation in contig coordinates. Prefer the cenSat BED
# when provided; otherwise build one from the per-haplotype dna-brnn satellite
# BEDs used for masking (decompress both and recompress as one bgzip stream).
if [ -n "${CENSAT_BED}" ]; then
    ANNOTATION=${CENSAT_BED}
else
    zcat "${HAP1_SATELLITE}" "${HAP2_SATELLITE}" \
        | bgzip -c > "${WORK_DIR}/${NORMAL}.satellite.annotation.bed.gz"
    ANNOTATION=${WORK_DIR}/${NORMAL}.satellite.annotation.bed.gz
fi

Rscript "${SCRIPT_DIR}"/copynumber/plot_copy_number.R \
    -i "${OUTPUT_DIR}/${TUMOR}.hap1.copynumber.tsv" \
    -j "${OUTPUT_DIR}/${TUMOR}.hap2.copynumber.tsv" \
    -k "${OUTPUT_DIR}/${TUMOR}.hap1.cbs.split.txt" \
    -l "${OUTPUT_DIR}/${TUMOR}.hap2.cbs.split.txt" \
    -o "${OUTPUT_DIR}/${TUMOR}.copynumber.png" \
    -s "${TUMOR}" \
    -p "${PLOIDY_HAP1}" \
    -t "${PLOIDY_HAP2}" \
    -q "${OUTPUT_DIR}/${NORMAL}.hap1.ref.table" \
    -r "${OUTPUT_DIR}/${NORMAL}.hap2.ref.table" \
    -a "${ANNOTATION}" \
    -w "${BIN_WIDTH}" \
    --hap1_label "${HAP1_LABEL}" \
    --hap2_label "${HAP2_LABEL}" \
    --plot-sex-chrom "${PLOT_SEX_CHROM}" \
    --chm13_censat "${CHM13_CENSAT}" \
    --chm13_lengths "${CHM13_LENGTHS}"

echo "[copynumber_plot] done: ${TUMOR}"
