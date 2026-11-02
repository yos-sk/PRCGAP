#!/bin/bash
# Per-haplotype depth over the reference-table intervals for tumor and normal,
# reduced to windowed tumor/normal depth ratios.
#
# samtools depth writes per-base depth (-a: zero-depth positions included) over
# the reference-table intervals; copynumber_window.py then sums it into 50 kb
# windows and writes the five-column copynumber.tsv that cbs.R reads.
#
# Required positional args:
#   $1  TUMOR       tumor sample name
#   $2  NORMAL      normal sample name
#   $3  HAP         hap1 | hap2
#   $4  TUMOR_BAM   refined tumor BAM
#   $5  NORMAL_BAM  refined normal BAM
#   $6  REF_TABLE   published reference table for this haplotype
#   $7  WORK_DIR    scratch dir
#   $8  OUT_TSV     copynumber TSV to write
#   $9  SCRIPT_DIR  absolute path to workflow/scripts
#   $10 THREADS     samtools / bgzip threads

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
HAP=$3
TUMOR_BAM=$4
NORMAL_BAM=$5
REF_TABLE=$6
WORK_DIR=$7
OUT_TSV=$8
SCRIPT_DIR=$9
THREADS=${10:-4}

mkdir -p "${WORK_DIR}" "$(dirname "${OUT_TSV}")"

REF_BED="${WORK_DIR}/${NORMAL}.${HAP}.ref.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' "${REF_TABLE}" > "${REF_BED}"

for entry in "${TUMOR}:${TUMOR_BAM}" "${NORMAL}:${NORMAL_BAM}"; do
    SAMPLE=${entry%%:*}
    BAM=${entry#*:}
    DEPTH="${WORK_DIR}/${SAMPLE}.${HAP}.depth.bed"
    samtools depth -@ "${THREADS}" -a -b "${REF_BED}" -Q 40 "${BAM}" \
        | awk 'BEGIN{OFS="\t"} {print $1, $2 - 1, $2, $3}' > "${DEPTH}"
    bgzip -f -@ "${THREADS}" "${DEPTH}"
    tabix -f -p bed "${DEPTH}.gz"
done

python3 "${SCRIPT_DIR}"/copynumber/copynumber_window.py \
    -t "${WORK_DIR}/${TUMOR}.${HAP}.depth.bed.gz" \
    -n "${WORK_DIR}/${NORMAL}.${HAP}.depth.bed.gz" \
    -r "${REF_TABLE}" > "${OUT_TSV}"

echo "[copynumber_depth] done: ${TUMOR} ${HAP}"
