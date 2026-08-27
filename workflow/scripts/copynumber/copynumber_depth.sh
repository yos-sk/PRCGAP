#!/bin/bash
# Per-haplotype depth over the reference-table intervals for tumor and normal,
# reduced to windowed tumor/normal depth ratios.
#
# mosdepth measures the 50 kb windows directly, so neither the per-base
# intermediate nor the python aggregation is needed. --no-per-base is not just
# an output switch: with per-base enabled mosdepth traverses every contig in the
# BAM header (both haplotypes are present in a refined BAM) and writes GB-scale
# output, while with it the traversal is confined to the contigs in --by.
#
# mosdepth reports a region MEAN to two decimals; cbs.R thresholds the absolute
# per-base SUM over the window, so windows_to_tsv.py converts back and keeps the
# five-column copynumber.tsv format. The resulting quantization is 500 per 50 kb
# window; see workspace/mosdepth_bench for the validation against samtools depth.
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
#   $10 THREADS     mosdepth decompression threads

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

WINDOWS="${WORK_DIR}/${NORMAL}.${HAP}.windows.bed"
python3 "${SCRIPT_DIR}"/copynumber/make_windows.py "${REF_TABLE}" > "${WINDOWS}"

for entry in "${TUMOR}:${TUMOR_BAM}" "${NORMAL}:${NORMAL_BAM}"; do
    SAMPLE=${entry%%:*}
    BAM=${entry#*:}
    # Threads only decompress the BAM, so they keep returning less as they go
    # up. Peak memory is flat in the thread count -- it tracks the longest
    # contig instead.
    mosdepth -t "${THREADS}" -Q 40 --by "${WINDOWS}" --no-per-base \
        "${WORK_DIR}/${SAMPLE}.${HAP}" "${BAM}"
done

python3 "${SCRIPT_DIR}"/copynumber/windows_to_tsv.py \
    -t "${WORK_DIR}/${TUMOR}.${HAP}.regions.bed.gz" \
    -n "${WORK_DIR}/${NORMAL}.${HAP}.regions.bed.gz" \
    -r "${REF_TABLE}" > "${OUT_TSV}"

echo "[copynumber_depth] done: ${TUMOR} ${HAP}"
