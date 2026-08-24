#!/bin/bash
# Pileup processing script
# Used by Snakemake for mutation postprocessing
#
# Phase A change (2026-07):
#   * --no-BAQ defaults to true (BAQ is the dominant memory consumer for
#     long-read / high-depth pileup; disabling it cuts peak RSS ~26x and keeps
#     the read-base column format identical, ref-relative `.,`).
#   * each chunk is restricted to its contig(s) with `samtools mpileup -r`
#     (index seek + early termination) instead of scanning the whole BAM per
#     chunk. Requires the BAM index (.bai) to be present.
#   * optional -d/--max-depth via MAX_DEPTH (0 => samtools default).

set -xv
set -o errexit
set -o nounset
set -o pipefail

INPUT_BAM=$1
PILEUP_TASKS=$2
REFERENCE_FA=$3
SAMPLE=$4
OUTPUT_DIR=$5
THREADS=$6
MEM_MB=${7:-8000}    # Memory (MB) of upper bound  for sorting buffer
NO_BAQ=${8:-true} 
MAX_DEPTH=${9:-0}

mkdir -p ${OUTPUT_DIR}/pileup/workspace

BAQ_OPT=""
if [ "${NO_BAQ}" = "true" ]; then BAQ_OPT="--no-BAQ"; fi
DEPTH_OPT=""
if [ "${MAX_DEPTH}" -gt 0 ] 2>/dev/null; then DEPTH_OPT="-d ${MAX_DEPTH}"; fi
export INPUT_BAM REFERENCE_FA BAQ_OPT DEPTH_OPT

run_pileup() {
    local bed_file=$1
    local pileup_file=$2
    # iterate the distinct contigs in this chunk; -r seeks via the BAM index so
    # only the contig's records are read (no full-BAM scan). -l keeps the exact
    # candidate positions.
    for ctg in $(cut -f1 "${bed_file}" | sort -u); do
        samtools mpileup ${BAQ_OPT} ${DEPTH_OPT} -r "${ctg}" -l "${bed_file}" \
            -f "${REFERENCE_FA}" --output-QNAME "${INPUT_BAM}"
    done | awk '{print $1 "\t" $2 -1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > "${pileup_file}"
}
export -f run_pileup

xargs -a "${PILEUP_TASKS}" -n 2 -P "${THREADS}" bash -c 'run_pileup "$@"' _

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -S $(( MEM_MB / 2 < 8192 ? MEM_MB / 2 : 8192 ))M --parallel="${THREADS}" -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
bgzip -@ ${THREADS} -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz

echo ${?}
