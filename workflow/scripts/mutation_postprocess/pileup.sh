#!/bin/bash
# Pileup processing script
# Used by Snakemake for mutation postprocessing

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
MEM_MB=${7:-8000}   # 割当メモリ(MB)。sort のバッファ上限に使う
NO_BAQ=${8:-false}  # true なら samtools mpileup に --no-BAQ を付与(BAQ計算を省きメモリ/計算を削減)

mkdir -p ${OUTPUT_DIR}/pileup/workspace

BAQ_OPT=""
if [ "${NO_BAQ}" = "true" ]; then BAQ_OPT="--no-BAQ"; fi
export INPUT_BAM REFERENCE_FA BAQ_OPT

run_pileup() {
    local bed_file=$1
    local pileup_file=$2
    samtools mpileup ${BAQ_OPT} -l "${bed_file}" -f "${REFERENCE_FA}" --output-QNAME "${INPUT_BAM}" | \
        awk '{print $1 "\t" $2 -1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > "${pileup_file}"
}
export -f run_pileup

xargs -a "${PILEUP_TASKS}" -n 2 -P "${THREADS}" bash -c 'run_pileup "$@"' _

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -S $(( MEM_MB / 2 < 8192 ? MEM_MB / 2 : 8192 ))M --parallel="${THREADS}" -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
bgzip -@ ${THREADS} -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz

echo ${?}
