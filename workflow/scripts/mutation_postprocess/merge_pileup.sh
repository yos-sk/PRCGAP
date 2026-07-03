#!/bin/bash
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
THREADS=$3
MEM_MB=${4:-8000}   # 割当メモリ(MB)。sort のバッファ上限に使う

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -S $(( MEM_MB / 2 < 8192 ? MEM_MB / 2 : 8192 ))M --parallel="${THREADS}" -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
bgzip -@ ${THREADS} -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz

echo ${?}