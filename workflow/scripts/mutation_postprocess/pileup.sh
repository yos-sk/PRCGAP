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

mkdir -p ${OUTPUT_DIR}/pileup/workspace

while IFS=$'\t' read bed_file pileup_file; do
    samtools mpileup -l ${bed_file} -f ${REFERENCE_FA} --output-QNAME ${INPUT_BAM} | \
        awk '{print $1 "\t" $2 -1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > ${pileup_file}
done < ${PILEUP_TASKS}

cat ${OUTPUT_DIR}/pileup/workspace/*.pileup | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
bgzip -@ ${THREADS} -f ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed
tabix -p bed ${OUTPUT_DIR}/pileup/${SAMPLE}_pileup.bed.gz

echo ${?}
