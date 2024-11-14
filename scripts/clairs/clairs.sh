#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR_BAM=$1
CONTROL_BAM=$2
OUTPUT_DIR=$3
ASSEMBLY_HAP1=$4
ASSEMBLY_HAP2=$5
KMER_RATIO=$6

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa
samtools faidx ${OUTPUT_DIR}/reference.fa

awk '{if ($4 > 0.8) print $1}' ${KMER_RATIO} > ${OUTPUT_DIR}/read_list.txt
samtools view -@ 16 -Shb ${TUMOR_BAM} -N ${OUTPUT_DIR}/read_list.txt > ${OUTPUT_DIR}/${SAMPLE}_filted.unsorted
samtools sort -@ 16 -o ${OUTPUT_DIR}/${SAMPLE}_filtered.sorted.bam ${OUTPUT_DIR}/${SAMPLE}_filted.unsorted
samtools index -@ 16 ${OUTPUT_DIR}/${SAMPLE}_filtered.sorted.bam
rm ${OUTPUT_DIR}/${SAMPLE}_filted.unsorted

/opt/bin/ClairS/run_clairs \
    -T ${OUTPUT_DIR}/${SAMPLE}_filtered.sorted.bam \
    -N ${CONTROL_BAM} \
    -R ${OUTPUT_DIR}/reference.fa \
    -o ${OUTPUT_DIR} \
    -t 16 \
    -p hifi_sequel2 \
    --include_all_ctgs 

#rm ${OUTPUT_DIR}/${SAMPLE}_filtered.sorted.bam

echo ${?}
