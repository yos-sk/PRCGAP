#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
ASSEMBLY_HAP1=$3
ASSEMBLY_HAP2=$4
TUMOR_BAM=$5
NORMAL_BAM=$6
REFERENCE=$7
WORK_DIR=$8
OUTPUT_DIR=$9
DNANN_MODEL=${10}
SEX=${11}

mkdir -p ${WORK_DIR}
if [ ${SEX} = "female" ]; then
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' ${REFERENCE} > ${WORK_DIR}/reference_noY.fa
    REFERENCE=${WORK_DIR}/reference_noY.fa
fi

for hap in hap1 hap2
do
    if [ $hap = "hap1" ]; then
        INPUT_FASTA=${ASSEMBLY_HAP1}
    else
        INPUT_FASTA=${ASSEMBLY_HAP2}
    fi
        
    # 1. Mask alpha satelite in contigs
    dna-brnn \
        -Ai ${DNANN_MODEL} \
        -t16 ${INPUT_FASTA} > ${WORK_DIR}/${NORMAL}.${hap}_dna-brnn.bed
    gzip -f ${WORK_DIR}/${NORMAL}.${hap}_dna-brnn.bed
    CN_utils mask \
        -i ${INPUT_FASTA} \
        -b ${WORK_DIR}/${NORMAL}.${hap}_dna-brnn.bed.gz \
        > ${WORK_DIR}/${NORMAL}.${hap}.masked.fa

    # 2. Align ref to contigs
    minimap2 -cx asm5 -t 16 ${WORK_DIR}/${NORMAL}.${hap}.masked.fa ${REFERENCE} > ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.paf
    grep -v 'tp:A:S' ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.paf > ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.rmsec.paf

    # 3. Make correspondence table between contigs and ref
    python3 /tools/CN_utils/scripts/create_correspo_table.py \
        -i ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.rmsec.paf \
        > ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table

    # 4. Calculate depth
    awk '{print $1 "\t" $2 "\t" $3}' ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table > ${WORK_DIR}/${NORMAL}.${hap}.ref.bed
    samtools depth -@ 16 -a -b ${WORK_DIR}/${NORMAL}.${hap}.ref.bed ${tumor_bamfile} -Q 40 > ${WORK_DIR}/${TUMOR}.${hap}.depth
    samtools depth -@ 16 -a -b ${WORK_DIR}/${NORMAL}.${hap}.ref.bed ${normal_bamfile} -Q 40 > ${WORK_DIR}/${NORMAL}.${hap}.depth
    gzip -f ${WORK_DIR}/${TUMOR}.${hap}.depth
    gzip -f ${WORK_DIR}/${NORMAL}.${hap}.depth
    CN_utils copynumber \
        -t ${WORK_DIR}/${TUMOR}.${hap}.depth.gz \
        -c ${WORK_DIR}/${NORMAL}.${hap}.depth.gz \
        -i ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table > ${OUTPUT_DIR}/${TUMOR}.${hap}.copynumber.tsv
done

echo ${?}
