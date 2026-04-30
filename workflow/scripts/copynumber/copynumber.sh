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
HAP1_SATELLITE=${10}
HAP2_SATELLITE=${11}
SEX=${12}
SCRIPT_DIR=${13:-./scripts}
THREAD=${14:-8}

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}
if [ ${SEX} = "female" ]; then
    awk '/^>/ {p = ($0 !~ /^>chrY/)} p' ${REFERENCE} > ${WORK_DIR}/reference_noY.fa
    REFERENCE=${WORK_DIR}/reference_noY.fa
fi

for hap in hap1 hap2
do
    if [ $hap = "hap1" ]; then
        INPUT_FASTA=${ASSEMBLY_HAP1}
        SATELLITE=${HAP1_SATELLITE}
    else
        INPUT_FASTA=${ASSEMBLY_HAP2}
        SATELLITE=${HAP2_SATELLITE}
    fi

    bedtools maskfasta \
        -fi ${INPUT_FASTA} \
        -bed ${SATELLITE} \
        -fo ${WORK_DIR}/${NORMAL}.${hap}.masked.fa

    # 2. Align ref to contigs
    minimap2 -cx asm5 -t ${THREAD} ${WORK_DIR}/${NORMAL}.${hap}.masked.fa ${REFERENCE} > ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.paf
    grep -v 'tp:A:S' ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.paf > ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.rmsec.paf

    # 3. Make correspondence table between contigs and ref
    python3 "${SCRIPT_DIR}"/copynumber/make_reference_table.py \
        -i ${WORK_DIR}/${NORMAL}.${hap}.masked_ref.rmsec.paf \
    > ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table
done

# 4. Postprocess sex chromosomes for male samples
if [ ${SEX} = "male" ]; then
    python3 "${SCRIPT_DIR}"/copynumber/postprocess_sex_chrom.py \
        --hap1 ${OUTPUT_DIR}/${NORMAL}.hap1.ref.table \
        --hap2 ${OUTPUT_DIR}/${NORMAL}.hap2.ref.table \
        --out1 ${OUTPUT_DIR}/${NORMAL}.hap1.ref.table \
        --out2 ${OUTPUT_DIR}/${NORMAL}.hap2.ref.table
fi

# 5. Calculate depth
for hap in hap1 hap2
do
    awk '{print $1 "\t" $2 "\t" $3}' ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table > ${WORK_DIR}/${NORMAL}.${hap}.ref.bed
    samtools depth -@ ${THREAD} -a -b ${WORK_DIR}/${NORMAL}.${hap}.ref.bed ${TUMOR_BAM} -Q 40 | awk '{print $1 "\t" $2 - 1 "\t" $2 "\t" $3}' >  ${WORK_DIR}/${TUMOR}.${hap}.depth.bed
    samtools depth -@ ${THREAD} -a -b ${WORK_DIR}/${NORMAL}.${hap}.ref.bed ${NORMAL_BAM} -Q 40 | awk '{print $1 "\t" $2 - 1 "\t" $2 "\t" $3}' >  ${WORK_DIR}/${NORMAL}.${hap}.depth.bed
    bgzip -f -@ ${THREAD} ${WORK_DIR}/${TUMOR}.${hap}.depth.bed
    tabix -p bed ${WORK_DIR}/${TUMOR}.${hap}.depth.bed.gz
    bgzip -f -@ ${THREAD} ${WORK_DIR}/${NORMAL}.${hap}.depth.bed
    tabix -p bed ${WORK_DIR}/${NORMAL}.${hap}.depth.bed.gz
    python3 "${SCRIPT_DIR}"/copynumber/copynumber_window.py \
        -t ${WORK_DIR}/${TUMOR}.${hap}.depth.bed.gz \
        -n ${WORK_DIR}/${NORMAL}.${hap}.depth.bed.gz \
        -r ${OUTPUT_DIR}/${NORMAL}.${hap}.ref.table > ${OUTPUT_DIR}/${TUMOR}.${hap}.copynumber.tsv
done

echo ${?}
