#!/bin/bash
set -xv
set -o errexit
set -o nounset
set -o pipefail

CLAIRS_RESULTS=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
WORK_DIR=$4
OUTPUT_DIR=$5
CHAINFILE=$6

mkdir -p ${WORK_DIR}
mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${WORK_DIR}/reference.fa
samtools faidx ${WORK_DIR}/reference.fa

python3 /tools/mutation_postprocess/scripts/parse_vcf.py \
    ${CLAIRS_RESULTS} > ${WORK_DIR}/output.parsed.bed

coordconv \
    -b ${WORK_DIR}/output.parsed.bed \
    -c ${CHAINFILE} \
> ${WORK_DIR}/output.coordconv_chm13.bed

awk '{if ($9 == "PASS") print}' ${WORK_DIR}/output.coordconv_chm13.bed > ${WORK_DIR}/output.coordconv_chm13_pass.bed
awk '{if ($5 != "Match") print}' ${WORK_DIR}/output.coordconv_chm13_pass.bed > ${WORK_DIR}/output.coordconv_chm13_failed.bed


mutation_postprocess extract-seq \
    -i ${WORK_DIR}/output.coordconv_chm13_failed.bed \
    -f ${WORK_DIR}/reference.fa \
> ${WORK_DIR}/output.coordconv_chm13_failed_seq.txt

mutation_postprocess realignment \
    -i ${WORK_DIR}/output.coordconv_chm13_failed_seq.txt > ${WORK_DIR}/pair_info.txt

mutation_postprocess group \
    -i ${WORK_DIR}/pair_info.txt > ${WORK_DIR}/group_info.txt


python3 /tools/mutation_postprocess/scripts/remove_duplicates.py \
    ${WORK_DIR}/output.coordconv_chm13_pass.bed \
    ${WORK_DIR}/group_info.txt > ${OUTPUT_DIR}/output.coordconv_chm13_rmdup.bed
    
python3 /tools/mutation_postprocess/scripts/filter.py \
    ${OUTPUT_DIR}/output.coordconv_chm13_rmdup.bed \
    ${CLAIRS_RESULTS} \
> ${OUTPUT_DIR}/${CLAIRS_RESULTS%.vcf.gz}.filtered.vcf

bgzip -f ${OUTPUT_DIR}/${CLAIRS_RESULTS%.vcf.gz}.filtered.vcf > ${OUTPUT_DIR}/${CLAIRS_RESULTS%.vcf.gz}.filtered.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${CLAIRS_RESULTS%.vcf.gz}.filtered.vcf.gz
