#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
SIMPLE_REPEAT=$3
BAM_FILE=$4
REFERENCE=${OUTPUT_DIR}/reference.fa

if [ -n $REFERENCE ]; then
    exit(1)
fi


python3 /tools/nanomonsv/misc/add_simple_repeat.py \
    ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.txt \
    ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.txt \
    ${SIMPLE_REPEAT}

if [ -e ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt ]; then
    rm  ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt
fi

head -n 1 ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.txt >> ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt
tail -n +2 ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.txt | grep -v "Simple_repeat" | grep -v "Too_small_size" >> ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt

awk '{if (NR != 1) print $1 "\t" $2 - 1 "\t" $2 "\t" $8 "\n" $4 "\t" $5 - 1 "\t" $5 "\t" $8}' \
   ${OUTPUT_DIR}//${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt \
> ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.bed

nanomonsv_postprocess extract-seq \
    -b ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.bed \
    -f ${REFERENCE} > ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.seq.txt

nanomonsv_postprocess realignment \
    -i ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.seq.txt \
    -s ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.supporting_read.txt \
    -b ${BAM_FILE} \
    -d 98 \
    -l 160 \
    1>${OUTPUT_DIR}/${SAMPLE}.nanomonsv.group_info.txt 2>${OUTPUT_DIR}/${SAMPLE}.nanomonsv.pair_info.txt

nanomonsv_postprocess filt \
    -i ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.group_info.txt \
    -n ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.result.filt.pass_low_vaf.txt \
    -s ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.supporting_read.txt \
    -b ${BAM_FILE} \
    1>${OUTPUT_DIR}/${SAMPLE}.nanomonsv.new_result.txt 2>${OUTPUT_DIR}/${SAMPLE}.nanomonsv_postprocess_filt.log


python3 /tools/nanomonsv_postprocess/script/sv_type.py \
    ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.new_result.txt \
    ${OUTPUT_DIR}/${SAMPLE}.nanomonsv.new_result.sv_typed.txt

echo ${?}
