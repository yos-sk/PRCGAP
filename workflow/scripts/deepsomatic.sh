#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7
THREAD=${8:-16}

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa
samtools faidx ${OUTPUT_DIR}/reference.fa

run_deepsomatic \
    --model_type=PACBIO \
    --ref=${OUTPUT_DIR}/reference.fa \
    --reads_normal=${NORMAL_BAM} \
    --reads_tumor=${TUMOR_BAM} \
    --output_vcf=${OUTPUT_DIR}/output.vcf.gz \
    --sample_name_tumor=${TUMOR} \
    --sample_name_normal=${NORMAL} \
    --num_shards=${THREAD} \
    --logging_dir=${OUTPUT_DIR}/logs \
    --intermediate_results_dir=${OUTPUT_DIR}/intermediate_results_dir 

echo ${?}
