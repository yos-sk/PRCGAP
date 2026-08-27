#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

# Cap glibc malloc arenas. On many-core hosts (e.g. 192 cores) the default
# (cores * 8) arenas reserve hundreds of GB of virtual memory upfront, which
# trips SGE's s_vmem limit even though resident memory is small.
export MALLOC_ARENA_MAX=4

TUMOR=$1
NORMAL=$2
TUMOR_BAM=$3
NORMAL_BAM=$4
OUTPUT_DIR=$5
ASSEMBLY_HAP1=$6
ASSEMBLY_HAP2=$7
THREAD=${8:-16}
# DeepSomatic model type: PACBIO for HiFi, ONT for ONT.
MODEL_TYPE=${9:-PACBIO}

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa
samtools faidx ${OUTPUT_DIR}/reference.fa

run_deepsomatic \
    --model_type=${MODEL_TYPE} \
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
