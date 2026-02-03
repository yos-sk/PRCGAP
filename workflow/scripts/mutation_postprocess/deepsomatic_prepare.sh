#!/bin/bash
# DeepSomatic postprocess prepare script
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

DEEPSOMATIC_DIR=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
OUTPUT_DIR=$4
OUTPUT_VCF=$5
OUTPUT_REF=$6

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_REF}
samtools faidx ${OUTPUT_REF}
cp ${DEEPSOMATIC_DIR}/output.vcf.gz ${OUTPUT_VCF}

echo ${?}
