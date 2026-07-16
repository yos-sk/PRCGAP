#!/bin/bash
# ClairS postprocess prepare script
# Used by Snakemake for mutation postprocessing

set -xv
set -o errexit
set -o nounset
set -o pipefail

CLAIRS_DIR=$1
ASSEMBLY_HAP1=$2
ASSEMBLY_HAP2=$3
OUTPUT_DIR=$4
OUTPUT_VCF=$5
OUTPUT_REF=$6

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_REF}
samtools faidx ${OUTPUT_REF}

# Combine ClairS VCF (if indel.vcf.gz exists). Strip only the VCF header
# (lines starting with '#'); anchor to '^#' so records on contigs whose names
# contain '#' (e.g. hifiasm HiC 'sample#hap#contig') are NOT dropped.
if [ -f ${CLAIRS_DIR}/indel.vcf.gz ]; then
    zcat ${CLAIRS_DIR}/output.vcf.gz | gzip -c > ${OUTPUT_VCF}
    zgrep -v "^#" ${CLAIRS_DIR}/indel.vcf.gz | gzip -f -c >> ${OUTPUT_VCF}
else
    cp ${CLAIRS_DIR}/output.vcf.gz ${OUTPUT_VCF}
fi

echo ${?}
