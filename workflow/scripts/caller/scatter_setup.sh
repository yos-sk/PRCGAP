#!/bin/bash
# Shared inputs for a scattered caller run: the concatenated hap1+hap2
# reference and the tumor/normal BAMs the chunk jobs read.
#
# hifiasm HiC assemblies name contigs "sample#hap#contig". The '#' trips
# ClairS's internal Clair3 phasing (longphase emits out-of-order phased VCFs
# and the downstream tabix fails), so names are sanitised to '_HASH_' here and
# restored by merge_caller_vcf.sh. `samtools reheader` rewrites the whole BAM,
# so it must happen once here rather than per chunk. Without '#' contigs the
# BAMs are symlinked instead of copied.

set -xv
set -o errexit
set -o nounset
set -o pipefail

TUMOR_BAM=$1
NORMAL_BAM=$2
ASSEMBLY_HAP1=$3
ASSEMBLY_HAP2=$4
OUTPUT_DIR=$5
THREAD=${6:-8}

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa

if grep -qm1 '^>.*#' ${OUTPUT_DIR}/reference.fa; then
    sed -i '/^>/ s/#/_HASH_/g' ${OUTPUT_DIR}/reference.fa

    for ROLE in tumor normal; do
        if [ ${ROLE} = "tumor" ]; then SRC=${TUMOR_BAM}; else SRC=${NORMAL_BAM}; fi
        samtools view -H ${SRC} | sed '/^@SQ/ s/#/_HASH_/g' > ${OUTPUT_DIR}/${ROLE}.header.sam
        samtools reheader ${OUTPUT_DIR}/${ROLE}.header.sam ${SRC} > ${OUTPUT_DIR}/${ROLE}.bam
        samtools index -@ ${THREAD} ${OUTPUT_DIR}/${ROLE}.bam
        rm -f ${OUTPUT_DIR}/${ROLE}.header.sam
    done

    echo 1 > ${OUTPUT_DIR}/sanitized.flag
else
    # Snakemake runs with cwd=output_dir, so the link targets must be absolute.
    for ROLE in tumor normal; do
        if [ ${ROLE} = "tumor" ]; then SRC=${TUMOR_BAM}; else SRC=${NORMAL_BAM}; fi
        ln -sf "$(realpath ${SRC})" ${OUTPUT_DIR}/${ROLE}.bam
        ln -sf "$(realpath ${SRC}.bai)" ${OUTPUT_DIR}/${ROLE}.bam.bai
    done

    echo 0 > ${OUTPUT_DIR}/sanitized.flag
fi

samtools faidx ${OUTPUT_DIR}/reference.fa

echo ${?}
