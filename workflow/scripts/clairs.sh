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
THREAD=${6:-16}
# ClairS platform model. One of: hifi_sequel2, hifi_revio,
# ont_r10_dorado_sup_5khz_ssrs, ont_r10_dorado_sup_4khz.
MODEL=${7:-hifi_sequel2}

mkdir -p ${OUTPUT_DIR}
cat ${ASSEMBLY_HAP1} ${ASSEMBLY_HAP2} > ${OUTPUT_DIR}/reference.fa

# hifiasm HiC assemblies name contigs "sample#hap#contig". The '#' trips ClairS's
# internal Clair3 phasing: longphase emits the per-contig phased VCFs with
# out-of-order positions, so the downstream `tabix` fails
# ("tbx_index_build failed") and ClairS STEP 4 aborts. Work around it by running
# ClairS against sanitized contig names ('#' -> '_HASH_') and restoring the
# original names in the output VCFs afterwards. '_HASH_' is collision-free for
# these assemblies; the whole block is a no-op for assemblies without '#'.
#   - reference: sed the FASTA headers.
#   - BAMs: `samtools reheader` only — reads reference sequences by index, so the
#     read records are untouched and no realignment is needed.
if grep -qm1 '^>.*#' ${OUTPUT_DIR}/reference.fa; then
    sed '/^>/ s/#/_HASH_/g' ${OUTPUT_DIR}/reference.fa > ${OUTPUT_DIR}/reference.sanitized.fa
    mv ${OUTPUT_DIR}/reference.sanitized.fa ${OUTPUT_DIR}/reference.fa

    samtools view -H ${TUMOR_BAM} | sed '/^@SQ/ s/#/_HASH_/g' > ${OUTPUT_DIR}/tumor.header.sam
    samtools reheader ${OUTPUT_DIR}/tumor.header.sam ${TUMOR_BAM} > ${OUTPUT_DIR}/tumor.sanitized.bam
    samtools index -@ ${THREAD} ${OUTPUT_DIR}/tumor.sanitized.bam
    TUMOR_USE=${OUTPUT_DIR}/tumor.sanitized.bam

    samtools view -H ${CONTROL_BAM} | sed '/^@SQ/ s/#/_HASH_/g' > ${OUTPUT_DIR}/normal.header.sam
    samtools reheader ${OUTPUT_DIR}/normal.header.sam ${CONTROL_BAM} > ${OUTPUT_DIR}/normal.sanitized.bam
    samtools index -@ ${THREAD} ${OUTPUT_DIR}/normal.sanitized.bam
    NORMAL_USE=${OUTPUT_DIR}/normal.sanitized.bam

    SANITIZED=1
else
    TUMOR_USE=${TUMOR_BAM}
    NORMAL_USE=${CONTROL_BAM}
    SANITIZED=0
fi

samtools faidx ${OUTPUT_DIR}/reference.fa

/opt/bin/run_clairs \
    -T ${TUMOR_USE} \
    -N ${NORMAL_USE} \
    -R ${OUTPUT_DIR}/reference.fa \
    -o ${OUTPUT_DIR} \
    -t ${THREAD} \
    -p ${MODEL} \
    --include_all_ctgs \
    --enable_indel_calling

# Restore the original '#' contig names in the ClairS output VCFs consumed by
# clairs_prepare.sh (output.vcf.gz, indel.vcf.gz), so coordinates match the
# '#'-named assembly/BAM used downstream. Then drop the sanitized intermediates.
if [ ${SANITIZED} -eq 1 ]; then
    for f in output indel; do
        if [ -f ${OUTPUT_DIR}/${f}.vcf.gz ]; then
            zcat ${OUTPUT_DIR}/${f}.vcf.gz | sed 's/_HASH_/#/g' | bgzip -c > ${OUTPUT_DIR}/${f}.vcf.gz.renamed
            mv ${OUTPUT_DIR}/${f}.vcf.gz.renamed ${OUTPUT_DIR}/${f}.vcf.gz
            tabix -f -p vcf ${OUTPUT_DIR}/${f}.vcf.gz
        fi
    done
    rm -f ${OUTPUT_DIR}/tumor.sanitized.bam ${OUTPUT_DIR}/tumor.sanitized.bam.bai \
          ${OUTPUT_DIR}/normal.sanitized.bam ${OUTPUT_DIR}/normal.sanitized.bam.bai \
          ${OUTPUT_DIR}/tumor.header.sam ${OUTPUT_DIR}/normal.header.sam
fi

echo ${?}
