#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
hap1_contig=$2
hap2_contig=$3
FASTQ=$4
OPTION_SPLIT=$5
WORK_DIR=$6
OUTPUT_DIR=$7
THREAD=$8
DATA=$9

# Step1: Mapping
mkdir -p ${WORK_DIR} 
OUTPUT_BAM_PREFIX=${WORK_DIR}/${SAMPLE}

cat ${hap1_contig} ${hap2_contig} > ${WORK_DIR}/reference.fa
if [ ${DATA} == "hifi" ]; then
    minimap2 -t ${THREAD} -ax asm5 --MD ${WORK_DIR}/reference.fa ${FASTQ} | samtools view -Shb > ${OUTPUT_BAM_PREFIX}.unsorted
else
    minimap2 -t ${THREAD} -ax asm10 --MD ${WORK_DIR}/reference.fa ${FASTQ} | samtools view -Shb > ${OUTPUT_BAM_PREFIX}.unsorted
fi
samtools sort -@ ${THREAD} -m 2G -n ${OUTPUT_BAM_PREFIX}.unsorted -o ${OUTPUT_BAM_PREFIX}.bam
rm ${OUTPUT_BAM_PREFIX}.unsorted

# Step2: Extract haplotype-specific unique k-mer
mkdir -p ${WORK_DIR}/meryl
meryl count k=21 threads=${THREAD} ${hap1_contig} output ${WORK_DIR}/meryl/hap1.meryl
meryl count k=21 threads=${THREAD} ${hap2_contig} output ${WORK_DIR}/meryl/hap2.meryl

meryl difference ${WORK_DIR}/meryl/hap1.meryl ${WORK_DIR}/meryl/hap2.meryl output ${WORK_DIR}/meryl/hap1.uniq.meryl
meryl difference ${WORK_DIR}/meryl/hap2.meryl ${WORK_DIR}/meryl/hap1.meryl output ${WORK_DIR}/meryl/hap2.uniq.meryl

meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap1.uniq.meryl > ${WORK_DIR}/meryl/hap1.uniq.tsv
meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap2.uniq.meryl > ${WORK_DIR}/meryl/hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 == 1) print}' ${WORK_DIR}/meryl/${hap}.uniq.tsv > ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.uniq.tsv
done


mkdir -p ${OUTPUT_DIR}
bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap1.cnt.uniq.tsv.gz \
    -f ${hap1_contig} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz

bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap2.cnt.uniq.tsv.gz \
    -f ${hap2_contig} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz

# Step 3; List up contig names
grep ">" ${hap1_contig} | sed s/\>// > ${OUTPUT_DIR}/hap1_list.txt
gzip -f ${OUTPUT_DIR}/hap1_list.txt
grep ">" ${hap2_contig} | sed s/\>// > ${OUTPUT_DIR}/hap2_list.txt
gzip -f ${OUTPUT_DIR}/hap2_list.txt

# Step 4: Refine BAM file
if [ $OPTION_SPLIT = "true" ]
then
    mkdir -p ${WORK_DIR}/split
    SIZE=`split_bam size --input-file ${OUTPUT_BAM_PREFIX}.bam`
    split_bam split \
        --input-file ${OUTPUT_BAM_PREFIX}.bam \
        --output-dir ${WORK_DIR}/split \
        --input-size ${SIZE} \
        --num-split ${THREAD}
   

    for i in $(seq 0 $(( ${THREAD} - 1))); do
        bam_refiner refine \
            --input-bam ${WORK_DIR}/split/${i}.bam \
            --output-bam ${WORK_DIR}/split/${i}.refined.bam \
            --hap1-tabix ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz \
            --hap2-tabix ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz \
            --hap1-list ${OUTPUT_DIR}/hap1_list.txt.gz \
            --hap2-list ${OUTPUT_DIR}/hap2_list.txt.gz \
            --kmer-size 21 \
            1>${WORK_DIR}/split/${i}.bam_refiner.tsv 2>${WORK_DIR}/split/${i}.bam_refiner.log &
    done
    wait
   
    cat ${WORK_DIR}/split/*.bam_refiner.tsv > ${OUTPUT_DIR}/bam_refiner_result.tsv
    cat ${WORK_DIR}/split/*.bam_refiner.log > ${OUTPUT_DIR}/bam_refiner.log
    samtools merge \
        -@ ${THREAD} \
        -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
        ${WORK_DIR}/split/*.refined.bam
else
    bam_refiner \
        --input-bam ${OUTPUT_BAM_PREFIX}.bam \
        --output-bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
        --hap1-tabix ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz \
        --hap2-tabix ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz \
        --hap1-list ${OUTPUT_DIR}/hap1_list.txt.gz \
        --hap2-list ${OUTPUT_DIR}/hap2_list.txt.gz \
        --kmer-size 21 \
        1>${OUTPUT_DIR}/bam_refiner_result.tsv 2>${OUTPUT_DIR}/bam_refiner.log
fi

samtools sort \
    -@ ${THREAD} \
    -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam 
samtools index ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam 
rm ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam

gzip -f ${OUTPUT_DIR}/bam_refiner_result.tsv
gzip -f ${OUTPUT_DIR}/bam_refiner.log

rm -rf ${WORK_DIR}


bam_refiner kmer-ratio \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    --threads ${THREAD} \
    > ${OUTPUT_DIR}/${SAMPLE}_kmer_ratio.txt

echo ${?}
