#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G
#SBATCH -p mjobs,rjobs

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
INPUT_VCF_DIR=$2
HAP1_FASTA=$3
HAP2_FASTA=$4
OUTPUT_DIR=$5
INPUT_BAM=$6
TOOL=$7
IMAGE=$8


mkdir -p ${OUTPUT_DIR}/realign
cat ${HAP1_FASTA} ${HAP2_FASTA} > ${OUTPUT_DIR}/reference.fa
singularity exec ${IMAGE} samtools faidx ${OUTPUT_DIR}/reference.fa

if [ ${TOOL} = "ClairS" ]; then
    cat <(zcat ${INPUT_VCF_DIR}/output.vcf.gz) <(zgrep -v "#" ${INPUT_VCF_DIR}/indel.vcf.gz) | gzip -f -c > ${OUTPUT_DIR}/output.vcf.gz 
else
    cp ${INPUT_VCF_DIR}/output.vcf.gz ${OUTPUT_DIR}
fi

# 1. Realignment by bwa mem 
singularity exec ${IMAGE} python3 ./scripts/mutation_postprocess/parse_vcf.py \
    -i ${OUTPUT_DIR}/output.vcf.gz \
    -f ${OUTPUT_DIR}/reference.fa \
    -o ${OUTPUT_DIR}/realign/parsed_vcf.bed \
    -g ${OUTPUT_DIR}/realign/realign_ref.fasta \
    -t ${TOOL}
   
REALIGN_JOB=$(sbatch -J realign_${SAMPLE} -o ${OUTPUT_DIR}/realign/log/realign.out -e ${OUTPUT_DIR}/realign/log/realign.err ./scripts/mutation_postprocess/realign.sh \
              ${OUTPUT_DIR}/output.vcf.gz \
              ${HAP1_FASTA} \
              ${HAP2_FASTA} \
              ${OUTPUT_DIR} \
              ${TOOL} \
              ${IMAGE})
REALIGN_JOB_ID=$(echo "$REALIGN_JOB" | awk '{print $4}')

# 2. samtools mpileup
mkdir -p ${OUTPUT_DIR}/pileup/workspace
split -n l/16 ${OUTPUT_DIR}/realign/parsed_vcf.bed ${OUTPUT_DIR}/pileup/workspace/input_

find ${OUTPUT_DIR}/pileup/workspace/ -type f -name input_*.pileup | while read f; do 
    rm  ${f}
done

ls ${OUTPUT_DIR}/pileup/workspace/input_* | while read f; do
    echo -e "${f}\t${f}.pileup"
done > ${OUTPUT_DIR}/pileup/workspace/pileup_tasks.txt

PILEUP_JOB_ID=$(sbatch --dependency=afterok:${REALIGN_JOB_ID} -J pileup_${SAMPLE} -o ${OUTPUT_DIR}/pileup/log/pileup.out -e ${OUTPUT_DIR}/pileup/log/pileup.err \
               --parsable -a 1-$(wc -l ${OUTPUT_DIR}/pileup/workspace/pileup_tasks.txt | awk '{print $1}'):1 ./scripts/mutation_postprocess/pileup.sh \
               ${INPUT_BAM} \
               ${OUTPUT_DIR}/pileup/workspace/pileup_tasks.txt \
               ${OUTPUT_DIR}/reference.fa \
               ${IMAGE})

MERGE_JOB=$(sbatch --dependency=afterok:${PILEUP_JOB_ID} -J merge_pileup_${SAMPLE} -o ${OUTPUT_DIR}/pileup/log/merge_pileup.out -e ${OUTPUT_DIR}/pileup/log/merge_pileup.err \
              ./scripts/mutation_postprocess/merge_pileup.sh ${SAMPLE} ${OUTPUT_DIR} ${IMAGE})

MERGE_JOB_ID=$(echo "$MERGE_JOB" | awk '{print $4}')


# 3. haplotyping 
sbatch --dependency=afterok:${MERGE_JOB_ID} -J haplotyping_${SAMPLE} -o ${OUTPUT_DIR}/log/haplotyping.out -e ${OUTPUT_DIR}/log/haplotyping.err ./scripts/mutation_postprocess/haplotyping.sh \
    ${SAMPLE} \
    ${OUTPUT_DIR} \
    ${INPUT_BAM} \
    ${IMAGE}

echo ${?}
