#!/bin/bash

set -ex

CONF=$1
source ${CONF}


# 1. bam_refiner_hifi
for SAMPLE in ${TUMOR} ${CONTROL}; do
    if [ ${SAMPLE} = ${TUMOR} ]; then
        INPUT_FASTQ=${TUMOR_HIFI}
    else
        INPUT_FASTQ=${CONTROL_HIFI}
    fi
    sbatch -J ${SAMPLE}_bam_refiner_hifi -e log/${SAMPLE}_bam_refiner_hifi.err -o log/${SAMPLE}_bam_refiner_hifi.out \
        scripts/bam_refiner/run_bam_refiner.sh \
        ${SAMPLE} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${INPUT_FASTQ} \
        ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/workspace \
        ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output \
        hifi
done

# 2. bam_refiner_hifi for methylation (optional) 
if [ ${TUMOR_METHYLATION_HIFI} != ${TUMOR_HIFI} ]; then
    sbatch -J ${TUMOR}_bam_refiner_methylation_hifi -e log/${TUMOR}_bam_refiner_methylation_hifi.err -o log/${TUMOR}_bam_refiner_methylation_hifi.out \
        scripts/bam_refiner/run_bam_refiner.sh \
        ${TUMOR} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${TUMOR_METHYLATION_HIFI}\
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/methylation_hifi/workspace \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/methylation_hifi/output \
        hifi
fi    

if [ ${CONTROL_METHYLATION_HIFI} != ${CONTROL_HIFI} ]; then
    sbatch -J ${CONTROL}_bam_refiner_methylation_hifi -e log/${CONTROL}_bam_refiner_methylation_hifi.err -o log/${CONTROL}_bam_refiner_methylation_hifi.out \
        scripts/bam_refiner/run_bam_refiner.sh \
        ${CONTROL} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${CONTROL_METHYLATION_HIFI} \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/methylation_hifi/workspace \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/methylation_hifi/output \
        hifi
fi    

# 3. bam_refiner_ont
for SAMPLE in ${TUMOR} ${CONTROL}; do
    if [ ${SAMPLE} = ${TUMOR} ]; then
        INPUT_FASTQ=${TUMOR_ONT}
    else
        INPUT_FASTQ=${CONTROL_ONT}
    fi
    sbatch -J ${SAMPLE}_bam_refiner_ont -e log/${SAMPLE}_bam_refiner_ont.err -o log/${SAMPLE}_bam_refiner_ont.out \
        scripts/bam_refiner/run_bam_refiner.sh \
        ${SAMPLE} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${INPUT_FASTQ} \
        ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/workspace \
        ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output \
        ont
done


# 4. copynumber
TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_bam_refiner_hifi) | cut -d ' ' -f 2) 
NORMAL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_bam_refiner_hifi) | cut -d ' ' -f 2)
#sbatch -J ${TUMOR}_copynumber -e log/${TUMOR}_copynumber.err -o log/${TUMOR}_copynumber.out \
sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${NORMAL_BAM_REFINER_JOBID} -J ${TUMOR}_copynumber -e log/${TUMOR}_copynumber.err -o log/${TUMOR}_copynumber.out \
    scripts/copynumber/run_copynumber.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
        ${REFERENCE} \
        ${OUTPUT_PREFIX}/copynumber/${TUMOR}/workspace \
        ${OUTPUT_PREFIX}/copynumber/${TUMOR}/output \
        ${SEX}


# 5. nanomonsv parse hifi
for SAMPLE in ${TUMOR} ${CONTROL}; do
    BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_hifi) | cut -d ' ' -f 2)
    sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_nanomonsv_parse_hifi -e log/${SAMPLE}_nanomonsv_parse_hifi.err -o log/${SAMPLE}_nanomonsv_parse_hifi.out \
        scripts/nanomonsv/run_nanomonsv_parse.sh \
            ${SAMPLE} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam \
            ${OUTPUT_PREFIX}/nanomonsv/hifi
done
        
# 6. nanomonsv get hifi
TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_nanomonsv_get_hifi -e log/${TUMOR}_nanomonsv_get_hifi.err -o log/${TUMOR}_nanomonsv_get_hifi.out \
    scripts/nanomonsv/run_nanomonsv_get.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/nanomonsv/hifi \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        hifi \
        ${SIMPLE_REPEAT} 

# 7. nanomonsv postprocess hifi
NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_get_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_nanomonsv_postprocess_hifi -e ./log/${TUMOR}_nanomonsv_postprocess_hifi.err -o ./log/${TUMOR}_nanomonsv_postprocess_hifi.out \
    scripts/nanomonsv/run_nanomonsv_postprocess.sh \
        ${TUMOR} \
        ${OUTPUT_PREFIX}/nanomonsv/hifi \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam

# 8. nanomonsv connect hifi
NANOMONSV_POSTPROCESS_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_postprocess_hifi) | cut -d ' ' -f 2)
#sbatch -J ${TUMOR}_nanomonsv_connect_hifi -e ./log/${TUMOR}_nanomonsv_connect_hifi.err -o ./log/${TUMOR}_nanomonsv_connect_hifi.out \
sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_JOBID} -J ${TUMOR}_nanomonsv_connect_hifi -e ./log/${TUMOR}_nanomonsv_connect_hifi.err -o ./log/${TUMOR}_nanomonsv_connect_hifi.out \
    scripts/nanomonsv/run_nanomonsv_connect.sh \
        ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
        ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.supporting_read.txt \
        ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR} 

# 9. nanomonsv parse ont
for SAMPLE in ${TUMOR} ${CONTROL}; do
    BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_ont) | cut -d ' ' -f 2)
    sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_nanomonsv_parse_ont -e log/${SAMPLE}_nanomonsv_parse_ont.err -o log/${SAMPLE}_nanomonsv_parse_ont.out \
        scripts/nanomonsv/run_nanomonsv_parse.sh \
            ${SAMPLE} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
            ${OUTPUT_PREFIX}/nanomonsv/ont
done

# 10. nanomonsv get ont
TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_nanomonsv_get_ont -e log/${TUMOR}_nanomonsv_get_ont.err -o log/${TUMOR}_nanomonsv_get_ont.out \
    scripts/nanomonsv/run_nanomonsv_get.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/ont/output/${CONTROL}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/nanomonsv/ont \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ont \
        ${SIMPLE_REPEAT} 

# 11. nanomonsv postprocess ont
NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_get_ont) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_nanomonsv_postprocess_ont -e ./log/${TUMOR}_nanomonsv_postprocess_ont.err -o ./log/${TUMOR}_nanomonsv_postprocess_ont.out \
    scripts/nanomonsv/run_nanomonsv_postprocess.sh \
        ${TUMOR} \
        ${OUTPUT_PREFIX}/nanomonsv/ont \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.sorted.bam

# 12. nanomonsv connect ont
NANOMONSV_POSTPROCESS_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_postprocess_ont) | cut -d ' ' -f 2)
#sbatch -J ${TUMOR}_nanomonsv_connect_ont -e ./log/${TUMOR}_nanomonsv_connect_ont.err -o ./log/${TUMOR}_nanomonsv_connect_ont.out \
sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_JOBID} -J ${TUMOR}_nanomonsv_connect_ont -e ./log/${TUMOR}_nanomonsv_connect_ont.err -o ./log/${TUMOR}_nanomonsv_connect_ont.out \
    scripts/nanomonsv/run_nanomonsv_connect.sh \
        ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
        ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.supporting_read.txt \
        ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR} 

# 13. merge nanomonsv results
NANOMONSV_POSTPROCESS_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_postprocess_hifi) | cut -d ' ' -f 2)
NANOMONSV_POSTPROCESS_ONT_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_postprocess_ont) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_HIFI_JOBID},${NANOMONSV_POSTPROCESS_ONT_JOBID} -J ${TUMOR}_nanomonsv_merge -e log/${TUMOR}_nanomonsv_merge_ont.err -o log/${TUMOR}_nanomonsv_merge_ont.out \
    scripts/nanomonsv/run_merge.sh \
        ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
        ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
        ${OUTPUT_PREFIX}/nanomonsv/${TUMOR}.nanomonsv.result.merged.txt \


# 14. clairs
TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_bam_refiner_hifi) | cut -d ' ' -f 2)
CONTROL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_bam_refiner_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${CONTROL_BAM_REFINER_JOBID} -J ${TUMOR}_clairs -e ./log/${TUMOR}_clairs.err -o ./log/${TUMOR}_clairs.out \
    scripts/clairs/run_clairs.sh \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${OUTPUT_PREFIX}/clairs/${TUMOR}


# TODO: write filter.py
# 15. clairs postprocess
CLAIRS_JOB_ID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_clairs) | cut -d ' ' -f 2)
#sbatch -J ${TUMOR}_clairs_postprocess -e ./log/${TUMOR}_clairs_postprocess.err -o ./log/${TUMOR}_clairs_postprocess.out \
sbatch --dependency=afterok:${CLAIRS_JOB_ID} -J ${TUMOR}_clairs_postprocess -e ./log/${TUMOR}_clairs_postprocess.err -o ./log/${TUMOR}_clairs_postprocess.out \
    scripts/clairs/run_postprocess.sh \
        ${OUTPUT_PREFIX}/clairs/${TUMOR}/output.vcf.gz \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${OUTPUT_PREFIX}/clairs/${TUMOR}/postprocess/workspace \
        ${OUTPUT_PREFIX}/clairs/${TUMOR}/postprocess/ouptut \
        ${CHAIN_FILE}


# 16. methylation:
for SAMPLE in ${TUMOR} ${CONTROL}; do
    if [ ${SAMPLE} = ${TUMOR} ]; then
        INPUT_MODBAM_HIFI=${TUMOR_MODBAM_HIFI}
        INPUT_MODBAM_ONT=${TUMOR_MODBAM_ONT}
    else
        INPUT_MODBAM_HIFI=${CONTROL_MODBAM_HIFI}
        INPUT_MODBAM_ONT=${CONTROL_MODBAM_ONT}
    fi

    if [ ${TUMOR_HIFI} != ${TUMOR_METHYLATION_HIFI} ]; then
        BAM_REFINER_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_hifi) | cut -d ' ' -f 2)
        INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/methylation_hifi/output/${SAMPLE}_bam_refined.sorted.bam
    else
        BAM_REFINER_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_methylation_hifi) | cut -d ' ' -f 2)
        INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam
    fi

    #sbatch -J ${SAMPLE}_methylation_hifi -e log/${SAMPLE}_methylation_hifi.err -o log/${SAMPLE}_methylation_hifi.out \
    sbatch --dependency=afterok:${BAM_REFINER_HIFI_JOBID} -J ${SAMPLE}_methylation_hifi -e log/${SAMPLE}_methylation_hifi.err -o log/${SAMPLE}_methylation_hifi.out \
        scripts/methylation/run_methylation.sh \
            ${SAMPLE} \
            ${INPUT_BAM} \
            ${INPUT_MODBAM_HIFI} \
            ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/workspace \
            ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/output \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            HiFi

    BAM_REFINER_ONT_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_ont) | cut -d ' ' -f 2)
    sbatch --dependency=afterok:${BAM_REFINER_ONT_JOBID} -J ${SAMPLE}_methylation_ont -e log/${SAMPLE}_methylation_ont.err -o log/${SAMPLE}_methylation_ont.out \
        scripts/methylation/run_methylation.sh \
            ${SAMPLE} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
            ${INPUT_MODBAM_ONT} \
            ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/workspace \
            ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/output \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            ONT
done

