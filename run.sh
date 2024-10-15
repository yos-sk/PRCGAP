#!/bin/bash

source ./conf.sh


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

# 2. bam_refiner_ont
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

<<_
# 3. copynumber
TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_bam_refiner_hifi) | cut -d ' ' -f 2)
NORMAL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_bam_refiner_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${NORMAL_BAM_REFINER_JOBID} -J ${TUMOR}_copynumber -e log/${TUMOR}_copynumber.err -o log/${TUMOR}_copynumber.out \
    scripts/copynumber/run_copynumber.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.bam \
        ${CHM13} \
        ${OUTPUT_PREFIX}/copynumber/${TUMOR}/workspace \
        ${OUTPUT_PREFIX}/copynumber/${TUMOR}/output
_

# 4. nanomonsv parse hifi
for SAMPLE in ${TUMOR} ${CONTROL}; do
    BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_hifi) | cut -d ' ' -f 2)
    sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_nanomonsv_parse_hifi -e log/${SAMPLE}_nanomonsv_parse_hifi.err -o log/${SAMPLE}_nanomonsv_parse_hifi.out \
        scripts/nanomonsv/run_nanomonsv_parse.sh \
            ${SAMPLE} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.bam \
            ${OUTPUT_PREFIX}/nanomonsv/hifi
done

        
# 5. nanomonsv get hifi
TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_nanomonsv_get_hifi -e log/${TUMOR}_nanomonsv_get_hifi.err -o log/${TUMOR}_nanomonsv_get_hifi.out \
    scripts/nanomonsv/run_nanomonsv_get.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.bam \
        ${OUTPUT_PREFIX}/nanomonsv/hifi \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        hifi

# 6. nanomonsv postprocess hifi
NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_get_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_nanomonsv_postprocess_hifi -e ./log/${TUMOR}_nanomonsv_postprocess_hifi.err -o ./log/${TUMOR}_nanomonsv_postprocess_hifi.out \
    scripts/nanomonsv/run_nanomonsv_postprocess.sh \
        ${TUMOR} \
        ${OUTPUT_PREFIX}/nanomonsv/hifi \
        ${SIMPLE_REPEAT} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.bam
    
# 7. nanomonsv parse ont
for SAMPLE in ${TUMOR} ${CONTROL}; do
    BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_bam_refiner_ont) | cut -d ' ' -f 2)
    sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_nanomonsv_parse_ont -e log/${SAMPLE}_nanomonsv_parse_ont.err -o log/${SAMPLE}_nanomonsv_parse_ont.out \
        scripts/nanomonsv/run_nanomonsv_parse.sh \
            ${SAMPLE} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.bam \
            ${OUTPUT_PREFIX}/nanomonsv/ont
done
        
# 8. nanomonsv get ont
TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_nanomonsv_get_ont -e log/${TUMOR}_nanomonsv_get_ont.err -o log/${TUMOR}_nanomonsv_get_ont.out \
    scripts/nanomonsv/run_nanomonsv_get.sh \
        ${TUMOR} \
        ${CONTROL} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/ont/output/${CONTROL}_bam_refined.bam \
        ${OUTPUT_PREFIX}/nanomonsv/ont \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ont

# 9. nanomonsv postprocess ont
NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_nanomonsv_get_ont) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_nanomonsv_postprocess_ont -e ./log/${TUMOR}_nanomonsv_postprocess_ont.err -o ./log/${TUMOR}_nanomonsv_postprocess_ont.out \
    scripts/nanomonsv/run_nanomonsv_postprocess.sh \
        ${TUMOR} \
        ${OUTPUT_PREFIX}/nanomonsv/ont \
        ${SIMPLE_REPEAT} \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.bam

# 10. merge nanomonsv results

# 11. clairs
TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_bam_refiner_hifi) | cut -d ' ' -f 2)
CONTROL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_bam_refiner_hifi) | cut -d ' ' -f 2)
sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${CONTROL_BAM_REFINER_JOBID} -J ${TUMOR}_clairs -e ./log/${TUMOR}_clairs.err -o ./log/${TUMOR}_clairs.out \
    scripts/clairs/run_clairs.sh \
        ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.bam \
        ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.bam \
        ${ASSEMBLY_HAP1} \
        ${ASSEMBLY_HAP2} \
        ${OUTPUT_PREFIX}/clairs/${TUMOR}
        
# 12. methylation:

