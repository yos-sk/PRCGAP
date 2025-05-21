#!/bin/bash

set -ex

CONF=$1
source ${CONF}

DATE=`date +"%Y%m%d"`

# 1. bam_refiner
if [ ${BAM_REFINER} = "TRUE" ]; then
    # 1-1. bam_refiner for hifi
    for SAMPLE in ${TUMOR} ${CONTROL}; do
        if [ ${SAMPLE} = ${TUMOR} ]; then
            INPUT=${TUMOR_HIFI}
        else
            INPUT=${CONTROL_HIFI}
        fi
        sbatch -J ${SAMPLE}_${ASSEMBLER}_bam_refiner_hifi -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_bam_refiner_hifi.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_bam_refiner_hifi.out \
            scripts/bam_refiner/run_bam_refiner.sh \
            ${SAMPLE} \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            ${INPUT} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output \
            hifi
    done

    # 1-2. bam_refiner_hifi for methylation (optional) 
    if [ ${TUMOR_METHYLATION_HIFI} != ${TUMOR_HIFI} ]; then
        sbatch -J ${TUMOR}_${ASSEMBLER}_bam_refiner_methylation_hifi -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_bam_refiner_methylation_hifi.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_bam_refiner_methylation_hifi.out \
            scripts/bam_refiner/run_bam_refiner.sh \
            ${TUMOR} \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            ${TUMOR_METHYLATION_HIFI}\
            ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/methylation_hifi/output \
            hifi
    fi    

    if [ ${CONTROL_METHYLATION_HIFI} != ${CONTROL_HIFI} ]; then
        sbatch -J ${CONTROL}_${ASSEMBLER}_bam_refiner_methylation_hifi -e log/${DATE}/${OUTPUT}/${CONTROL}_${ASSEMBLER}_bam_refiner_methylation_hifi.err -o log/${DATE}/${OUTPUT}/${CONTROL}_${ASSEMBLER}_bam_refiner_methylation_hifi.out \
            scripts/bam_refiner/run_bam_refiner.sh \
            ${CONTROL} \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            ${CONTROL_METHYLATION_HIFI} \
            ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/methylation_hifi/output \
            hifi
    fi    

    # 1-3. bam_refiner_ont
    for SAMPLE in ${TUMOR} ${CONTROL}; do
        if [ ${SAMPLE} = ${TUMOR} ]; then
            INPUT=${TUMOR_ONT}
        else
            INPUT=${CONTROL_ONT}
        fi
        sbatch -J ${SAMPLE}_${ASSEMBLER}_bam_refiner_ont -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_bam_refiner_ont.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_bam_refiner_ont.out \
            scripts/bam_refiner/run_bam_refiner.sh \
            ${SAMPLE} \
            ${ASSEMBLY_HAP1} \
            ${ASSEMBLY_HAP2} \
            ${INPUT} \
            ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output \
            ont
    done
fi

# 2. methylation:
if [ ${METHYLATION} = "TRUE" ]; then
    for SAMPLE in ${TUMOR} ${CONTROL}; do
        if [ ${BAM_REFINER} = "TRUE" ]; then
            if [ ${TUMOR_HIFI} != ${TUMOR_METHYLATION_HIFI} ]; then
                BAM_REFINER_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_${ASSEMBLER}_bam_refiner_methylation_hifi) | cut -d ' ' -f 2)
                INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/methylation_hifi/output/${SAMPLE}_bam_refined.sorted.bam
            else
                BAM_REFINER_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
                INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam
            fi
            
            sbatch --dependency=afterok:${BAM_REFINER_HIFI_JOBID} -J ${SAMPLE}_${ASSEMBLER}_methylation_hifi -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_hifi.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_hifi.out \
                scripts/methylation/run_methylation.sh \
                    ${SAMPLE} \
                    ${INPUT_BAM} \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/workspace \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/output \
                    ${ASSEMBLY_HAP1} \
                    ${ASSEMBLY_HAP2} \
                    HiFi

            BAM_REFINER_ONT_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_${ASSEMBLER}_bam_refiner_ont) | cut -d ' ' -f 2)
            sbatch --dependency=afterok:${BAM_REFINER_ONT_JOBID} -J ${SAMPLE}_${ASSEMBLER}_methylation_ont -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_ont.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_ont.out \
                scripts/methylation/run_methylation.sh \
                    ${SAMPLE} \
                    ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/workspace \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/output \
                    ${ASSEMBLY_HAP1} \
                    ${ASSEMBLY_HAP2} \
                    ONT
        else
            if [ ${TUMOR_HIFI} != ${TUMOR_METHYLATION_HIFI} ]; then
                INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/methylation_hifi/output/${SAMPLE}_bam_refined.sorted.bam
            else
                INPUT_BAM=${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam
            fi
            
            sbatch -J ${SAMPLE}_${ASSEMBLER}_methylation_hifi -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_hifi.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_hifi.out \
                scripts/methylation/run_methylation.sh \
                    ${SAMPLE} \
                    ${INPUT_BAM} \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/workspace \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/hifi/output \
                    ${ASSEMBLY_HAP1} \
                    ${ASSEMBLY_HAP2} \
                    HiFi

            sbatch -J ${SAMPLE}_${ASSEMBLER}_methylation_ont -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_ont.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_methylation_ont.out \
                scripts/methylation/run_methylation.sh \
                    ${SAMPLE} \
                    ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/workspace \
                    ${OUTPUT_PREFIX}/methylation/${SAMPLE}/ont/output \
                    ${ASSEMBLY_HAP1} \
                    ${ASSEMBLY_HAP2} \
                    ONT
        fi

    done
fi

# 3. copynumber
if [ ${COPYNUMBER} = "TRUE" ]; then
    if [ ${BAM_REFINER} = "TRUE" ]; then
        TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2) 
        NORMAL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${NORMAL_BAM_REFINER_JOBID} -J ${TUMOR}_${ASSEMBLER}_copynumber -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_copynumber.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_copynumber.out \
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
                ${HAP1_SATELLITE} \
                ${HAP2_SATELLITE} \
                ${SEX}
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_copynumber -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_copynumber.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_copynumber.out \
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
                ${HAP1_SATELLITE} \
                ${HAP2_SATELLITE} \
                ${SEX}
    fi
fi

# 4. nanomonsv parse 
if [ ${NANOMONSV_PARSE} = "TRUE" ]; then
    for SAMPLE in ${TUMOR} ${CONTROL}; do
        if [ ${BAM_REFINER} = "TRUE" ]; then
            BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
            sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi.out \
                scripts/nanomonsv/run_nanomonsv_parse.sh \
                    ${SAMPLE} \
                    ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam \
                    ${OUTPUT_PREFIX}/nanomonsv/hifi

            BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${SAMPLE}_${ASSEMBLER}_bam_refiner_ont) | cut -d ' ' -f 2)
            sbatch --dependency=afterok:${BAM_REFINER_JOBID} -J ${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont.out \
                scripts/nanomonsv/run_nanomonsv_parse.sh \
                ${SAMPLE} \
                ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/nanomonsv/ont
        else
            sbatch -J ${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_hifi.out \
                scripts/nanomonsv/run_nanomonsv_parse.sh \
                    ${SAMPLE} \
                    ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/hifi/output/${SAMPLE}_bam_refined.sorted.bam \
                    ${OUTPUT_PREFIX}/nanomonsv/hifi
                    
            sbatch -J ${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont -e log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont.err -o log/${DATE}/${OUTPUT}/${SAMPLE}_${ASSEMBLER}_nanomonsv_parse_ont.out \
                scripts/nanomonsv/run_nanomonsv_parse.sh \
                ${SAMPLE} \
                ${OUTPUT_PREFIX}/bam_refiner/${SAMPLE}/ont/output/${SAMPLE}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/nanomonsv/ont
        fi
    done
fi
            
# 5. nanomonsv get 
if [ ${NANOMONSV_GET} = "TRUE" ]; then
    if [ ${NANOMONSV_PARSE} = "TRUE" ]; then
        TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
        CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_${ASSEMBLER}_nanomonsv_parse_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi.out \
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

        TUMOR_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
        CONTROL_NANOMONSV_PARSE_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_${ASSEMBLER}_nanomonsv_parse_ont) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${TUMOR_NANOMONSV_PARSE_JOBID},${CONTROL_NANOMONSV_PARSE_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont.out \
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
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi.out \
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

        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont.out \
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
    fi
fi

# 6. nanomonsv postprocess
if [ ${NANOMONSV_POST} = "TRUE" ]; then
    if [ ${NANOMONSV_GET} = "TRUE" ]; then
        NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi.out \
            scripts/nanomonsv/run_nanomonsv_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/nanomonsv/hifi \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam

        NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont.out \
            scripts/nanomonsv/run_nanomonsv_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/nanomonsv/ont \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.sorted.bam
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi.out \
            scripts/nanomonsv/run_nanomonsv_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/nanomonsv/hifi \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam
                
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont.out \
            scripts/nanomonsv/run_nanomonsv_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/nanomonsv/ont \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/ont/output/${TUMOR}_bam_refined.sorted.bam
    fi
fi

# 7. nanomonsv insert_classify
if [ ${NANOMONSV_INSERT_CLASSIFY} = "TRUE" ]; then
    if [ ${NANOMONSV_GET} = "TRUE" ]; then
        NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_get_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi.log -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi.log \
            scripts/nanomonsv/run_nanomonsv_insert_classify.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi \
                ${TUMOR}.nanomonsv.new_result.sv_typed.insert_classified.txt \
                ${GTF_FILE} \
                ${LINE1_BED} 

        NANOMONSV_GET_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_get_ont) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_GET_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont.log -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont.log \
            scripts/nanomonsv/run_nanomonsv_insert_classify.sh \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont \
                ${TUMOR}.nanomonsv.new_result.sv_typed.insert_classified.txt \
                ${GTF_FILE} \
                ${LINE1_BED} 
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi.log -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_hifi.log \
            scripts/nanomonsv/run_nanomonsv_insert_classify.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi \
                ${TUMOR}.nanomonsv.new_result.sv_typed.insert_classified.txt \
                ${GTF_FILE} \
                ${LINE1_BED} 

        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont.log -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_insert_classify_ont.log \
            scripts/nanomonsv/run_nanomonsv_insert_classify.sh \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont \
                ${TUMOR}.nanomonsv.new_result.sv_typed.insert_classified.txt \
                ${GTF_FILE} \
                ${LINE1_BED} 
    fi
fi


# 8. nanomonsv connect
if [ ${NANOMONSV_CONNECT} = "TRUE" ]; then
    if [ ${NANOMONSV_POST} = "TRUE" ]; then
        NANOMONSV_POSTPROCESS_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi.out \
            scripts/nanomonsv/run_nanomonsv_connect.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.supporting_read.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR} 

        NANOMONSV_POSTPROCESS_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont.out \
            scripts/nanomonsv/run_nanomonsv_connect.sh \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.supporting_read.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_hifi.out \
            scripts/nanomonsv/run_nanomonsv_connect.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.supporting_read.txt \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR} 

        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_connect_ont.out \
            scripts/nanomonsv/run_nanomonsv_connect.sh \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.supporting_read.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}
    fi
fi

# 9. merge nanomonsv results
if [ ${NANOMONSV_MERGE} = "TRUE" ]; then
    if [ ${NANOMONSV_POST} = "TRUE" ]; then
        NANOMONSV_POSTPROCESS_HIFI_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_hifi) | cut -d ' ' -f 2)
        NANOMONSV_POSTPROCESS_ONT_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_nanomonsv_postprocess_ont) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${NANOMONSV_POSTPROCESS_HIFI_JOBID},${NANOMONSV_POSTPROCESS_ONT_JOBID} -J ${TUMOR}_${ASSEMBLER}_nanomonsv_merge -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_merge.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_merge.out \
            scripts/nanomonsv/run_merge.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/${TUMOR}.nanomonsv.result.merged.txt 
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_nanomonsv_merge -e log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_merge.err -o log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_nanomonsv_merge.out \
            scripts/nanomonsv/run_merge.sh \
                ${OUTPUT_PREFIX}/nanomonsv/hifi/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/ont/${TUMOR}.nanomonsv.new_result.sv_typed.txt \
                ${OUTPUT_PREFIX}/nanomonsv/${TUMOR}.nanomonsv.result.merged.txt 
    fi
fi

# 10. clairs
if [ ${CLAIRS} = "TRUE" ]; then
    if [ ${BAM_REFINER} = "TRUE" ]; then
        TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
        CONTROL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${CONTROL_BAM_REFINER_JOBID} -J ${TUMOR}_${ASSEMBLER}_clairs -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs.out \
            scripts/clairs/run_clairs.sh \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/clairs/${TUMOR} 
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_clairs -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs.out \
            scripts/clairs/run_clairs.sh \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/clairs/${TUMOR}
    fi
fi

# 11. deeosomatic
if [ ${DEEPSOMATIC} = "TRUE" ]; then
    if [ ${BAM_REFINER} = "TRUE" ]; then
        TUMOR_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
        CONTROL_BAM_REFINER_JOBID=$(echo $(squeue -noheader --format %i --name ${CONTROL}_${ASSEMBLER}_bam_refiner_hifi) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${TUMOR_BAM_REFINER_JOBID},${CONTROL_BAM_REFINER_JOBID} -J ${TUMOR}_${ASSEMBLER}_deepsomatic -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic.out \
            scripts/deepsomatic/run_deepsomatic.sh \
                ${TUMOR} \
                ${CONTROL} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/deepsomatic/${TUMOR} \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} 
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_deepsomatic -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic.out \
            scripts/deepsomatic/run_deepsomatic.sh \
                ${TUMOR} \
                ${CONTROL} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/bam_refiner/${CONTROL}/hifi/output/${CONTROL}_bam_refined.sorted.bam \
                ${OUTPUT_PREFIX}/deepsomatic/${TUMOR} \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} 
    fi
fi

# 12. clairs postprocess
if [ ${CLAIRS_POST} = "TRUE" ]; then
    if [ ${CLAIRS} = "TRUE" ]; then
        CLAIRS_JOB_ID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_clairs) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${CLAIRS_JOB_ID} -J ${TUMOR}_${ASSEMBLER}_clairs_postprocess -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs_postprocess.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs_postprocess.out \
            scripts/mutation_postprocess/run_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/clairs/${TUMOR} \ ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/clairs_post/${TUMOR} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ClairS \
                ./images/mutation_postprocess-v0.1.2b.sif
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_clairs_postprocess -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs_postprocess.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_clairs_postprocess.out \
            scripts/mutation_postprocess/run_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/clairs/${TUMOR} \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/clairs_post/${TUMOR} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                ClairS \
                ./images/mutation_postprocess-v0.1.2b.sif
    fi
fi

# 13. deepsomatic postprocess
if [ ${DEEPSOMATIC_POST} = "TRUE" ]; then
    if [ ${CLAIRS} = "TRUE" ]; then
        DEEPSOMATIC_JOB_ID=$(echo $(squeue -noheader --format %i --name ${TUMOR}_${ASSEMBLER}_deepsomatic) | cut -d ' ' -f 2)
        sbatch --dependency=afterok:${DEEPSOMATIC_JOB_ID} -J ${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess.out \
            scripts/mutation_postprocess/run_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/deepsomatic/${TUMOR} \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/deepsomatic_post/${TUMOR} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                DeepSomatic \
                ./images/mutation_postprocess-v0.1.2b.sif
    else
        sbatch -J ${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess -e ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess.err -o ./log/${DATE}/${OUTPUT}/${TUMOR}_${ASSEMBLER}_deepsomatic_postprocess.out \
            scripts/mutation_postprocess/run_postprocess.sh \
                ${TUMOR} \
                ${OUTPUT_PREFIX}/deepsomatic/${TUMOR} \
                ${ASSEMBLY_HAP1} \
                ${ASSEMBLY_HAP2} \
                ${OUTPUT_PREFIX}/deepsomatic_post/${TUMOR} \
                ${OUTPUT_PREFIX}/bam_refiner/${TUMOR}/hifi/output/${TUMOR}_bam_refined.sorted.bam \
                DeepSomatic \
                ./images/mutation_postprocess-v0.1.2b.sif
    fi
fi
