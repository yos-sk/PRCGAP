import os

SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")

# ====================================================================
# 1. BAM REFINER
# ====================================================================

rule bam_refiner_hifi:
    input:
        file=lambda wc: samples.loc[wc.sample, "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        bam="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam",
        bai="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam.bai"
    message:
        "--- Running bam_refiner for {wildcards.sample} HiFi data"
    params:
        sample="{sample}",
        output_dir="bam_refiner/{sample}/hifi"
    threads:
        get_threads("bam_refiner", 16)
    resources:
        mem_mb=get_mem_mb("bam_refiner", 64000)
    log:
        "logs/bam_refiner/{sample}_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("bam_refiner", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/bam_refiner.sh \
            {params.sample} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {input.file} \
            {params.output_dir} \
            {threads} \
            hifi &> {log}
        """

rule bam_refiner_ont:
    input:
        file=lambda wc: samples.loc[wc.sample, "ont"],
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        bam="bam_refiner/{sample}/ont/{sample}_bam_refined.sorted.bam",
        bai="bam_refiner/{sample}/ont/{sample}_bam_refined.sorted.bam.bai"
    message:
        "--- Running bam_refiner for {wildcards.sample} ONT data"
    params:
        sample="{sample}",
        output_dir="bam_refiner/{sample}/ont"
    threads:
        get_threads("bam_refiner", 16)
    resources:
        mem_mb=get_mem_mb("bam_refiner", 64000)
    log:
        "logs/bam_refiner/{sample}_ont.log"
    singularity:
        config.get("singularity_images", {}).get("bam_refiner", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/bam_refiner.sh \
            {params.sample} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {input.file} \
            {params.output_dir} \
            {threads} \
            ont &> {log}
        """

# ====================================================================
# 2. METHYLATION
# ====================================================================

rule methylation_hifi:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample) if steps["bam_refiner"] else samples.loc[wc.sample, "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        directory("methylation/{sample}/hifi/output")
    message:
        "--- Running methylation for {wildcards.sample} HiFi data"
    params:
        sample="{sample}",
        work_dir="methylation/{sample}/hifi/workspace",
        output_dir="methylation/{sample}/hifi/output",
        type="HiFi"
    threads:
        get_threads("methylation", 8)
    resources:
        mem_mb=get_mem_mb("methylation", 32000)
    log:
        "logs/methylation/{sample}_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("methylation", "")
    shell:
        """
        mkdir -p {params.work_dir} {params.output_dir}
        /bin/bash {SCRIPTS_DIR}/methylation.sh \
            {params.sample} \
            {input.bam} \
            {params.work_dir} \
            {params.output_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.type} \
            {threads} &> {log}
        """

rule methylation_ont:
    input:
        bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample) if steps["bam_refiner"] else samples.loc[wc.sample, "ont"],
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        directory("methylation/{sample}/ont/output")
    message:
        "--- Running methylation for {wildcards.sample} ONT data"
    params:
        sample="{sample}",
        work_dir="methylation/{sample}/ont/workspace",
        output_dir="methylation/{sample}/ont/output",
        type="ONT"
    threads:
        get_threads("methylation", 8)
    resources:
        mem_mb=get_mem_mb("methylation", 32000)
    log:
        "logs/methylation/{sample}_ont.log"
    singularity:
        config.get("singularity_images", {}).get("methylation", "")
    shell:
        """
        mkdir -p {params.work_dir} {params.output_dir}
        /bin/bash {SCRIPTS_DIR}/methylation.sh \
            {params.sample} \
            {input.bam} \
            {params.work_dir} \
            {params.output_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.type} \
            {threads} &> {log}
        """

# ====================================================================
# 3. COPY NUMBER
# ====================================================================

rule copynumber:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)) if steps["bam_refiner"] else samples.loc[get_paired_normal(wc.tumor), "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        reference=config.get("reference", ""),
    output:
        directory("copynumber/{tumor}/output")
    message:
        "--- Running copy number analysis for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        work_dir="copynumber/{tumor}/workspace",
        output_dir="copynumber/{tumor}/output",
        hap1_satellite=config.get("hap1_satellite", ""),
        hap2_satellite=config.get("hap2_satellite", ""),
        sex=config.get("sex", "female")
    threads:
        get_threads("copynumber", 8)
    resources:
        mem_mb=get_mem_mb("copynumber", 32000)
    log:
        "logs/copynumber/{tumor}.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation", "")
    shell:
        """
        mkdir -p {params.work_dir} {params.output_dir}
        SCRIPT_DIR={SCRIPTS_DIR} /bin/bash {SCRIPTS_DIR}/copynumber.sh \
            {params.tumor} \
            {params.normal} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {input.tumor_bam} \
            {input.normal_bam} \
            {input.reference} \
            {params.work_dir} \
            {params.output_dir} \
            {params.hap1_satellite} \
            {params.hap2_satellite} \
            {params.sex} \
            {SCRIPTS_DIR} \
            {threads} &> {log}
        """

# ====================================================================
# 4. NANOMONSV PARSE
# ====================================================================

rule nanomonsv_parse_hifi:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample) if steps["bam_refiner"] else samples.loc[wc.sample, "hifi"],
    output:
        bp_info="nanomonsv/hifi/{sample}.bp_info.sorted.bed.gz",
        bp_info_tbi="nanomonsv/hifi/{sample}.bp_info.sorted.bed.gz.tbi",
        deletion="nanomonsv/hifi/{sample}.deletion.sorted.bed.gz",
        deletion_tbi="nanomonsv/hifi/{sample}.deletion.sorted.bed.gz.tbi",
        insertion="nanomonsv/hifi/{sample}.insertion.sorted.bed.gz",
        insertion_tbi="nanomonsv/hifi/{sample}.insertion.sorted.bed.gz.tbi",
        rearrangement="nanomonsv/hifi/{sample}.rearrangement.sorted.bedpe.gz",
        rearrangement_tbi="nanomonsv/hifi/{sample}.rearrangement.sorted.bedpe.gz.tbi",
    message:
        "--- Running NanoMonSV parse for {wildcards.sample} HiFi data"
    params:
        sample="{sample}",
        output_dir="nanomonsv/hifi"
    threads:
        get_threads("nanomonsv_parse", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_parse", 16000)
    log:
        "logs/nanomonsv/{sample}_parse_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        mkdir -p {params.output_dir}
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_parse.sh \
            {params.sample} \
            {input.bam} \
            {params.output_dir} &> {log}
        """

rule nanomonsv_parse_ont:
    input:
        bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample) if steps["bam_refiner"] else samples.loc[wc.sample, "ont"],
    output:
        bp_info="nanomonsv/ont/{sample}.bp_info.sorted.bed.gz",
        bp_info_tbi="nanomonsv/ont/{sample}.bp_info.sorted.bed.gz.tbi",
        deletion="nanomonsv/ont/{sample}.deletion.sorted.bed.gz",
        deletion_tbi="nanomonsv/ont/{sample}.deletion.sorted.bed.gz.tbi",
        insertion="nanomonsv/ont/{sample}.insertion.sorted.bed.gz",
        insertion_tbi="nanomonsv/ont/{sample}.insertion.sorted.bed.gz.tbi",
        rearrangement="nanomonsv/ont/{sample}.rearrangement.sorted.bedpe.gz",
        rearrangement_tbi="nanomonsv/ont/{sample}.rearrangement.sorted.bedpe.gz.tbi",
    message:
        "--- Running NanoMonSV parse for {wildcards.sample} ONT data"
    params:
        sample="{sample}",
        output_dir="nanomonsv/ont"
    threads:
        get_threads("nanomonsv_parse", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_parse", 16000)
    log:
        "logs/nanomonsv/{sample}_parse_ont.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        mkdir -p {params.output_dir}
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_parse.sh \
            {params.sample} \
            {input.bam} \
            {params.output_dir} &> {log}
        """

# ====================================================================
# 5. NANOMONSV GET
# ====================================================================

rule nanomonsv_get_hifi:
    input:
        tumor_rearrangement=lambda wc: "nanomonsv/hifi/{}.rearrangement.sorted.bedpe.gz".format(wc.tumor) if steps["nanomonsv_parse"] else [],
        normal_rearrangement=lambda wc: "nanomonsv/hifi/{}.rearrangement.sorted.bedpe.gz".format(get_paired_normal(wc.tumor)) if steps["nanomonsv_parse"] else [],
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)) if steps["bam_refiner"] else samples.loc[get_paired_normal(wc.tumor), "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        result="nanomonsv/hifi/{tumor}.nanomonsv.result.txt",
        supporting="nanomonsv/hifi/{tumor}.nanomonsv.supporting_read.txt"
    message:
        "--- Running NanoMonSV get for {wildcards.tumor} HiFi data"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        output_dir="nanomonsv/hifi",
        simple_repeat=config.get("simple_repeat", "")
    threads:
        get_threads("nanomonsv_get", 8)
    resources:
        mem_mb=get_mem_mb("nanomonsv_get", 32000)
    log:
        "logs/nanomonsv/{tumor}_get_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_get.sh \
            {params.tumor} \
            {params.normal} \
            {input.tumor_bam} \
            {input.normal_bam} \
            {params.output_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            hifi \
            {params.simple_repeat} \
            {threads} &> {log}
        """

rule nanomonsv_get_ont:
    input:
        tumor_rearrangement=lambda wc: "nanomonsv/ont/{}.rearrangement.sorted.bedpe.gz".format(wc.tumor) if steps["nanomonsv_parse"] else [],
        normal_rearrangement=lambda wc: "nanomonsv/ont/{}.rearrangement.sorted.bedpe.gz".format(get_paired_normal(wc.tumor)) if steps["nanomonsv_parse"] else [],
        tumor_bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "ont"],
        normal_bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)) if steps["bam_refiner"] else samples.loc[get_paired_normal(wc.tumor), "ont"],
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        result="nanomonsv/ont/{tumor}.nanomonsv.result.txt",
        supporting="nanomonsv/ont/{tumor}.nanomonsv.supporting_read.txt"
    message:
        "--- Running NanoMonSV get for {wildcards.tumor} ONT data"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        output_dir="nanomonsv/ont",
        simple_repeat=config.get("simple_repeat", "")
    threads:
        get_threads("nanomonsv_get", 8)
    resources:
        mem_mb=get_mem_mb("nanomonsv_get", 32000)
    log:
        "logs/nanomonsv/{tumor}_get_ont.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_get.sh \
            {params.tumor} \
            {params.normal} \
            {input.tumor_bam} \
            {input.normal_bam} \
            {params.output_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            ont \
            {params.simple_repeat} \
            {threads} &> {log}
        """

# ====================================================================
# 6. NANOMONSV POSTPROCESS
# ====================================================================

rule nanomonsv_postprocess_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.result.txt",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt"
    message:
        "--- Running NanoMonSV postprocess for {wildcards.tumor} HiFi data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/hifi"
    threads:
        get_threads("nanomonsv_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_postprocess", 16000)
    log:
        "logs/nanomonsv/{tumor}_postprocess_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_postprocess.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} &> {log}
        touch {output}
        """

rule nanomonsv_postprocess_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.result.txt",
        bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "ont"],
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt"
    message:
        "--- Running NanoMonSV postprocess for {wildcards.tumor} ONT data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/ont"
    threads:
        get_threads("nanomonsv_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_postprocess", 16000)
    log:
        "logs/nanomonsv/{tumor}_postprocess_ont.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_postprocess.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} \
            {SCRIPTS_DIR} &> {log}
        touch {output}
        """

# ====================================================================
# 7. NANOMONSV INSERT CLASSIFY
# ====================================================================

rule nanomonsv_insert_classify_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt"
    message:
        "--- Running NanoMonSV insert classification for {wildcards.tumor} HiFi data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/hifi",
        gtf_file=config.get("gtf_file", ""),
        line1_bed=config.get("line1_bed", "")
    threads:
        get_threads("nanomonsv_insert_classify", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_insert_classify", 16000)
    log:
        "logs/nanomonsv/{tumor}_insert_classify_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_insert_classify.sh \
            {input.result} \
            {params.output_dir} \
            {params.tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt \
            {params.gtf_file} \
            {params.line1_bed} &> {log}
        """

rule nanomonsv_insert_classify_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt",
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt"
    message:
        "--- Running NanoMonSV insert classification for {wildcards.tumor} ONT data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/ont",
        gtf_file=config.get("gtf_file", ""),
        line1_bed=config.get("line1_bed", "")
    threads:
        get_threads("nanomonsv_insert_classify", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_insert_classify", 16000)
    log:
        "logs/nanomonsv/{tumor}_insert_classify_ont.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_insert_classify.sh \
            {input.result} \
            {params.output_dir} \
            {params.tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt \
            {params.gtf_file} \
            {params.line1_bed} &> {log}
        """

# ====================================================================
# 8. NANOMONSV CONNECT
# ====================================================================

rule nanomonsv_connect_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
        supporting="nanomonsv/hifi/{tumor}.nanomonsv.supporting_read.txt"
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.connected.txt"
    message:
        "--- Running NanoMonSV connect for {wildcards.tumor} HiFi data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/hifi"
    threads:
        get_threads("nanomonsv_connect", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_connect", 16000)
    log:
        "logs/nanomonsv/{tumor}_connect_hifi.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_connect.sh \
            {input.result} \
            {input.supporting} \
            {params.output_dir} &> {log}
        touch {output}
        """

rule nanomonsv_connect_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt",
        supporting="nanomonsv/ont/{tumor}.nanomonsv.supporting_read.txt"
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.connected.txt"
    message:
        "--- Running NanoMonSV connect for {wildcards.tumor} ONT data"
    params:
        tumor="{tumor}",
        output_dir="nanomonsv/ont"
    threads:
        get_threads("nanomonsv_connect", 4)
    resources:
        mem_mb=get_mem_mb("nanomonsv_connect", 16000)
    log:
        "logs/nanomonsv/{tumor}_connect_ont.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/nanomonsv_connect.sh \
            {input.result} \
            {input.supporting} \
            {params.output_dir} &> {log}
        touch {output}
        """

# ====================================================================
# 9. NANOMONSV MERGE
# ====================================================================

rule nanomonsv_merge:
    input:
        hifi="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
        ont="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt"
    output:
        "nanomonsv/{tumor}.nanomonsv.result.merged.txt"
    message:
        "--- Running NanoMonSV merge for {wildcards.tumor}"
    params:
        tumor="{tumor}"
    threads:
        get_threads("nanomonsv_merge", 2)
    resources:
        mem_mb=get_mem_mb("nanomonsv_merge", 8000)
    log:
        "logs/nanomonsv/{tumor}_merge.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/nanomonsv/merge.sh \
            {input.hifi} \
            {input.ont} \
            {output} \
            {SCRIPTS_DIR} &> {log}
        """

# ====================================================================
# 10. CLAIRS
# ====================================================================

rule clairs:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)) if steps["bam_refiner"] else samples.loc[get_paired_normal(wc.tumor), "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        directory("clairs/{tumor}")
    message:
        "--- Running ClairS for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor)
    threads:
        get_threads("clairs", 16)
    resources:
        mem_mb=get_mem_mb("clairs", 64000)
    log:
        "logs/clairs/{tumor}.log"
    singularity:
        config.get("singularity_images", {}).get("clairs", "")
    shell:
        """
        mkdir -p {output}
        /bin/bash {SCRIPTS_DIR}/clairs.sh \
            {input.tumor_bam} \
            {input.normal_bam} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {output} \
            {threads} &> {log}
        """

# ====================================================================
# 11. DEEPSOMATIC
# ====================================================================

rule deepsomatic:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)) if steps["bam_refiner"] else samples.loc[get_paired_normal(wc.tumor), "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        directory("deepsomatic/{tumor}")
    message:
        "--- Running DeepSomatic for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor)
    threads:
        get_threads("deepsomatic", 16)
    resources:
        mem_mb=get_mem_mb("deepsomatic", 64000)
    log:
        "logs/deepsomatic/{tumor}.log"
    singularity:
        config.get("singularity_images", {}).get("deepsomatic", "")
    shell:
        """
        mkdir -p {output}
        /bin/bash {SCRIPTS_DIR}/deepsomatic.sh \
            {params.tumor} \
            {params.normal} \
            {input.tumor_bam} \
            {input.normal_bam} \
            {output} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {threads} &> {log}
        """

# ====================================================================
# 12. CLAIRS POSTPROCESS PIPELINE
# ====================================================================

rule clairs_postprocess_prepare:
    input:
        clairs_dir="clairs/{tumor}",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        vcf="clairs_post/{tumor}/output.vcf.gz",
        reference_fa="clairs_post/{tumor}/reference.fa",
        reference_fai="clairs_post/{tumor}/reference.fa.fai",
    message:
        "--- Preparing ClairS results for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
        output_dir="clairs_post/{tumor}",
    threads:
        get_threads("clairs_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_prepare.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.assembly_hap1} {input.assembly_hap2} > {output.reference_fa}
        samtools faidx {output.reference_fa}

        # Combine ClairS VCF (if indel.vcf.gz exists)
        if [ -f {input.clairs_dir}/indel.vcf.gz ]; then
            zcat {input.clairs_dir}/output.vcf.gz | gzip -c > {output.vcf}
            zgrep -v "#" {input.clairs_dir}/indel.vcf.gz | gzip -f -c >> {output.vcf}
        else
            cp {input.clairs_dir}/output.vcf.gz {output.vcf}
        fi
        &> {log}
        """

rule clairs_postprocess_parse_vcf:
    input:
        vcf="clairs_post/{tumor}/output.vcf.gz",
        reference_fa="clairs_post/{tumor}/reference.fa",
    output:
        parsed_bed="clairs_post/{tumor}/realign/parsed_vcf.bed",
        realign_fasta="clairs_post/{tumor}/realign/realign_ref.fasta",
    message:
        "--- Parsing ClairS VCF for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
        tool="ClairS",
    threads:
        get_threads("clairs_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_parse_vcf.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}/realign
        python3 {SCRIPTS_DIR}/mutation_postprocess/parse_vcf.py \
            -i {input.vcf} \
            -f {input.reference_fa} \
            -o {output.parsed_bed} \
            -g {output.realign_fasta} \
            -t {params.tool} &> {log}
        """

rule clairs_postprocess_realign:
    input:
        parsed_bed="clairs_post/{tumor}/realign/parsed_vcf.bed",
        realign_fasta="clairs_post/{tumor}/realign/realign_ref.fasta",
    output:
        realign_bam="clairs_post/{tumor}/realign/realign.bam",
        realign_bai="clairs_post/{tumor}/realign/realign.bam.bai",
    message:
        "--- Realigning ClairS results for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("clairs_postprocess_realign", 16)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_realign", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_realign.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        bwa index {input.realign_fasta}
        bwa mem -a -k 50 -c 1000000 -t {threads} {input.realign_fasta} {input.realign_fasta} | \
            samtools view -Shb > {params.output_dir}/realign/realign.unsorted
        samtools sort -@ {threads} {params.output_dir}/realign/realign.unsorted > {output.realign_bam}
        samtools index -@ {threads} {output.realign_bam}
        rm {params.output_dir}/realign/realign.unsorted &> {log}
        """

rule clairs_postprocess_split_bed:
    input:
        parsed_bed="clairs_post/{tumor}/realign/parsed_vcf.bed",
    output:
        pileup_tasks="clairs_post/{tumor}/pileup/workspace/pileup_tasks.txt",
    message:
        "--- Splitting ClairS BED for pileup for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
    threads:
        get_threads("clairs_postprocess_split", 1)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_split", 4000)
    shell:
        """
        mkdir -p {params.output_dir}/pileup/workspace
        split -n l/16 {input.parsed_bed} {params.output_dir}/pileup/workspace/input_
        ls {params.output_dir}/pileup/workspace/input_* | while read f; do
            echo -e "${{f}}\t${{f}}.pileup"
        done > {output.pileup_tasks}
        """

rule clairs_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        pileup_tasks="clairs_post/{tumor}/pileup/workspace/pileup_tasks.txt",
        reference_fa="clairs_post/{tumor}/reference.fa",
    output:
        pileup="clairs_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
        pileup_tbi="clairs_post/{tumor}/pileup/{tumor}_pileup.bed.gz.tbi",
    message:
        "--- Running pileup for ClairS postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("clairs_postprocess_pileup", 16)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_pileup", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_pileup.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}/pileup/workspace
        while IFS=$'\\t' read bed_file pileup_file; do
            samtools mpileup -l ${{bed_file}} -f {input.reference_fa} --output-QNAME {input.bam} | \
                awk '{{print $1 "\\t" $2 -1 "\\t" $2 "\\t" $4 "\\t" $5 "\\t" $6 "\\t" $7}}' > ${{pileup_file}}
        done < {input.pileup_tasks}

        cat {params.output_dir}/pileup/workspace/*.pileup | sort -k 1,1 -k 2,2n > {params.output_dir}/pileup/{params.tumor}_pileup.bed
        bgzip -@ {threads} -f {params.output_dir}/pileup/{params.tumor}_pileup.bed
        tabix -p bed {output.pileup} &> {log}
        """

rule clairs_postprocess_haplotype:
    input:
        parsed_bed="clairs_post/{tumor}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        realign_bam="clairs_post/{tumor}/realign/realign.bam",
        pileup="clairs_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
    output:
        haplotyped="clairs_post/{tumor}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for ClairS postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("clairs_postprocess_haplotype", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_haplotype", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_haplotype.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mutation_postprocess haplotype \
            -i {input.parsed_bed} \
            -b {input.bam} \
            -r {input.realign_bam} \
            -p {input.pileup} \
            -d 0.98 \
            1>{output.haplotyped} 2>>{log}
        """

# ====================================================================
# 13. DEEPSOMATIC POSTPROCESS PIPELINE
# ====================================================================

rule deepsomatic_postprocess_prepare:
    input:
        deepsomatic_dir="deepsomatic/{tumor}",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        vcf="deepsomatic_post/{tumor}/output.vcf.gz",
        reference_fa="deepsomatic_post/{tumor}/reference.fa",
        reference_fai="deepsomatic_post/{tumor}/reference.fa.fai",
    message:
        "--- Preparing DeepSomatic results for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
        output_dir="deepsomatic_post/{tumor}",
    threads:
        get_threads("deepsomatic_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_prepare.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.assembly_hap1} {input.assembly_hap2} > {output.reference_fa}
        samtools faidx {output.reference_fa}
        cp {input.deepsomatic_dir}/output.vcf.gz {output.vcf}
        &> {log}
        """

rule deepsomatic_postprocess_parse_vcf:
    input:
        vcf="deepsomatic_post/{tumor}/output.vcf.gz",
        reference_fa="deepsomatic_post/{tumor}/reference.fa",
    output:
        parsed_bed="deepsomatic_post/{tumor}/realign/parsed_vcf.bed",
        realign_fasta="deepsomatic_post/{tumor}/realign/realign_ref.fasta",
    message:
        "--- Parsing DeepSomatic VCF for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
        tool="DeepSomatic",
    threads:
        get_threads("deepsomatic_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_parse_vcf.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}/realign
        python3 {SCRIPTS_DIR}/mutation_postprocess/parse_vcf.py \
            -i {input.vcf} \
            -f {input.reference_fa} \
            -o {output.parsed_bed} \
            -g {output.realign_fasta} \
            -t {params.tool} &> {log}
        """

rule deepsomatic_postprocess_realign:
    input:
        parsed_bed="deepsomatic_post/{tumor}/realign/parsed_vcf.bed",
        realign_fasta="deepsomatic_post/{tumor}/realign/realign_ref.fasta",
    output:
        realign_bam="deepsomatic_post/{tumor}/realign/realign.bam",
        realign_bai="deepsomatic_post/{tumor}/realign/realign.bam.bai",
    message:
        "--- Realigning DeepSomatic results for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("deepsomatic_postprocess_realign", 16)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_realign", 32000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_realign.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        bwa index {input.realign_fasta}
        bwa mem -a -k 50 -c 1000000 -t {threads} {input.realign_fasta} {input.realign_fasta} | \
            samtools view -Shb > {params.output_dir}/realign/realign.unsorted
        samtools sort -@ {threads} {params.output_dir}/realign/realign.unsorted > {output.realign_bam}
        samtools index -@ {threads} {output.realign_bam}
        rm {params.output_dir}/realign/realign.unsorted &> {log}
        """

rule deepsomatic_postprocess_split_bed:
    input:
        parsed_bed="deepsomatic_post/{tumor}/realign/parsed_vcf.bed",
    output:
        pileup_tasks="deepsomatic_post/{tumor}/pileup/workspace/pileup_tasks.txt",
    message:
        "--- Splitting DeepSomatic BED for pileup for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
    threads:
        get_threads("deepsomatic_postprocess_split", 1)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_split", 4000)
    shell:
        """
        mkdir -p {params.output_dir}/pileup/workspace
        split -n l/16 {input.parsed_bed} {params.output_dir}/pileup/workspace/input_
        ls {params.output_dir}/pileup/workspace/input_* | while read f; do
            echo -e "${{f}}\t${{f}}.pileup"
        done > {output.pileup_tasks}
        """

rule deepsomatic_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        pileup_tasks="deepsomatic_post/{tumor}/pileup/workspace/pileup_tasks.txt",
        reference_fa="deepsomatic_post/{tumor}/reference.fa",
    output:
        pileup="deepsomatic_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
        pileup_tbi="deepsomatic_post/{tumor}/pileup/{tumor}_pileup.bed.gz.tbi",
    message:
        "--- Running pileup for DeepSomatic postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("deepsomatic_postprocess_pileup", 16)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_pileup", 32000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_pileup.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mkdir -p {params.output_dir}/pileup/workspace
        while IFS=$'\\t' read bed_file pileup_file; do
            samtools mpileup -l ${{bed_file}} -f {input.reference_fa} --output-QNAME {input.bam} | \
                awk '{{print $1 "\\t" $2 -1 "\\t" $2 "\\t" $4 "\\t" $5 "\\t" $6 "\\t" $7}}' > ${{pileup_file}}
        done < {input.pileup_tasks}

        cat {params.output_dir}/pileup/workspace/*.pileup | sort -k 1,1 -k 2,2n > {params.output_dir}/pileup/{params.tumor}_pileup.bed
        bgzip -@ {threads} -f {params.output_dir}/pileup/{params.tumor}_pileup.bed
        tabix -p bed {output.pileup} &> {log}
        """

rule deepsomatic_postprocess_haplotype:
    input:
        parsed_bed="deepsomatic_post/{tumor}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor) if steps["bam_refiner"] else samples.loc[wc.tumor, "hifi"],
        realign_bam="deepsomatic_post/{tumor}/realign/realign.bam",
        pileup="deepsomatic_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
    output:
        haplotyped="deepsomatic_post/{tumor}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for DeepSomatic postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
        mutation_postprocess_image=config.get("singularity_images", {}).get("point_mutation_postprocess", ""),
    threads:
        get_threads("deepsomatic_postprocess_haplotype", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_haplotype", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_haplotype.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        mutation_postprocess haplotype \
            -i {input.parsed_bed} \
            -b {input.bam} \
            -r {input.realign_bam} \
            -p {input.pileup} \
            -d 0.98 \
            1>{output.haplotyped} 2>>{log}
        """
