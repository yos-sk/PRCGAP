# ====================================================================
# METHYLATION
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
