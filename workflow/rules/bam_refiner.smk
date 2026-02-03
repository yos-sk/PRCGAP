# ====================================================================
# BAM REFINER
# ====================================================================

rule bam_refiner_hifi:
    input:
        file=lambda wc: samples.loc[wc.sample, "hifi"],
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        bam="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam",
        bai="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam.bai",
        kmer_ratio="bam_refiner/{sample}/hifi/{sample}_kmer_ratio.txt"
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
        bai="bam_refiner/{sample}/ont/{sample}_bam_refined.sorted.bam.bai",
        kmer_ratio="bam_refiner/{sample}/ont/{sample}_kmer_ratio.txt"
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
