# ====================================================================
# COPY NUMBER
# ====================================================================

rule copynumber:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        reference=config.get("chm13_fasta", ""),
    output:
        "copynumber/{tumor}/output/{tumor}.copynumber.png"
    message:
        "--- Running copy number analysis for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        work_dir="copynumber/{tumor}/workspace",
        output_dir="copynumber/{tumor}/output",
        hap1_satellite=config.get("hap1_satellite", ""),
        hap2_satellite=config.get("hap2_satellite", ""),
        sex=config.get("sex", "female"),
        binwidth=config.get("copynumber_binwidth", 0.05),
        hap1_label=config.get("copynumber_hap1_label", "Haplotype1"),
        hap2_label=config.get("copynumber_hap2_label", "Haplotype2"),
        ploidy_hap1=config.get("copynumber_ploidy_hap1", ""),
        ploidy_hap2=config.get("copynumber_ploidy_hap2", ""),
        censat_bed=config.get("censat_bed", ""),
        plot_sex_chrom=str(config.get("copynumber_plot_sex_chrom", True)),
        chm13_censat=config.get("chm13_censat", "")
    threads:
        get_threads("copynumber", 8)
    resources:
        mem_mb=get_mem_mb("copynumber", 32000)
    log:
        "logs/copynumber/{tumor}.log"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber.sh \
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
            {threads} \
            {params.binwidth} \
            {params.hap1_label} \
            {params.hap2_label} \
            "{params.ploidy_hap1}" \
            "{params.ploidy_hap2}" \
            "{params.censat_bed}" \
            "{params.plot_sex_chrom}" \
            "{params.chm13_censat}" &> {log}
        """
