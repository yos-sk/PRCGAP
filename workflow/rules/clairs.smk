# ====================================================================
# CLAIRS
# ====================================================================
#
# ClairS runs per available sequencing type ({seqtype} = hifi / ont). The
# per-seqtype tumor/normal BAMs come from bam_refiner. The ClairS platform
# model is selected via config["clairs_model"] (default hifi_sequel2); valid
# values: hifi_sequel2, hifi_revio, ont_r10_dorado_sup_5khz_ssrs,
# ont_r10_dorado_sup_4khz.

rule clairs:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), wc.seqtype, get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        directory("clairs/{tumor}/{seqtype}")
    message:
        "--- Running ClairS for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        model=config.get("clairs_model", "hifi_sequel2"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs", 16)
    resources:
        mem_mb=get_mem_mb("clairs", 64000)
    log:
        "logs/clairs/{tumor}_{seqtype}.log"
    singularity:
        config.get("singularity_images", {}).get("clairs", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/clairs.sh \
            {input.tumor_bam} \
            {input.normal_bam} \
            {output} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {threads} \
            {params.model} &> {log}
        """

# ====================================================================
# CLAIRS POSTPROCESS PIPELINE
# ====================================================================

rule clairs_postprocess_prepare:
    input:
        clairs_dir="clairs/{tumor}/{seqtype}",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        vcf="clairs_post/{tumor}/{seqtype}/output.vcf.gz",
        reference_fa="clairs_post/{tumor}/{seqtype}/reference.fa",
        reference_fai="clairs_post/{tumor}/{seqtype}/reference.fa.fai",
    message:
        "--- Preparing ClairS results for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_prepare.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/clairs_prepare.sh \
            {input.clairs_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.output_dir} \
            {output.vcf} \
            {output.reference_fa} &> {log}
        """

rule clairs_postprocess_parse_vcf:
    input:
        vcf="clairs_post/{tumor}/{seqtype}/output.vcf.gz",
        reference_fa="clairs_post/{tumor}/{seqtype}/reference.fa",
    output:
        parsed_bed="clairs_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        realign_fasta="clairs_post/{tumor}/{seqtype}/realign/realign_ref.fasta",
    message:
        "--- Parsing ClairS VCF for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
        tool="ClairS",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_parse_vcf.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/parse_vcf.sh \
            {SCRIPTS_DIR} \
            {input.vcf} \
            {input.reference_fa} \
            {output.parsed_bed} \
            {output.realign_fasta} \
            {params.tool} \
            {params.output_dir} &> {log}
        """

rule clairs_postprocess_realign:
    input:
        parsed_bed="clairs_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        realign_fasta="clairs_post/{tumor}/{seqtype}/realign/realign_ref.fasta",
    output:
        realign_bam="clairs_post/{tumor}/{seqtype}/realign/realign.bam",
        realign_bai="clairs_post/{tumor}/{seqtype}/realign/realign.bam.bai",
    message:
        "--- Realigning ClairS results for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_realign", 16)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_realign", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_realign.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/realign.sh \
            {input.realign_fasta} \
            {params.output_dir} \
            {threads} &> {log}
        """

rule clairs_postprocess_split_bed:
    input:
        parsed_bed="clairs_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
    output:
        pileup_tasks="clairs_post/{tumor}/{seqtype}/pileup/workspace/pileup_tasks.txt",
    message:
        "--- Splitting ClairS BED for pileup for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_split", 1)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_split", 4000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_split_bed.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/split_bed.sh \
            {input.parsed_bed} \
            {params.output_dir} \
            {output.pileup_tasks} &> {log}
        """

rule clairs_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        pileup_tasks="clairs_post/{tumor}/{seqtype}/pileup/workspace/pileup_tasks.txt",
        reference_fa="clairs_post/{tumor}/{seqtype}/reference.fa",
    output:
        pileup="clairs_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz",
        pileup_tbi="clairs_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz.tbi",
    message:
        "--- Running pileup for ClairS postprocess {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
        no_baq=config.get("pileup_no_baq", "false"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_pileup", 16)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_pileup", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_pileup.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/pileup.sh \
            {input.bam} \
            {input.pileup_tasks} \
            {input.reference_fa} \
            {params.tumor} \
            {params.output_dir} \
            {threads} \
            {resources.mem_mb} \
            {params.no_baq} &> {log}
        """

rule clairs_postprocess_haplotype:
    input:
        parsed_bed="clairs_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        realign_bam="clairs_post/{tumor}/{seqtype}/realign/realign.bam",
        pileup="clairs_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz",
        kmer_ratio=lambda wc: "bam_refiner/{}/{}/{}_kmer_ratio.txt".format(wc.tumor, wc.seqtype, wc.tumor),
    output:
        haplotyped="clairs_post/{tumor}/{seqtype}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for ClairS postprocess {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_haplotype", 4)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_haplotype", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_haplotype.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/haplotyping.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} \
            {input.kmer_ratio} &> {log}
        """
