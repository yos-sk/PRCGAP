# ====================================================================
# CLAIRS
# ====================================================================

rule clairs:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)),
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
        /bin/bash {SCRIPTS_DIR}/clairs.sh \
            {input.tumor_bam} \
            {input.normal_bam} \
            {output} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {threads} &> {log}
        """

# ====================================================================
# CLAIRS POSTPROCESS PIPELINE
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/realign.sh \
            {input.realign_fasta} \
            {params.output_dir} \
            {threads} &> {log}
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
    log:
        "logs/clairs_postprocess/{tumor}_split_bed.log"
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
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/pileup.sh \
            {input.bam} \
            {input.pileup_tasks} \
            {input.reference_fa} \
            {params.tumor} \
            {params.output_dir} \
            {threads} &> {log}
        """

rule clairs_postprocess_haplotype:
    input:
        parsed_bed="clairs_post/{tumor}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        realign_bam="clairs_post/{tumor}/realign/realign.bam",
        pileup="clairs_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
        kmer_ratio=lambda wc: "bam_refiner/{}/hifi/{}_kmer_ratio.txt".format(wc.tumor, wc.tumor)
    output:
        haplotyped="clairs_post/{tumor}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for ClairS postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="clairs_post/{tumor}",
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/haplotyping.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} \
            {input.kmer_ratio} &> {log}
        """
