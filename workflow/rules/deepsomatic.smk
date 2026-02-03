# ====================================================================
# DEEPSOMATIC
# ====================================================================

rule deepsomatic:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)),
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
# DEEPSOMATIC POSTPROCESS PIPELINE
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/deepsomatic_prepare.sh \
            {input.deepsomatic_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.output_dir} \
            {output.vcf} \
            {output.reference_fa} &> {log}
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/parse_vcf.sh \
            {SCRIPTS_DIR} \
            {input.vcf} \
            {input.reference_fa} \
            {output.parsed_bed} \
            {output.realign_fasta} \
            {params.tool} \
            {params.output_dir} &> {log}
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/realign.sh \
            {input.realign_fasta} \
            {params.output_dir} \
            {threads} &> {log}
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
    log:
        "logs/deepsomatic_postprocess/{tumor}_split_bed.log"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/split_bed.sh \
            {input.parsed_bed} \
            {params.output_dir} \
            {output.pileup_tasks} &> {log}
        """

rule deepsomatic_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/pileup.sh \
            {input.bam} \
            {input.pileup_tasks} \
            {input.reference_fa} \
            {params.tumor} \
            {params.output_dir} \
            {threads} &> {log}
        """

rule deepsomatic_postprocess_haplotype:
    input:
        parsed_bed="deepsomatic_post/{tumor}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        realign_bam="deepsomatic_post/{tumor}/realign/realign.bam",
        pileup="deepsomatic_post/{tumor}/pileup/{tumor}_pileup.bed.gz",
        kmer_ratio=lambda wc: "bam_refiner/{}/hifi/{}_kmer_ratio.txt".format(wc.tumor, wc.tumor),
    output:
        haplotyped="deepsomatic_post/{tumor}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for DeepSomatic postprocess {wildcards.tumor}"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}",
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
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/haplotyping.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} \
            {input.kmer_ratio} &> {log}
        """
