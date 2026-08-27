# ====================================================================
# DEEPSOMATIC (contig scatter -> gather)
# ====================================================================
#
# DeepSomatic runs per available sequencing type ({seqtype} = hifi / ont),
# scattered over contig chunks: every contig >= config["caller_solo_contig_min_bp"]
# gets its own job and the remaining small contigs share one. Per-chunk threads
# and memory come from config["resources"].
#
# --model_type is derived from the seqtype: PACBIO for hifi, ONT for ont.

rule deepsomatic_scatter_setup:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        tumor_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(wc.tumor, wc.seqtype, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), wc.seqtype, get_paired_normal(wc.tumor)),
        normal_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(get_paired_normal(wc.tumor), wc.seqtype, get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        reference_fa="deepsomatic_scatter/{tumor}/{seqtype}/reference.fa",
        reference_fai="deepsomatic_scatter/{tumor}/{seqtype}/reference.fa.fai",
        tumor_bam="deepsomatic_scatter/{tumor}/{seqtype}/tumor.bam",
        tumor_bai="deepsomatic_scatter/{tumor}/{seqtype}/tumor.bam.bai",
        normal_bam="deepsomatic_scatter/{tumor}/{seqtype}/normal.bam",
        normal_bai="deepsomatic_scatter/{tumor}/{seqtype}/normal.bam.bai",
        sanitized="deepsomatic_scatter/{tumor}/{seqtype}/sanitized.flag",
    message:
        "--- Preparing DeepSomatic scatter inputs for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        output_dir="deepsomatic_scatter/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_scatter_setup", 8)
    resources:
        mem_mb=get_mem_mb("deepsomatic_scatter_setup", 16000)
    log:
        "logs/deepsomatic/{tumor}_{seqtype}_scatter_setup.log"
    benchmark:
        "benchmarks/deepsomatic/{tumor}_{seqtype}_scatter_setup.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/caller/scatter_setup.sh \
            {input.tumor_bam} \
            {input.normal_bam} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.output_dir} \
            {threads} &> {log}
        """

# Checkpoint: the chunk count depends on the assembly, so the gather rule below
# discovers it after this has run.
checkpoint deepsomatic_chunks:
    input:
        reference_fai="deepsomatic_scatter/{tumor}/{seqtype}/reference.fa.fai",
    output:
        chunk_dir=directory("deepsomatic_scatter/{tumor}/{seqtype}/chunks"),
    message:
        "--- Building DeepSomatic contig chunks for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        min_bp=config.get("caller_solo_contig_min_bp", 1000000),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_chunks", 1)
    resources:
        mem_mb=get_mem_mb("deepsomatic_chunks", 4000)
    log:
        "logs/deepsomatic/{tumor}_{seqtype}_chunks.log"
    benchmark:
        "benchmarks/deepsomatic/{tumor}_{seqtype}_chunks.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/caller/make_chunks.sh \
            {input.reference_fai} \
            {output.chunk_dir} \
            {params.min_bp} &> {log}
        """

rule deepsomatic_chunk:
    input:
        chunk_dir=lambda wc: checkpoints.deepsomatic_chunks.get(
            tumor=wc.tumor, seqtype=wc.seqtype).output.chunk_dir,
        reference_fa="deepsomatic_scatter/{tumor}/{seqtype}/reference.fa",
        tumor_bam="deepsomatic_scatter/{tumor}/{seqtype}/tumor.bam",
        normal_bam="deepsomatic_scatter/{tumor}/{seqtype}/normal.bam",
    output:
        vcf="deepsomatic_scatter/{tumor}/{seqtype}/out/{chunk}/output.vcf.gz",
    message:
        "--- Running DeepSomatic chunk {wildcards.chunk} for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        output_dir="deepsomatic_scatter/{tumor}/{seqtype}/out/{chunk}",
        regions_bed="deepsomatic_scatter/{tumor}/{seqtype}/chunks/{chunk}.bed",
        model_type=lambda wc: "ONT" if wc.seqtype == "ont" else "PACBIO",
        postprocess_cpus=config.get("deepsomatic_postprocess_variants_cpus", 1),
    wildcard_constraints:
        seqtype="hifi|ont",
        chunk=r"\d+",
    threads:
        get_threads("deepsomatic", 8)
    resources:
        mem_mb=get_mem_mb("deepsomatic", 32000)
    log:
        "logs/deepsomatic/{tumor}_{seqtype}_chunk_{chunk}.log"
    benchmark:
        "benchmarks/deepsomatic/{tumor}_{seqtype}_chunk_{chunk}.tsv"
    singularity:
        config.get("singularity_images", {}).get("deepsomatic", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/deepsomatic/chunk.sh \
            {params.tumor} \
            {params.normal} \
            {input.tumor_bam} \
            {input.normal_bam} \
            {input.reference_fa} \
            {params.regions_bed} \
            {params.output_dir} \
            {threads} \
            {params.model_type} \
            {params.postprocess_cpus} &> {log}
        """

def _deepsomatic_chunk_files(wc, basename):
    chunk_dir = checkpoints.deepsomatic_chunks.get(
        tumor=wc.tumor, seqtype=wc.seqtype).output.chunk_dir
    chunks = sorted(glob_wildcards(os.path.join(chunk_dir, "{chunk}.bed")).chunk)
    return expand(
        "deepsomatic_scatter/{tumor}/{seqtype}/out/{chunk}/" + basename,
        tumor=wc.tumor, seqtype=wc.seqtype, chunk=chunks)

rule deepsomatic_merge:
    input:
        vcfs=lambda wc: _deepsomatic_chunk_files(wc, "output.vcf.gz"),
        reference_fai="deepsomatic_scatter/{tumor}/{seqtype}/reference.fa.fai",
        sanitized="deepsomatic_scatter/{tumor}/{seqtype}/sanitized.flag",
    output:
        vcf="deepsomatic/{tumor}/{seqtype}/output.vcf.gz",
        vcf_tbi="deepsomatic/{tumor}/{seqtype}/output.vcf.gz.tbi",
    message:
        "--- Merging DeepSomatic chunks for {wildcards.tumor} ({wildcards.seqtype})"
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_merge", 1)
    resources:
        mem_mb=get_mem_mb("deepsomatic_merge", 16000),
        disk_mb=8000,
    log:
        "logs/deepsomatic/{tumor}_{seqtype}_merge.log"
    benchmark:
        "benchmarks/deepsomatic/{tumor}_{seqtype}_merge.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/caller/merge_vcf.sh \
            {input.reference_fai} \
            {input.sanitized} \
            {output.vcf} \
            $(( {resources.mem_mb} / 2 )) \
            {input.vcfs} &> {log}
        """

# ====================================================================
# DEEPSOMATIC POSTPROCESS PIPELINE
# ====================================================================

rule deepsomatic_postprocess_prepare:
    input:
        deepsomatic_vcf="deepsomatic/{tumor}/{seqtype}/output.vcf.gz",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        vcf="deepsomatic_post/{tumor}/{seqtype}/output.vcf.gz",
        reference_fa="deepsomatic_post/{tumor}/{seqtype}/reference.fa",
        reference_fai="deepsomatic_post/{tumor}/{seqtype}/reference.fa.fai",
    message:
        "--- Preparing DeepSomatic results for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
        deepsomatic_dir="deepsomatic/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_prepare.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_prepare.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/deepsomatic_prepare.sh \
            {params.deepsomatic_dir} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.output_dir} \
            {output.vcf} \
            {output.reference_fa} &> {log}
        """

rule deepsomatic_postprocess_parse_vcf:
    input:
        vcf="deepsomatic_post/{tumor}/{seqtype}/output.vcf.gz",
        reference_fa="deepsomatic_post/{tumor}/{seqtype}/reference.fa",
    output:
        parsed_bed="deepsomatic_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        realign_fasta="deepsomatic_post/{tumor}/{seqtype}/realign/realign_ref.fasta",
    message:
        "--- Parsing DeepSomatic VCF for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
        tool="DeepSomatic",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_parse_vcf.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_parse_vcf.tsv"
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
        parsed_bed="deepsomatic_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        realign_fasta="deepsomatic_post/{tumor}/{seqtype}/realign/realign_ref.fasta",
    output:
        realign_bam="deepsomatic_post/{tumor}/{seqtype}/realign/realign.bam",
        realign_bai="deepsomatic_post/{tumor}/{seqtype}/realign/realign.bam.bai",
    message:
        "--- Realigning DeepSomatic results for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess_realign", 8)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_realign", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_realign.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_realign.tsv"
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
        parsed_bed="deepsomatic_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
    output:
        pileup_tasks="deepsomatic_post/{tumor}/{seqtype}/pileup/workspace/pileup_tasks.txt",
    message:
        "--- Splitting DeepSomatic BED for pileup for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
        num_chunks=config.get("pileup_num_chunks", 16),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess_split", 1)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_split", 4000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_split_bed.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_split_bed.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/split_bed.sh \
            {input.parsed_bed} \
            {params.output_dir} \
            {output.pileup_tasks} \
            {params.num_chunks} &> {log}
        """

rule deepsomatic_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        bam_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(wc.tumor, wc.seqtype, wc.tumor),
        pileup_tasks="deepsomatic_post/{tumor}/{seqtype}/pileup/workspace/pileup_tasks.txt",
        reference_fa="deepsomatic_post/{tumor}/{seqtype}/reference.fa",
    output:
        pileup="deepsomatic_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz",
        pileup_tbi="deepsomatic_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz.tbi",
    message:
        "--- Running pileup for DeepSomatic postprocess {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
        no_baq=config.get("pileup_no_baq", "true"),
        max_depth=config.get("pileup_max_depth", 0),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess_pileup", 8)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_pileup", 32000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_pileup.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_pileup.tsv"
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
            {params.no_baq} \
            {params.max_depth} &> {log}
        """

rule deepsomatic_postprocess_haplotype:
    input:
        parsed_bed="deepsomatic_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        realign_bam="deepsomatic_post/{tumor}/{seqtype}/realign/realign.bam",
        pileup="deepsomatic_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz",
        kmer_ratio=lambda wc: "bam_refiner/{}/{}/{}_kmer_ratio.txt".format(wc.tumor, wc.seqtype, wc.tumor),
        hap1_list=lambda wc: "bam_refiner/{}/kmer/hap1_list.txt.gz".format(get_kmer_source(wc.tumor)),
        hap2_list=lambda wc: "bam_refiner/{}/kmer/hap2_list.txt.gz".format(get_kmer_source(wc.tumor)),
    output:
        haplotyped="deepsomatic_post/{tumor}/{seqtype}/{tumor}.haplotyped.bed",
    message:
        "--- Running haplotyping for DeepSomatic postprocess {wildcards.tumor} ({wildcards.seqtype})"
    params:
        tumor="{tumor}",
        output_dir="deepsomatic_post/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("deepsomatic_postprocess_haplotype", 4)
    resources:
        mem_mb=get_mem_mb("deepsomatic_postprocess_haplotype", 16000)
    log:
        "logs/deepsomatic_postprocess/{tumor}_{seqtype}_haplotype.log"
    benchmark:
        "benchmarks/deepsomatic_postprocess/{tumor}_{seqtype}_haplotype.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/haplotyping.sh \
            {params.tumor} \
            {params.output_dir} \
            {input.bam} \
            {input.kmer_ratio} \
            {input.hap1_list} \
            {input.hap2_list} &> {log}
        """
