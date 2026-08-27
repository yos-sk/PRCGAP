# ====================================================================
# CLAIRS (contig scatter -> gather)
# ====================================================================
#
# ClairS runs per available sequencing type ({seqtype} = hifi / ont), scattered
# over contig chunks: every contig >= config["caller_solo_contig_min_bp"] gets its own
# job and the remaining small contigs share one. Per-chunk threads and memory
# come from config["resources"].
#
# The per-seqtype tumor/normal BAMs come from bam_refiner. The ClairS platform
# model is selected via config["clairs_model"] (default hifi_sequel2); valid
# values: hifi_sequel2, hifi_revio, ont_r10_dorado_sup_5khz_ssrs,
# ont_r10_dorado_sup_4khz.

rule clairs_scatter_setup:
    input:
        tumor_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        tumor_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(wc.tumor, wc.seqtype, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), wc.seqtype, get_paired_normal(wc.tumor)),
        normal_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(get_paired_normal(wc.tumor), wc.seqtype, get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        reference_fa="clairs_scatter/{tumor}/{seqtype}/reference.fa",
        reference_fai="clairs_scatter/{tumor}/{seqtype}/reference.fa.fai",
        tumor_bam="clairs_scatter/{tumor}/{seqtype}/tumor.bam",
        tumor_bai="clairs_scatter/{tumor}/{seqtype}/tumor.bam.bai",
        normal_bam="clairs_scatter/{tumor}/{seqtype}/normal.bam",
        normal_bai="clairs_scatter/{tumor}/{seqtype}/normal.bam.bai",
        sanitized="clairs_scatter/{tumor}/{seqtype}/sanitized.flag",
    message:
        "--- Preparing ClairS scatter inputs for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        output_dir="clairs_scatter/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_scatter_setup", 8)
    resources:
        mem_mb=get_mem_mb("clairs_scatter_setup", 16000)
    log:
        "logs/clairs/{tumor}_{seqtype}_scatter_setup.log"
    benchmark:
        "benchmarks/clairs/{tumor}_{seqtype}_scatter_setup.tsv"
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
checkpoint clairs_chunks:
    input:
        reference_fai="clairs_scatter/{tumor}/{seqtype}/reference.fa.fai",
    output:
        chunk_dir=directory("clairs_scatter/{tumor}/{seqtype}/chunks"),
    message:
        "--- Building ClairS contig chunks for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        min_bp=config.get("caller_solo_contig_min_bp", 1000000),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_chunks", 1)
    resources:
        mem_mb=get_mem_mb("clairs_chunks", 8000)
    log:
        "logs/clairs/{tumor}_{seqtype}_chunks.log"
    benchmark:
        "benchmarks/clairs/{tumor}_{seqtype}_chunks.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/caller/make_chunks.sh \
            {input.reference_fai} \
            {output.chunk_dir} \
            {params.min_bp} &> {log}
        """

rule clairs_chunk:
    input:
        chunk_dir=lambda wc: checkpoints.clairs_chunks.get(
            tumor=wc.tumor, seqtype=wc.seqtype).output.chunk_dir,
        reference_fa="clairs_scatter/{tumor}/{seqtype}/reference.fa",
        tumor_bam="clairs_scatter/{tumor}/{seqtype}/tumor.bam",
        normal_bam="clairs_scatter/{tumor}/{seqtype}/normal.bam",
    output:
        vcf="clairs_scatter/{tumor}/{seqtype}/out/{chunk}/output.vcf.gz",
        indel="clairs_scatter/{tumor}/{seqtype}/out/{chunk}/indel.vcf.gz",
    message:
        "--- Running ClairS chunk {wildcards.chunk} for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        output_dir="clairs_scatter/{tumor}/{seqtype}/out/{chunk}",
        ctg_file="clairs_scatter/{tumor}/{seqtype}/chunks/{chunk}.ctg",
        model=config.get("clairs_model", "hifi_sequel2"),
    wildcard_constraints:
        seqtype="hifi|ont",
        chunk=r"\d+",
    threads:
        get_threads("clairs", 8)
    resources:
        mem_mb=get_mem_mb("clairs", 32000)
    log:
        "logs/clairs/{tumor}_{seqtype}_chunk_{chunk}.log"
    benchmark:
        "benchmarks/clairs/{tumor}_{seqtype}_chunk_{chunk}.tsv"
    singularity:
        config.get("singularity_images", {}).get("clairs", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/clairs/chunk.sh \
            {input.tumor_bam} \
            {input.normal_bam} \
            {input.reference_fa} \
            {params.ctg_file} \
            {params.output_dir} \
            {threads} \
            {params.model} &> {log}
        """

def _clairs_chunk_files(wc, basename):
    chunk_dir = checkpoints.clairs_chunks.get(
        tumor=wc.tumor, seqtype=wc.seqtype).output.chunk_dir
    chunks = sorted(glob_wildcards(os.path.join(chunk_dir, "{chunk}.ctg")).chunk)
    return expand(
        "clairs_scatter/{tumor}/{seqtype}/out/{chunk}/" + basename,
        tumor=wc.tumor, seqtype=wc.seqtype, chunk=chunks)

rule clairs_merge:
    input:
        vcfs=lambda wc: _clairs_chunk_files(wc, "output.vcf.gz"),
        indels=lambda wc: _clairs_chunk_files(wc, "indel.vcf.gz"),
        reference_fai="clairs_scatter/{tumor}/{seqtype}/reference.fa.fai",
        sanitized="clairs_scatter/{tumor}/{seqtype}/sanitized.flag",
    output:
        vcf="clairs/{tumor}/{seqtype}/output.vcf.gz",
        vcf_tbi="clairs/{tumor}/{seqtype}/output.vcf.gz.tbi",
        indel="clairs/{tumor}/{seqtype}/indel.vcf.gz",
        indel_tbi="clairs/{tumor}/{seqtype}/indel.vcf.gz.tbi",
    message:
        "--- Merging ClairS chunks for {wildcards.tumor} ({wildcards.seqtype})"
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_merge", 1)
    resources:
        mem_mb=get_mem_mb("clairs_merge", 16000),
        disk_mb=8000,
    log:
        "logs/clairs/{tumor}_{seqtype}_merge.log"
    benchmark:
        "benchmarks/clairs/{tumor}_{seqtype}_merge.tsv"
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
        /bin/bash {SCRIPTS_DIR}/caller/merge_vcf.sh \
            {input.reference_fai} \
            {input.sanitized} \
            {output.indel} \
            $(( {resources.mem_mb} / 2 )) \
            {input.indels} &>> {log}
        """

# ====================================================================
# CLAIRS POSTPROCESS PIPELINE
# ====================================================================

rule clairs_postprocess_prepare:
    input:
        clairs_vcf="clairs/{tumor}/{seqtype}/output.vcf.gz",
        clairs_indel="clairs/{tumor}/{seqtype}/indel.vcf.gz",
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
        clairs_dir="clairs/{tumor}/{seqtype}",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess", 1)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_prepare.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_prepare.tsv"
    singularity:
        config.get("singularity_images", {}).get("point_mutation_postprocess", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/mutation_postprocess/clairs_prepare.sh \
            {params.clairs_dir} \
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
        get_threads("clairs_postprocess", 1)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_parse_vcf.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_parse_vcf.tsv"
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
        get_threads("clairs_postprocess_realign", 8)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_realign", 16000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_realign.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_realign.tsv"
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
        num_chunks=config.get("pileup_num_chunks", 16),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_split", 1)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_split", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_split_bed.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_split_bed.tsv"
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

rule clairs_postprocess_pileup:
    input:
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        bam_bai=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam.bai".format(wc.tumor, wc.seqtype, wc.tumor),
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
        no_baq=config.get("pileup_no_baq", "true"),
        max_depth=config.get("pileup_max_depth", 0),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("clairs_postprocess_pileup", 8)
    resources:
        mem_mb=get_mem_mb("clairs_postprocess_pileup", 32000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_pileup.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_pileup.tsv"
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

rule clairs_postprocess_haplotype:
    input:
        parsed_bed="clairs_post/{tumor}/{seqtype}/realign/parsed_vcf.bed",
        bam=lambda wc: "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(wc.tumor, wc.seqtype, wc.tumor),
        realign_bam="clairs_post/{tumor}/{seqtype}/realign/realign.bam",
        pileup="clairs_post/{tumor}/{seqtype}/pileup/{tumor}_pileup.bed.gz",
        kmer_ratio=lambda wc: "bam_refiner/{}/{}/{}_kmer_ratio.txt".format(wc.tumor, wc.seqtype, wc.tumor),
        hap1_list=lambda wc: "bam_refiner/{}/kmer/hap1_list.txt.gz".format(get_kmer_source(wc.tumor)),
        hap2_list=lambda wc: "bam_refiner/{}/kmer/hap2_list.txt.gz".format(get_kmer_source(wc.tumor)),
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
        mem_mb=get_mem_mb("clairs_postprocess_haplotype", 64000)
    log:
        "logs/clairs_postprocess/{tumor}_{seqtype}_haplotype.log"
    benchmark:
        "benchmarks/clairs_postprocess/{tumor}_{seqtype}_haplotype.tsv"
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
