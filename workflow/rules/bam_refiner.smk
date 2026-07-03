# ====================================================================
# BAM REFINER
# ====================================================================
#
# bam_refiner is split into two phases:
#   1. bam_refiner_kmer  — extract haplotype-specific k-mers from the
#      assembly. Depends only on the assembly, so it runs once per sample
#      and is shared across HiFi/ONT.
#   2. bam_refiner_hifi / bam_refiner_ont — map reads against the assembly
#      and refine the alignment using the precomputed k-mer files.

rule bam_refiner_kmer:
    input:
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        hap1_kmer="bam_refiner/{sample}/kmer/hap1_cnt_kmerposition.bed.gz",
        hap1_kmer_tbi="bam_refiner/{sample}/kmer/hap1_cnt_kmerposition.bed.gz.tbi",
        hap2_kmer="bam_refiner/{sample}/kmer/hap2_cnt_kmerposition.bed.gz",
        hap2_kmer_tbi="bam_refiner/{sample}/kmer/hap2_cnt_kmerposition.bed.gz.tbi",
        hap1_list="bam_refiner/{sample}/kmer/hap1_list.txt.gz",
        hap2_list="bam_refiner/{sample}/kmer/hap2_list.txt.gz",
    message:
        "--- Extracting haplotype-specific k-mers for {wildcards.sample}"
    params:
        sample="{sample}",
        output_dir="bam_refiner/{sample}/kmer"
    threads:
        get_threads("bam_refiner_kmer", 8)
    resources:
        mem_mb=get_mem_mb("bam_refiner_kmer", 32000)
    log:
        "logs/bam_refiner/{sample}_kmer.log"
    singularity:
        config.get("singularity_images", {}).get("bam_refiner", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/bam_refiner/kmer_extract.sh \
            {params.sample} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.output_dir} \
            {threads} \
            {resources.mem_mb} &> {log}
        """

rule bam_refiner_hifi:
    input:
        files=lambda wc: get_sample_files_list(wc.sample, "hifi"),
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
        hap1_kmer=lambda wc: "bam_refiner/{}/kmer/hap1_cnt_kmerposition.bed.gz".format(get_kmer_source(wc.sample)),
        hap2_kmer=lambda wc: "bam_refiner/{}/kmer/hap2_cnt_kmerposition.bed.gz".format(get_kmer_source(wc.sample)),
        hap1_list=lambda wc: "bam_refiner/{}/kmer/hap1_list.txt.gz".format(get_kmer_source(wc.sample)),
        hap2_list=lambda wc: "bam_refiner/{}/kmer/hap2_list.txt.gz".format(get_kmer_source(wc.sample)),
    output:
        bam="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam",
        bai="bam_refiner/{sample}/hifi/{sample}_bam_refined.sorted.bam.bai",
        kmer_ratio="bam_refiner/{sample}/hifi/{sample}_kmer_ratio.txt"
    message:
        "--- Refining HiFi BAM for {wildcards.sample}"
    params:
        sample="{sample}",
        output_dir="bam_refiner/{sample}/hifi",
        kmer_dir=lambda wc: "bam_refiner/{}/kmer".format(get_kmer_source(wc.sample)),
        input_files=lambda wc: ",".join(get_sample_files_list(wc.sample, "hifi"))
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
        /bin/bash {SCRIPTS_DIR}/bam_refiner/refine.sh \
            {params.sample} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.input_files} \
            {params.output_dir} \
            {threads} \
            hifi \
            {params.kmer_dir} &> {log}
        """

rule bam_refiner_ont:
    input:
        files=lambda wc: get_sample_files_list(wc.sample, "ont"),
        assembly_hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
        hap1_kmer=lambda wc: "bam_refiner/{}/kmer/hap1_cnt_kmerposition.bed.gz".format(get_kmer_source(wc.sample)),
        hap2_kmer=lambda wc: "bam_refiner/{}/kmer/hap2_cnt_kmerposition.bed.gz".format(get_kmer_source(wc.sample)),
        hap1_list=lambda wc: "bam_refiner/{}/kmer/hap1_list.txt.gz".format(get_kmer_source(wc.sample)),
        hap2_list=lambda wc: "bam_refiner/{}/kmer/hap2_list.txt.gz".format(get_kmer_source(wc.sample)),
    output:
        bam="bam_refiner/{sample}/ont/{sample}_bam_refined.sorted.bam",
        bai="bam_refiner/{sample}/ont/{sample}_bam_refined.sorted.bam.bai",
        kmer_ratio="bam_refiner/{sample}/ont/{sample}_kmer_ratio.txt"
    message:
        "--- Refining ONT BAM for {wildcards.sample}"
    params:
        sample="{sample}",
        output_dir="bam_refiner/{sample}/ont",
        kmer_dir=lambda wc: "bam_refiner/{}/kmer".format(get_kmer_source(wc.sample)),
        input_files=lambda wc: ",".join(get_sample_files_list(wc.sample, "ont"))
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
        /bin/bash {SCRIPTS_DIR}/bam_refiner/refine.sh \
            {params.sample} \
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.input_files} \
            {params.output_dir} \
            {threads} \
            ont \
            {params.kmer_dir} &> {log}
        """
