# ====================================================================
# NANOMONSV PARSE
# ====================================================================

rule nanomonsv_parse_hifi:
    input:
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample),
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
        "--- Running nanomonsv parse for {wildcards.sample} HiFi data"
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
        bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.sample, wc.sample),
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
        "--- Running nanomonsv parse for {wildcards.sample} ONT data"
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
# NANOMONSV GET
# ====================================================================

rule nanomonsv_get_hifi:
    input:
        tumor_rearrangement=lambda wc: "nanomonsv/hifi/{}.rearrangement.sorted.bedpe.gz".format(wc.tumor),
        normal_rearrangement=lambda wc: "nanomonsv/hifi/{}.rearrangement.sorted.bedpe.gz".format(get_paired_normal(wc.tumor)),
        tumor_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        result="nanomonsv/hifi/{tumor}.nanomonsv.result.txt",
        supporting="nanomonsv/hifi/{tumor}.nanomonsv.supporting_read.txt"
    message:
        "--- Running nanomonsv get for {wildcards.tumor} HiFi data"
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
        tumor_rearrangement=lambda wc: "nanomonsv/ont/{}.rearrangement.sorted.bedpe.gz".format(wc.tumor),
        normal_rearrangement=lambda wc: "nanomonsv/ont/{}.rearrangement.sorted.bedpe.gz".format(get_paired_normal(wc.tumor)),
        tumor_bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
        normal_bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(get_paired_normal(wc.tumor), get_paired_normal(wc.tumor)),
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        result="nanomonsv/ont/{tumor}.nanomonsv.result.txt",
        supporting="nanomonsv/ont/{tumor}.nanomonsv.supporting_read.txt"
    message:
        "--- Running nanomonsv get for {wildcards.tumor} ONT data"
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
# NANOMONSV POSTPROCESS
# ====================================================================

rule nanomonsv_postprocess_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.result.txt",
        bam=lambda wc: "bam_refiner/{}/hifi/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt"
    message:
        "--- Running nanomonsv postprocess for {wildcards.tumor} HiFi data"
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
            {input.bam} \
            {SCRIPTS_DIR} &> {log}
        touch {output}
        """

rule nanomonsv_postprocess_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.result.txt",
        bam=lambda wc: "bam_refiner/{}/ont/{}_bam_refined.sorted.bam".format(wc.tumor, wc.tumor),
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt"
    message:
        "--- Running nanomonsv postprocess for {wildcards.tumor} ONT data"
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
# ASSEMBLY BWA INDEX (shared by insert_classify hifi/ont)
# ====================================================================
#
# bwa index writes .pac/.bwt/.sa next to the assembly fasta. If hifi and ont
# insert_classify rules both invoke `bwa index` concurrently they corrupt
# each other's intermediate files. Centralise indexing in a single rule and
# make insert_classify depend on its marker output. We use the kmer-source
# canonical sample pattern so each (hap1, hap2) pair is indexed only once.

rule assembly_bwa_index:
    input:
        hap1=lambda wc: samples.loc[wc.sample, "assembly_hap1"],
        hap2=lambda wc: samples.loc[wc.sample, "assembly_hap2"],
    output:
        marker=touch("bwa_index/{sample}/.indexed")
    message:
        "--- Building bwa index for assembly of {wildcards.sample}"
    threads:
        get_threads("assembly_bwa_index", 1)
    resources:
        mem_mb=get_mem_mb("assembly_bwa_index", 16000)
    log:
        "logs/bwa_index/{sample}.log"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        ( for asm in {input.hap1} {input.hap2}; do
              if [ ! -e "${{asm}}.sa" ]; then
                  bwa index "${{asm}}"
              fi
          done ) &> {log}
        """


# ====================================================================
# NANOMONSV INSERT CLASSIFY
# ====================================================================

rule nanomonsv_insert_classify_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        bwa_index_done=lambda wc: "bwa_index/{}/.indexed".format(get_kmer_source(wc.tumor)),
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt"
    message:
        "--- Running nanomonsv insert classification for {wildcards.tumor} HiFi data"
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
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.gtf_file} \
            {params.line1_bed} \
            {SCRIPTS_DIR} &> {log}
        """

rule nanomonsv_insert_classify_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt",
        assembly_hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        assembly_hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        bwa_index_done=lambda wc: "bwa_index/{}/.indexed".format(get_kmer_source(wc.tumor)),
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt"
    message:
        "--- Running nanomonsv insert classification for {wildcards.tumor} ONT data"
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
            {input.assembly_hap1} \
            {input.assembly_hap2} \
            {params.gtf_file} \
            {params.line1_bed} \
            {SCRIPTS_DIR} &> {log}
        """

# ====================================================================
# NANOMONSV CONNECT
# ====================================================================

rule nanomonsv_connect_hifi:
    input:
        result="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
        supporting="nanomonsv/hifi/{tumor}.nanomonsv.supporting_read.txt"
    output:
        "nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.connected.txt"
    message:
        "--- Running nanomonsv connect for {wildcards.tumor} HiFi data"
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
            {params.output_dir} \
            {SCRIPTS_DIR} &> {log}
        touch {output}
        """

rule nanomonsv_connect_ont:
    input:
        result="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt",
        supporting="nanomonsv/ont/{tumor}.nanomonsv.supporting_read.txt"
    output:
        "nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.connected.txt"
    message:
        "--- Running nanomonsv connect for {wildcards.tumor} ONT data"
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
            {params.output_dir} \
            {SCRIPTS_DIR} &> {log}
        touch {output}
        """

# ====================================================================
# NANOMONSV MERGE
# ====================================================================

rule nanomonsv_merge:
    input:
        hifi="nanomonsv/hifi/{tumor}.nanomonsv.new_result.sv_typed.txt",
        ont="nanomonsv/ont/{tumor}.nanomonsv.new_result.sv_typed.txt"
    output:
        "nanomonsv/{tumor}.nanomonsv.result.merged.txt"
    message:
        "--- Running nanomonsv merge for {wildcards.tumor}"
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
