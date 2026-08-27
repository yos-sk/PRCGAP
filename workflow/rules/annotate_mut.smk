# ====================================================================
# SNV / INDEL ANNOTATION (per tumor, per tool, per seqtype, per mode)
# ====================================================================
#
# Mirrors the per-tool snv/annotate_snv.sh and indel/annotate_indel.sh
# scripts under PRCGAP-paper/analysis/point_mutation_analysis/scripts/
# annotation/. The `compare_current_reference` /
# `filter_diff_references` germline-comparison steps are still deferred.
# `check_homozygous` is replaced by `check_unassigned`: a somatic mutation is
# essentially never homozygous and what looks homozygous is LOH, so the
# homozygous relabeling is gone, but its other half -- dropping the unassigned
# copy of a variant a haplotype-assigned record already covers at the same
# reference coordinate -- is kept. Coordinate annotation is kept:
#
#   prep_mut ──► coordconv (snv) ──┐
#            └─► bed2vcf ─► transanno_liftvcf (indel) ──┐
#                                                       ├──► annotate_mut_main
#                                                       │     (add_lift_coords →
#                                                       │      gene → rmsk → size →
#                                                       │      misa → cen → segdup →
#                                                       │      other →
#                                                       │      check_unassigned →
#                                                       │      gnomad)
#
# coordconv liftover is applied to SNV only — INDEL needs proper
# VCF-aware liftover that transforms ref/alt alleles, which transanno
# provides.
#
# Point-mutation calling runs per available sequencing type, so every path
# below carries a {seqtype} (hifi / ont) segment. Both callers share the same
# postprocess directory layout: <tool>_post/{tumor}/{seqtype}/... .


def _mut_other_vcf_path(wildcards):
    """The bgzipped+tabix-indexed VCF for the `other` cross-check, or None.

    The cross-check compares the two callers WITHIN the same seqtype, so it
    only applies when config mutation_caller selects both. With a single caller
    the `point_mutation_other` column is left empty rather than pulling the
    unselected caller's whole chain into the DAG.

    The raw `deepsomatic/{tumor}/{seqtype}` / `clairs/{tumor}/{seqtype}`
    directory outputs aren't tracked at the file level, so we always go
    through the postprocess + index pipeline:
      `<tool>_post/{tumor}/{seqtype}/output.vcf.gz` ──► `mut_vcf_index` rule
      ──► `<tool>_post/{tumor}/{seqtype}/output.tabix.vcf.gz` (+ `.tbi`)
    Both clairs and deepsomatic carry SNV+INDEL together at this stage
    (clairs_postprocess_prepare concatenates output.vcf + indel.vcf;
    deepsomatic_postprocess_prepare copies output.vcf as-is).
    """
    if len(mutation_callers()) < 2:
        return None
    other_tool = "deepsomatic" if wildcards.tool == "clairs" else "clairs"
    return "{}_post/{}/{}/output.tabix.vcf.gz".format(
        other_tool, wildcards.tumor, wildcards.seqtype)


def _mut_other_vcf(wildcards):
    """Input-list form of _mut_other_vcf_path: empty when there is no other caller."""
    path = _mut_other_vcf_path(wildcards)
    return [path] if path else []


def _per_tool_indel_vcf(wildcards):
    """The per-tool tabix-indexed somatic VCF whose INDEL records bed2vcf
    joins against. Same indexed-VCF target used for the `other` check.
    """
    return "{}_post/{}/{}/output.tabix.vcf.gz".format(
        wildcards.tool, wildcards.tumor, wildcards.seqtype)


def _coordconv_snv_input(wildcards, ref):
    if chain_file(wildcards.tumor, ref):
        return (
            "annotate_snv/" + wildcards.tumor + "/" + wildcards.tool + "/"
            + wildcards.seqtype + "/workspace/"
            + wildcards.tumor + "." + wildcards.tool + ".snv.coordconv_" + ref + ".bed"
        )
    return []


def _liftvcf_indel_input(wildcards, ref):
    if chain_file(wildcards.tumor, ref):
        return (
            "annotate_indel/" + wildcards.tumor + "/" + wildcards.tool + "/"
            + wildcards.seqtype + "/workspace/"
            + wildcards.tumor + "." + wildcards.tool + ".indel.liftvcf_" + ref + ".vcf.gz"
        )
    return []


def _lift_input(wildcards, ref):
    if wildcards.mode == "snv":
        return _coordconv_snv_input(wildcards, ref)
    return _liftvcf_indel_input(wildcards, ref)


def _lift_param(wildcards, ref):
    paths = _lift_input(wildcards, ref)
    if isinstance(paths, list):
        return ""
    return paths


# ====================================================================
# VCF INDEX: re-bgzip + tabix-index the postprocess output VCF.
#
# clairs_postprocess_prepare emits a regular-gzip concatenation of two
# gzip members (not bgzip), and deepsomatic_postprocess_prepare just
# `cp`s the raw VCF without a `.tbi`. Both downstream consumers
# (`bed2vcf.py` and `add_annotation_mut.py other`) use
# `pysam.TabixFile`, which needs bgzip + tbi. This rule normalizes
# both into a single tabix-ready form.
# ====================================================================

rule mut_vcf_index:
    input:
        vcf="{tool}_post/{tumor}/{seqtype}/output.vcf.gz",
    output:
        vcf="{tool}_post/{tumor}/{seqtype}/output.tabix.vcf.gz",
        tbi="{tool}_post/{tumor}/{seqtype}/output.tabix.vcf.gz.tbi",
    message:
        "--- bgzip+tabix index {wildcards.tool} VCF for {wildcards.tumor} ({wildcards.seqtype})"
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("mut_vcf_index", 1)
    resources:
        mem_mb=get_mem_mb("mut_vcf_index", 4000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_vcf_index.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_vcf_index.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        # `clairs_postprocess_prepare` appends INDEL records after SNV
        # records without re-sorting, so we sort body lines by chrom/pos
        # before bgzipping. deepsomatic_postprocess_prepare's output is
        # already sorted but re-sorting it is idempotent.
        """
        ( zcat {input.vcf} > {output.vcf}.unsorted
          ( grep '^#' {output.vcf}.unsorted
            grep -v '^#' {output.vcf}.unsorted | sort -S $(( {resources.mem_mb} / 2 < 8192 ? {resources.mem_mb} / 2 : 8192 ))M --parallel={threads} -k1,1 -k2,2n
          ) | bgzip -c > {output.vcf}
          tabix -p vcf {output.vcf}
          rm {output.vcf}.unsorted
        ) &> {log}
        """


# ====================================================================
# PREP: split haplotyped.bed into SNV / INDEL filtered BED (14-column,
# no header — exactly the input format expected by coordconv (snv) and
# bed2vcf.py (indel)).
# ====================================================================

rule prep_mut:
    input:
        haplotyped=lambda wc: "{}_post/{}/{}/{}.haplotyped.bed".format(wc.tool, wc.tumor, wc.seqtype, wc.tumor),
    output:
        prep_bed="annotate_{mode}/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.{mode}.bed",
    message:
        "--- Preparing {wildcards.mode} BED for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    params:
        sample="{tumor}",
        tool="{tool}",
        mode="{mode}",
        work_dir="annotate_{mode}/{tumor}/{tool}/{seqtype}/workspace",
    wildcard_constraints:
        mode="snv|indel",
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("prep_mut", 1)
    resources:
        mem_mb=get_mem_mb("prep_mut", 8000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_{mode}_prep.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_{mode}_prep.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOT_DIR}/prep_mut.sh \
            "{input.haplotyped}" \
            "{params.work_dir}" \
            "{params.sample}" \
            "{params.tool}" \
            "{params.mode}" &> {log}
        """


# ====================================================================
# SNV: coordconv liftover (GRCh38 / CHM13). Mirrors snv/annotate_snv.sh
# lines 37-45.
# ====================================================================

rule coordconv_snv_grch38:
    input:
        bed="annotate_snv/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.snv.bed",
        chain=lambda wc: as_input(chain_file(wc.tumor, "GRCh38")),
    output:
        "annotate_snv/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.snv.coordconv_GRCh38.bed"
    message:
        "--- coordconv (GRCh38) SNV for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    params:
        chain=lambda wc: chain_file(wc.tumor, "GRCh38"),
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_mut", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_mut", 8000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_snv_coordconv_GRCh38.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_snv_coordconv_GRCh38.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        coordconv -b {input.bed} -c {params.chain} > {output} 2> {log}
        """


rule coordconv_snv_chm13:
    input:
        bed="annotate_snv/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.snv.bed",
        chain=lambda wc: as_input(chain_file(wc.tumor, "chm13")),
    output:
        "annotate_snv/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.snv.coordconv_chm13.bed"
    message:
        "--- coordconv (chm13) SNV for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    params:
        chain=lambda wc: chain_file(wc.tumor, "chm13"),
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_mut", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_mut", 8000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_snv_coordconv_chm13.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_snv_coordconv_chm13.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        coordconv -b {input.bed} -c {params.chain} > {output} 2> {log}
        """


# ====================================================================
# INDEL: bed2vcf → samtools faidx (concat hap1+hap2) → transanno liftvcf
# Mirrors indel/annotate_indel.sh lines 40-89.
# ====================================================================

rule bed2vcf_indel:
    input:
        bed="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.bed",
        vcf=_per_tool_indel_vcf,
    output:
        "annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.vcf"
    message:
        "--- bed2vcf INDEL for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("bed2vcf_mut", 1)
    resources:
        mem_mb=get_mem_mb("bed2vcf_mut", 8000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_indel_bed2vcf.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_indel_bed2vcf.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        python3 {ANNOT_DIR}/bed2vcf.py \
            -i {input.bed} \
            -v {input.vcf} \
            -o {output} &> {log}
        """


rule indel_hap_reference:
    """Concatenate hap1+hap2 fasta and index with faidx. Required by
    transanno liftvcf as the *source* reference (the assembly the VCF
    was called against).
    """
    input:
        hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
    output:
        fa="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/reference.fa",
        fai="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/reference.fa.fai",
    message:
        "--- assembling concat hap1+hap2 reference for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("indel_hap_reference", 1)
    resources:
        mem_mb=get_mem_mb("indel_hap_reference", 4000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_indel_hap_reference.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_indel_hap_reference.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        cat {input.hap1} {input.hap2} > {output.fa} 2> {log}
        samtools faidx {output.fa} 2>> {log}
        """


rule liftvcf_indel_grch38:
    input:
        vcf="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.vcf",
        ref_fa="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/reference.fa",
        chain=lambda wc: as_input(chain_file(wc.tumor, "GRCh38")),
    output:
        vcf="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.liftvcf_GRCh38.vcf.gz",
        rej="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.liftvcf_GRCh38.reject.vcf.gz",
    message:
        "--- transanno liftvcf (GRCh38) INDEL for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    params:
        chain=lambda wc: chain_file(wc.tumor, "GRCh38"),
        query=_opt_path("grch38_fasta"),
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("liftvcf_mut", 1)
    resources:
        mem_mb=get_mem_mb("liftvcf_mut", 16000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_indel_liftvcf_GRCh38.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_indel_liftvcf_GRCh38.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        transanno liftvcf -m \
            --chain {params.chain} \
            -o {output.vcf} \
            --reference {input.ref_fa} \
            --query {params.query} \
            --vcf {input.vcf} \
            --fail {output.rej} &> {log}
        """


rule liftvcf_indel_chm13:
    input:
        vcf="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.vcf",
        ref_fa="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/reference.fa",
        chain=lambda wc: as_input(chain_file(wc.tumor, "chm13")),
    output:
        vcf="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.liftvcf_chm13.vcf.gz",
        rej="annotate_indel/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.indel.liftvcf_chm13.reject.vcf.gz",
    message:
        "--- transanno liftvcf (chm13) INDEL for {wildcards.tumor} {wildcards.tool} ({wildcards.seqtype})"
    params:
        chain=lambda wc: chain_file(wc.tumor, "chm13"),
        query=_opt_path("chm13_fasta"),
    wildcard_constraints:
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("liftvcf_mut", 1)
    resources:
        mem_mb=get_mem_mb("liftvcf_mut", 16000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_indel_liftvcf_chm13.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_indel_liftvcf_chm13.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        transanno liftvcf -m \
            --chain {params.chain} \
            -o {output.vcf} \
            --reference {input.ref_fa} \
            --query {params.query} \
            --vcf {input.vcf} \
            --fail {output.rej} &> {log}
        """


# ====================================================================
# MAIN: add_lift_coords → gene → rmsk → size → misa → cen → segdup →
# other → gnomad
# ====================================================================

rule annotate_mut_main:
    input:
        prep_bed="annotate_{mode}/{tumor}/{tool}/{seqtype}/workspace/{tumor}.{tool}.{mode}.bed",
        other_vcf=_mut_other_vcf,
        hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        grch38_lift=lambda wc: _lift_input(wc, "GRCh38"),
        chm13_lift=lambda wc: _lift_input(wc, "chm13"),
        gff=lambda wc: as_input(liftoff_gff(wc.tumor)),
    output:
        "annotate_{mode}/{tumor}/{tool}/{seqtype}/{tumor}.{tool}.{mode}.annotated.txt"
    message:
        "--- Annotating {wildcards.mode} ({wildcards.tool}) for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        sample="{tumor}",
        tool="{tool}",
        mode="{mode}",
        output_dir="annotate_{mode}/{tumor}/{tool}/{seqtype}",
        work_dir="annotate_{mode}/{tumor}/{tool}/{seqtype}/workspace",
        other_vcf=lambda wc: _mut_other_vcf_path(wc) or "",
        gff_file=lambda wc: liftoff_gff(wc.tumor),
        cgc=_opt_path("cancer_gene_census_tsv"),
        cmrg=_opt_path("cmrg_gene_tsv"),
        mane=_opt_path("mane_summary"),
        rmsk=_opt_path("repeat_masker_bed"),
        censat=_opt_path("censat_bed"),
        segdup=_opt_path("segdup_bed"),
        misa1=_opt_path("misassembly_hap1_bed"),
        misa2=_opt_path("misassembly_hap2_bed"),
        gnomad=_opt_path("gnomad_vcf"),
        grch38_lift=lambda wc: _lift_param(wc, "GRCh38"),
        chm13_lift=lambda wc: _lift_param(wc, "chm13"),
    wildcard_constraints:
        mode="snv|indel",
        tool="clairs|deepsomatic",
        seqtype="hifi|ont",
    threads:
        get_threads("annotate_mut", 1)
    resources:
        mem_mb=get_mem_mb("annotate_mut", 16000)
    log:
        "logs/annotate_mut/{tumor}_{tool}_{seqtype}_{mode}.log"
    benchmark:
        "benchmarks/annotate_mut/{tumor}_{tool}_{seqtype}_{mode}.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOT_DIR}/annotate_mut.sh \
            "{params.sample}" \
            "{params.tool}" \
            "{params.mode}" \
            "{input.prep_bed}" \
            "{params.other_vcf}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{params.output_dir}" \
            "{params.work_dir}" \
            "{ANNOT_DIR}" \
            "{params.gff_file}" \
            "{params.cgc}" \
            "{params.cmrg}" \
            "{params.mane}" \
            "{params.rmsk}" \
            "{params.censat}" \
            "{params.segdup}" \
            "{params.misa1}" \
            "{params.misa2}" \
            "{params.gnomad}" \
            "{params.grch38_lift}" \
            "{params.chm13_lift}" &> {log}
        """
