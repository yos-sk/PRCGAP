# ====================================================================
# COPY NUMBER
# ====================================================================
#
# Split per haplotype so hap1 and hap2 run as concurrent cluster jobs instead
# of being looped inside one script:
#
#   copynumber_reference (once per run)
#     └► copynumber_ref_table (per hap) ──► copynumber_ref_table_final
#                                              └► copynumber_depth (per hap)
#                                                   └► copynumber_segment (per hap)
#                                                        └► copynumber_plot
#
# Only two steps genuinely need both haplotypes: the male sex-chromosome
# reconciliation in copynumber_ref_table_final, and the plot.
#
# Copy number uses a single seqtype: HiFi if available, else ONT
# (primary_paired_seqtype = the preferred seqtype present for both the tumor
# and its paired normal).

# The per-step resource keys fall back to the legacy `copynumber` key, so a
# config that only sets `copynumber` still sizes the split rules.
COPYNUMBER_REF_DIR = "copynumber/reference"
COPYNUMBER_REFERENCE = COPYNUMBER_REF_DIR + "/reference.fa"
COPYNUMBER_CHM13_LENGTHS = COPYNUMBER_REF_DIR + "/chm13_chrom_length.txt"


wildcard_constraints:
    hap="hap1|hap2",


def _cn_bam(sample, tumor):
    """Refined BAM of `sample` in the seqtype copy number analyses `tumor` with."""
    return "bam_refiner/{}/{}/{}_bam_refined.sorted.bam".format(
        sample, primary_paired_seqtype(tumor), sample)


def _cn_ref_table(tumor, hap):
    """Published reference table path for a tumor's haplotype."""
    return "copynumber/{}/output/{}.{}.ref.table".format(
        tumor, get_paired_normal(tumor), hap)


def _cn_ploidy_override(hap):
    """Manual ploidy for this haplotype, honoured only when BOTH are set.

    Mirrors the original script: a single-haplotype override is ignored so the
    two haplotypes are never scaled by a mix of manual and estimated ploidy.
    """
    hap1 = str(config.get("copynumber_ploidy_hap1", "") or "")
    hap2 = str(config.get("copynumber_ploidy_hap2", "") or "")
    if not (hap1 and hap2):
        return ""
    return hap1 if hap == "hap1" else hap2


# ====================================================================
# Staged reference (shared by every tumor and haplotype)
# ====================================================================

rule copynumber_reference:
    input:
        reference=config.get("chm13_fasta", ""),
    output:
        fasta=COPYNUMBER_REFERENCE,
        lengths=COPYNUMBER_CHM13_LENGTHS,
    message:
        "--- Staging the CHM13 reference for copy number analysis"
    params:
        sex=config.get("sex", "female"),
    threads:
        get_threads("copynumber_reference", 1)
    resources:
        mem_mb=get_mem_mb("copynumber_reference", 8000)
    log:
        "logs/copynumber/reference.log"
    benchmark:
        "benchmarks/copynumber/reference.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_reference.sh \
            "{params.sex}" \
            "{input.reference}" \
            "{output.fasta}" \
            "{output.lengths}" \
            "{SCRIPTS_DIR}" &> {log}
        """


# ====================================================================
# Contig ↔ CHM13 correspondence table, one job per haplotype
# ====================================================================

rule copynumber_ref_table:
    input:
        assembly=lambda wc: samples.loc[wc.tumor, "assembly_" + wc.hap],
        reference=COPYNUMBER_REFERENCE,
        chm13_censat=as_input(config.get("chm13_censat", "")),
    output:
        table="copynumber/{tumor}/workspace/ref_table_{hap}.raw",
        # postprocess_sex_chrom.py reads both PAFs for a male sample, so the
        # alignment has to be a declared output rather than scratch.
        paf="copynumber/{tumor}/workspace/{hap}/{hap}_ref.paf",
    message:
        "--- Reference table for {wildcards.tumor} {wildcards.hap}"
    params:
        hap="{hap}",
        work_dir="copynumber/{tumor}/workspace/{hap}",
        chm13_censat=config.get("chm13_censat", ""),
    threads:
        get_threads("copynumber_ref_table", get_threads("copynumber", 8))
    resources:
        mem_mb=get_mem_mb("copynumber_ref_table", get_mem_mb("copynumber", 48000))
    log:
        "logs/copynumber/{tumor}_ref_table_{hap}.log"
    benchmark:
        "benchmarks/copynumber/{tumor}_ref_table_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_ref_table.sh \
            "{params.hap}" \
            "{input.assembly}" \
            "{input.reference}" \
            "{params.chm13_censat}" \
            "{params.work_dir}" \
            "{output.table}" \
            "{SCRIPTS_DIR}" \
            "{threads}" &> {log}
        """


rule copynumber_ref_table_final:
    input:
        hap1="copynumber/{tumor}/workspace/ref_table_hap1.raw",
        hap2="copynumber/{tumor}/workspace/ref_table_hap2.raw",
        paf1="copynumber/{tumor}/workspace/hap1/hap1_ref.paf",
        paf2="copynumber/{tumor}/workspace/hap2/hap2_ref.paf",
    output:
        hap1="copynumber/{tumor}/output/{normal}.hap1.ref.table",
        hap2="copynumber/{tumor}/output/{normal}.hap2.ref.table",
    message:
        "--- Publishing reference tables for {wildcards.tumor}"
    params:
        sex=config.get("sex", "female"),
    wildcard_constraints:
        normal="|".join(re.escape(str(s)) for s in normals) or "$^",
    threads:
        get_threads("copynumber_ref_table_final", 1)
    resources:
        mem_mb=get_mem_mb("copynumber_ref_table_final", 8000)
    log:
        "logs/copynumber/{tumor}_{normal}_ref_table_final.log"
    benchmark:
        "benchmarks/copynumber/{tumor}_{normal}_ref_table_final.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_ref_table_final.sh \
            "{params.sex}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{output.hap1}" \
            "{output.hap2}" \
            "{SCRIPTS_DIR}" \
            "{input.paf1}" \
            "{input.paf2}" &> {log}
        """


# ====================================================================
# Depth ratios, one job per haplotype
# ====================================================================

rule copynumber_depth:
    input:
        tumor_bam=lambda wc: _cn_bam(wc.tumor, wc.tumor),
        normal_bam=lambda wc: _cn_bam(get_paired_normal(wc.tumor), wc.tumor),
        ref_table=lambda wc: _cn_ref_table(wc.tumor, wc.hap),
    output:
        tsv="copynumber/{tumor}/output/{tumor}.{hap}.copynumber.tsv",
    message:
        "--- Depth ratios for {wildcards.tumor} {wildcards.hap}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        hap="{hap}",
        work_dir="copynumber/{tumor}/workspace",
    threads:
        get_threads("copynumber_depth", get_threads("copynumber", 8))
    resources:
        mem_mb=get_mem_mb("copynumber_depth", get_mem_mb("copynumber", 16000))
    log:
        "logs/copynumber/{tumor}_depth_{hap}.log"
    benchmark:
        "benchmarks/copynumber/{tumor}_depth_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_depth.sh \
            "{params.tumor}" \
            "{params.normal}" \
            "{params.hap}" \
            "{input.tumor_bam}" \
            "{input.normal_bam}" \
            "{input.ref_table}" \
            "{params.work_dir}" \
            "{output.tsv}" \
            "{SCRIPTS_DIR}" \
            "{threads}" &> {log}
        """


# ====================================================================
# CBS segmentation + gap splitting, one job per haplotype
# ====================================================================

rule copynumber_segment:
    input:
        tsv="copynumber/{tumor}/output/{tumor}.{hap}.copynumber.tsv",
        ref_table=lambda wc: _cn_ref_table(wc.tumor, wc.hap),
    output:
        cbs="copynumber/{tumor}/output/{tumor}.{hap}.cbs.txt",
        split="copynumber/{tumor}/output/{tumor}.{hap}.cbs.split.txt",
        ploidy="copynumber/{tumor}/output/{tumor}.{hap}.ploidy",
    message:
        "--- CBS segmentation for {wildcards.tumor} {wildcards.hap}"
    params:
        tumor="{tumor}",
        hap="{hap}",
        output_dir="copynumber/{tumor}/output",
        ploidy=lambda wc: _cn_ploidy_override(wc.hap),
        binwidth=config.get("copynumber_binwidth", 0.05),
    threads:
        get_threads("copynumber_segment", 1)
    resources:
        mem_mb=get_mem_mb("copynumber_segment", get_mem_mb("copynumber", 8000))
    log:
        "logs/copynumber/{tumor}_segment_{hap}.log"
    benchmark:
        "benchmarks/copynumber/{tumor}_segment_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_segment.sh \
            "{params.tumor}" \
            "{params.hap}" \
            "{input.tsv}" \
            "{input.ref_table}" \
            "{params.output_dir}" \
            "{params.ploidy}" \
            "{params.binwidth}" \
            "{SCRIPTS_DIR}" &> {log}
        """


# ====================================================================
# Plot (needs both haplotypes)
# ====================================================================

rule copynumber_plot:
    input:
        tsv=expand("copynumber/{{tumor}}/output/{{tumor}}.{hap}.copynumber.tsv",
                   hap=["hap1", "hap2"]),
        split=expand("copynumber/{{tumor}}/output/{{tumor}}.{hap}.cbs.split.txt",
                     hap=["hap1", "hap2"]),
        ploidy=expand("copynumber/{{tumor}}/output/{{tumor}}.{hap}.ploidy",
                      hap=["hap1", "hap2"]),
        ref_tables=lambda wc: [_cn_ref_table(wc.tumor, h) for h in ("hap1", "hap2")],
        satellites=lambda wc: as_input(satellite_bed(wc.tumor, "hap1"),
                                       satellite_bed(wc.tumor, "hap2")),
        lengths=COPYNUMBER_CHM13_LENGTHS,
    output:
        "copynumber/{tumor}/output/{tumor}.copynumber.png"
    message:
        "--- Plotting copy number for {wildcards.tumor}"
    params:
        tumor="{tumor}",
        normal=lambda wc: get_paired_normal(wc.tumor),
        output_dir="copynumber/{tumor}/output",
        work_dir="copynumber/{tumor}/workspace",
        censat_bed=config.get("censat_bed", ""),
        hap1_satellite=lambda wc: satellite_bed(wc.tumor, "hap1"),
        hap2_satellite=lambda wc: satellite_bed(wc.tumor, "hap2"),
        binwidth=config.get("copynumber_binwidth", 0.05),
        hap1_label=config.get("copynumber_hap1_label", "Haplotype1"),
        hap2_label=config.get("copynumber_hap2_label", "Haplotype2"),
        plot_sex_chrom=str(config.get("copynumber_plot_sex_chrom", True)),
        chm13_censat=config.get("chm13_censat", ""),
    threads:
        get_threads("copynumber_plot", 1)
    resources:
        mem_mb=get_mem_mb("copynumber_plot", get_mem_mb("copynumber", 8000))
    log:
        "logs/copynumber/{tumor}_plot.log"
    benchmark:
        "benchmarks/copynumber/{tumor}_plot.tsv"
    singularity:
        config.get("singularity_images", {}).get("copynumber", "")
    shell:
        """
        /bin/bash {SCRIPTS_DIR}/copynumber/copynumber_plot.sh \
            "{params.tumor}" \
            "{params.normal}" \
            "{params.output_dir}" \
            "{params.work_dir}" \
            "{input.lengths}" \
            "{params.censat_bed}" \
            "{params.hap1_satellite}" \
            "{params.hap2_satellite}" \
            "{params.binwidth}" \
            "{params.hap1_label}" \
            "{params.hap2_label}" \
            "{params.plot_sex_chrom}" \
            "{params.chm13_censat}" \
            "{SCRIPTS_DIR}" &> {log}
        """
