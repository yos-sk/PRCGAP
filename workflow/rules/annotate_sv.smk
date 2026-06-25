# ====================================================================
# SV ANNOTATION (per tumor, per seqtype) + reclassification
# ====================================================================
#
# Chain (per seqtype):
#   prep_sv  ──►  coordconv_sv_grch38 ─┐
#            └──► coordconv_sv_chm13  ─┤
#                                       ├──► annotate_sv_main  ──► reclassify_sv
#            └──► (PASS-filtered TSV) ──┘
#
# Every rule runs inside the single `annotation` singularity image
# (pysam + samtools + coordconv + transanno). Snakemake's `singularity:`
# directive wraps the shell automatically — no `singularity exec` calls
# inside the shell wrappers.
#
# coordconv rules only fire when the corresponding chain file is
# configured (chain_to_grch38 / chain_to_chm13). When unset, the
# downstream `annotate_sv_main` simply skips the conv step.

ANNOT_DIR = os.path.join(SCRIPTS_DIR, "annotate")


def _opt_path(key):
    return config.get(key, "") or ""


# ====================================================================
# Shared liftoff-derived gene BED.
#
# Reference convention:
#   SV       → `liftoff.bed.gz` (4-col-ish BED indexed with `tabix -p bed`)
#   SNV/INDEL → `liftoff.gff.gz` (full GFF indexed with `tabix -p gff`)
#
# We accept only the GFF from the user (`gff_file`) and derive the SV
# BED on the fly using the awk pipeline from PRCGAP-paper's assembly
# annotation step:
#   gene rows → drop pseudogenes → require gene_name → strip quotes/semis
#   → 7-col BED (chr, start-1, end, gene_name, strand, gene_id,
#   gene_biotype)
# Then sort + bgzip + `tabix -p bed`.
# ====================================================================

LIFTOFF_GENE_BED = "annotate_common/liftoff.gene.bed.gz"


def _liftoff_bed_input(_wildcards):
    return [LIFTOFF_GENE_BED] if _opt_path("gff_file") else []


def _liftoff_bed_param(_wildcards):
    return LIFTOFF_GENE_BED if _opt_path("gff_file") else ""


rule gff_to_bed:
    input:
        gff=_opt_path("gff_file") or [],
    output:
        bed=LIFTOFF_GENE_BED,
        tbi=LIFTOFF_GENE_BED + ".tbi",
    message:
        "--- liftoff GFF → sorted+tabix-indexed gene BED"
    threads:
        get_threads("gff_to_bed", 1)
    resources:
        mem_mb=get_mem_mb("gff_to_bed", 4000)
    log:
        "logs/annotate_common/gff_to_bed.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        # Awk pipeline mirrors the assembly-annotation step used to build
        # PRCGAP-paper's `*.liftoff.bed.gz`. Field offsets ($10, $14,
        # $18) correspond to Ensembl-liftoff GFF attribute ordering after
        # quote/semicolon stripping.
        """
        ( zcat -f {input.gff} \
            | awk '{{if ($3 == "gene") print}}' \
            | grep -v pseudogene \
            | grep gene_name \
            | awk '{{gsub(/[";]/, ""); print $1 "\\t" $4 - 1 "\\t" $5 "\\t" $14 "\\t" $7 "\\t" $10 "\\t" $18}}' \
            | sort -k1,1 -k2,2n \
            | bgzip -c > {output.bed}
          tabix -p bed {output.bed}
        ) &> {log}
        """


def _coordconv_input(wildcards, ref):
    """Return the coordconv BED path if the chain is configured, else []."""
    chain_key = "chain_to_grch38" if ref == "GRCh38" else "chain_to_chm13"
    if _opt_path(chain_key):
        return (
            "annotate_sv/" + wildcards.tumor + "/" + wildcards.seqtype + "/workspace/"
            + wildcards.tumor + ".coordconv_" + ref + ".bed"
        )
    return []


def _coordconv_param(wildcards, ref):
    """Return the coordconv BED path as a string for the shell wrapper.

    Mirrors _coordconv_input but always returns a string ("" when unset)
    so the shell wrapper can detect "skip this step" via empty arg.
    """
    chain_key = "chain_to_grch38" if ref == "GRCh38" else "chain_to_chm13"
    if _opt_path(chain_key):
        return (
            "annotate_sv/" + wildcards.tumor + "/" + wildcards.seqtype + "/workspace/"
            + wildcards.tumor + ".coordconv_" + ref + ".bed"
        )
    return ""


# ====================================================================
# PREP: filter + extract breakpoint BED
# ====================================================================

rule prep_sv:
    input:
        nanomonsv="nanomonsv/{seqtype}/{tumor}.nanomonsv.new_result.sv_typed.insert_classified.txt",
    output:
        filt_pass="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.filt.pass.txt",
        bp_bed="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.bp.bed",
    message:
        "--- Preparing SV annotation inputs for {wildcards.tumor} {wildcards.seqtype}"
    params:
        sample="{tumor}",
        work_dir="annotate_sv/{tumor}/{seqtype}/workspace",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("prep_sv", 1)
    resources:
        mem_mb=get_mem_mb("prep_sv", 8000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_prep.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOT_DIR}/prep_sv.sh \
            "{input.nanomonsv}" \
            "{params.work_dir}" \
            "{params.sample}" \
            "{ANNOT_DIR}" &> {log}
        """


# ====================================================================
# COORDCONV: assembly-coord BED → GRCh38 / CHM13 (optional)
# ====================================================================

rule coordconv_sv_grch38:
    input:
        bp_bed="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.bp.bed",
    output:
        "annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.coordconv_GRCh38.bed"
    message:
        "--- coordconv (GRCh38) for {wildcards.tumor} {wildcards.seqtype}"
    params:
        chain=_opt_path("chain_to_grch38"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_sv", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_sv", 8000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_coordconv_GRCh38.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        coordconv -b {input.bp_bed} -c {params.chain} > {output} 2> {log}
        """


rule coordconv_sv_chm13:
    input:
        bp_bed="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.bp.bed",
    output:
        "annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.coordconv_chm13.bed"
    message:
        "--- coordconv (chm13) for {wildcards.tumor} {wildcards.seqtype}"
    params:
        chain=_opt_path("chain_to_chm13"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_sv", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_sv", 8000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_coordconv_chm13.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        coordconv -b {input.bp_bed} -c {params.chain} > {output} 2> {log}
        """


# ====================================================================
# MAIN: gene/rmsk/size/conv/kmer/cen/segdup/other/misassembly/gnomad
# ====================================================================

rule annotate_sv_main:
    input:
        pass_txt="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.filt.pass.txt",
        nanomonsv_other=lambda wc: (
            "nanomonsv/ont/{}.nanomonsv.new_result.sv_typed.insert_classified.txt".format(wc.tumor)
            if wc.seqtype == "hifi"
            else "nanomonsv/hifi/{}.nanomonsv.new_result.sv_typed.insert_classified.txt".format(wc.tumor)
        ),
        support_reads="nanomonsv/{seqtype}/{tumor}.nanomonsv.supporting_read.txt",
        kmer_ratio=lambda wc: "bam_refiner/{}/{}/{}_kmer_ratio.txt".format(wc.tumor, wc.seqtype, wc.tumor),
        hap1=lambda wc: samples.loc[wc.tumor, "assembly_hap1"],
        hap2=lambda wc: samples.loc[wc.tumor, "assembly_hap2"],
        grch38_bed=lambda wc: _coordconv_input(wc, "GRCh38"),
        chm13_bed=lambda wc: _coordconv_input(wc, "chm13"),
        liftoff_bed=lambda wc: _liftoff_bed_input(wc),
    output:
        "annotate_sv/{tumor}/{seqtype}/{tumor}.PRCGAP.nanomonsv_results.annotated.txt"
    message:
        "--- Running SV annotation for {wildcards.tumor} {wildcards.seqtype}"
    params:
        sample="{tumor}",
        output_dir="annotate_sv/{tumor}/{seqtype}",
        work_dir="annotate_sv/{tumor}/{seqtype}/workspace",
        liftoff_bed=lambda wc: _liftoff_bed_param(wc),
        cgc=_opt_path("cancer_gene_census_tsv"),
        rmsk=_opt_path("repeat_masker_bed"),
        censat=_opt_path("censat_bed"),
        segdup=_opt_path("segdup_bed"),
        misa1=_opt_path("misassembly_hap1_bed"),
        misa2=_opt_path("misassembly_hap2_bed"),
        grch38_bed=lambda wc: _coordconv_param(wc, "GRCh38"),
        chm13_bed=lambda wc: _coordconv_param(wc, "chm13"),
        gnomad=_opt_path("gnomad_bed"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("annotate_sv", 1)
    resources:
        mem_mb=get_mem_mb("annotate_sv", 16000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOT_DIR}/annotate_sv.sh \
            "{params.sample}" \
            "{input.pass_txt}" \
            "{input.nanomonsv_other}" \
            "{input.support_reads}" \
            "{input.kmer_ratio}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{params.output_dir}" \
            "{params.work_dir}" \
            "{ANNOT_DIR}" \
            "{params.liftoff_bed}" \
            "{params.cgc}" \
            "{params.rmsk}" \
            "{params.censat}" \
            "{params.segdup}" \
            "{params.misa1}" \
            "{params.misa2}" \
            "{params.grch38_bed}" \
            "{params.chm13_bed}" \
            "{params.gnomad}" &> {log}
        """


# ====================================================================
# RECLASSIFY SV TYPE (per tumor; combines hifi+ont, requires copynumber ref.table)
# ====================================================================

rule reclassify_sv:
    input:
        hifi_annotated="annotate_sv/{tumor}/hifi/{tumor}.PRCGAP.nanomonsv_results.annotated.txt",
        ont_annotated="annotate_sv/{tumor}/ont/{tumor}.PRCGAP.nanomonsv_results.annotated.txt",
        copynumber_dir=lambda wc: "copynumber/{}/output".format(wc.tumor),
    output:
        hifi="annotate_sv/{tumor}/{tumor}.hifi.PRCGAP.nanomonsv_results.reclassified.txt",
        ont="annotate_sv/{tumor}/{tumor}.ont.PRCGAP.nanomonsv_results.reclassified.txt",
    message:
        "--- Reclassifying SV types for {wildcards.tumor}"
    params:
        normal=lambda wc: get_paired_normal(wc.tumor),
    threads:
        get_threads("reclassify_sv", 1)
    resources:
        mem_mb=get_mem_mb("reclassify_sv", 16000)
    log:
        "logs/annotate_sv/{tumor}_reclassify.log"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        ( python3 {ANNOT_DIR}/reclassify_sv_type.py \
              -i {input.hifi_annotated} \
              -o {output.hifi} \
              -r {input.copynumber_dir}/{params.normal}.hap1.ref.table {input.copynumber_dir}/{params.normal}.hap2.ref.table
          python3 {ANNOT_DIR}/reclassify_sv_type.py \
              -i {input.ont_annotated} \
              -o {output.ont} \
              -r {input.copynumber_dir}/{params.normal}.hap1.ref.table {input.copynumber_dir}/{params.normal}.hap2.ref.table
        ) &> {log}
        """
