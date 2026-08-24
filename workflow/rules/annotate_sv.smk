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


def _nanomonsv_other_input(wildcards):
    """The complementary seqtype's insert_classified result, for the HiFi↔ONT
    cross-check in annotate_sv_main. Returns [] when the other seqtype is not
    available for this tumor (HiFi-only / ONT-only), skipping the cross-check.
    """
    other = "ont" if wildcards.seqtype == "hifi" else "hifi"
    if other in paired_seqtypes(wildcards.tumor):
        return "nanomonsv/{}/{}.nanomonsv.new_result.sv_typed.insert_classified.txt".format(
            other, wildcards.tumor)
    return []


def _nanomonsv_other_param(wildcards):
    path = _nanomonsv_other_input(wildcards)
    return path if isinstance(path, str) else ""


# ====================================================================
# Shared liftoff-derived gene BED.
#
# Reference convention:
#   SV       → `liftoff.bed.gz` (4-col-ish BED indexed with `tabix -p bed`)
#   SNV/INDEL → `liftoff.gff.gz` (full GFF indexed with `tabix -p gff`)
#
# The GFF is either supplied by the user (`gff_file`) or built by the
# workflow's own liftoff rule (`run_liftoff`); either way the SV BED is
# derived on the fly using the awk pipeline from PRCGAP-paper's assembly
# annotation step:
#   gene rows → drop pseudogenes → require gene_name → strip quotes/semis
#   → 7-col BED (chr, start-1, end, gene_name, strand, gene_id,
#   gene_biotype)
# Then sort + bgzip + `tabix -p bed`.
#
# The BED is keyed on the assembly source sample so a workflow-generated,
# per-assembly GFF stays distinguishable; with a config-supplied GFF the
# samples sharing an assembly simply produce identical copies.
# ====================================================================

def _liftoff_gene_bed(src):
    return "annotate_common/{}/liftoff.gene.bed.gz".format(src)


def _liftoff_bed_input(wildcards):
    src = annotation_src(wildcards.tumor)
    return [_liftoff_gene_bed(src)] if liftoff_gff_src(src) else []


def _liftoff_bed_param(wildcards):
    src = annotation_src(wildcards.tumor)
    return _liftoff_gene_bed(src) if liftoff_gff_src(src) else ""


rule gff_to_bed:
    input:
        gff=lambda wc: as_input(liftoff_gff_src(wc.asmsrc)),
    output:
        bed="annotate_common/{asmsrc}/liftoff.gene.bed.gz",
        tbi="annotate_common/{asmsrc}/liftoff.gene.bed.gz.tbi",
    message:
        "--- liftoff GFF → sorted+tabix-indexed gene BED ({wildcards.asmsrc})"
    threads:
        get_threads("gff_to_bed", 1)
    resources:
        mem_mb=get_mem_mb("gff_to_bed", 4000)
    log:
        "logs/annotate_common/{asmsrc}_gff_to_bed.log"
    benchmark:
        "benchmarks/annotate_common/{asmsrc}_gff_to_bed.tsv"
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
            | sort -S $(( {resources.mem_mb} / 2 < 8192 ? {resources.mem_mb} / 2 : 8192 ))M --parallel={threads} -k1,1 -k2,2n \
            | bgzip -c > {output.bed}
          tabix -p bed {output.bed}
        ) &> {log}
        """


def _coordconv_input(wildcards, ref):
    """Return the coordconv BED path if a chain is available, else []."""
    if chain_file(wildcards.tumor, ref):
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
    if chain_file(wildcards.tumor, ref):
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
    benchmark:
        "benchmarks/annotate_sv/{tumor}_{seqtype}_prep.tsv"
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
        chain=lambda wc: as_input(chain_file(wc.tumor, "GRCh38")),
    output:
        "annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.coordconv_GRCh38.bed"
    message:
        "--- coordconv (GRCh38) for {wildcards.tumor} {wildcards.seqtype}"
    params:
        chain=lambda wc: chain_file(wc.tumor, "GRCh38"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_sv", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_sv", 8000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_coordconv_GRCh38.log"
    benchmark:
        "benchmarks/annotate_sv/{tumor}_{seqtype}_coordconv_GRCh38.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        coordconv -b {input.bp_bed} -c {params.chain} > {output} 2> {log}
        """


rule coordconv_sv_chm13:
    input:
        bp_bed="annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.bp.bed",
        chain=lambda wc: as_input(chain_file(wc.tumor, "chm13")),
    output:
        "annotate_sv/{tumor}/{seqtype}/workspace/{tumor}.coordconv_chm13.bed"
    message:
        "--- coordconv (chm13) for {wildcards.tumor} {wildcards.seqtype}"
    params:
        chain=lambda wc: chain_file(wc.tumor, "chm13"),
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("coordconv_sv", 1)
    resources:
        mem_mb=get_mem_mb("coordconv_sv", 8000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_coordconv_chm13.log"
    benchmark:
        "benchmarks/annotate_sv/{tumor}_{seqtype}_coordconv_chm13.tsv"
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
        nanomonsv_other=lambda wc: _nanomonsv_other_input(wc),
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
        nanomonsv_other=lambda wc: _nanomonsv_other_param(wc),
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
    benchmark:
        "benchmarks/annotate_sv/{tumor}_{seqtype}.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOT_DIR}/annotate_sv.sh \
            "{params.sample}" \
            "{input.pass_txt}" \
            "{params.nanomonsv_other}" \
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
# RECLASSIFY SV TYPE (per tumor, per seqtype; requires copynumber ref.table)
# ====================================================================

rule reclassify_sv:
    input:
        annotated="annotate_sv/{tumor}/{seqtype}/{tumor}.PRCGAP.nanomonsv_results.annotated.txt",
        # Depend on the copynumber rule's declared output (the .png) so Snakemake
        # knows how to build it; the ref.table files are produced alongside it in
        # the same output dir (see params.copynumber_dir).
        copynumber_png=lambda wc: "copynumber/{}/output/{}.copynumber.png".format(wc.tumor, wc.tumor),
    output:
        "annotate_sv/{tumor}/{tumor}.{seqtype}.PRCGAP.nanomonsv_results.reclassified.txt",
    message:
        "--- Reclassifying SV types for {wildcards.tumor} ({wildcards.seqtype})"
    params:
        normal=lambda wc: get_paired_normal(wc.tumor),
        copynumber_dir="copynumber/{tumor}/output",
    wildcard_constraints:
        seqtype="hifi|ont",
    threads:
        get_threads("reclassify_sv", 1)
    resources:
        mem_mb=get_mem_mb("reclassify_sv", 16000)
    log:
        "logs/annotate_sv/{tumor}_{seqtype}_reclassify.log"
    benchmark:
        "benchmarks/annotate_sv/{tumor}_{seqtype}_reclassify.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        python3 {ANNOT_DIR}/reclassify_sv_type.py \
            -i {input.annotated} \
            -o {output} \
            -r {params.copynumber_dir}/{params.normal}.hap1.ref.table {params.copynumber_dir}/{params.normal}.hap2.ref.table &> {log}
        """
