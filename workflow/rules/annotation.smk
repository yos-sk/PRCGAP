# ====================================================================
# ASSEMBLY ANNOTATION (on by default; switchable per step)
# ====================================================================
#
# The minimum annotation set PRCGAP needs to interpret its own calls, built
# from the assembly plus the CHM13/GRCh38 references so a run no longer has to
# import these from the assembly_workflow repo:
#
#   run_dna_brnn     ──► dna_brnn (per hap)          satellite BED per haplotype
#                                                    (copynumber masking /
#                                                     cenSat fallback)
#   run_liftoff      ──► liftoff_reference
#                          └► liftoff_hap (per hap) ──► liftoff_merge
#                                                    gene GFF + GTF
#                                                    (SV/SNV/INDEL gene
#                                                     annotation, nanomonsv
#                                                     insert_classify)
#   run_chain_files  ──► prepare_mask_regions
#                          └► chain_files_hap (per hap × ref) ──► chain_files_merge
#                                                    assembly → CHM13/GRCh38
#                                                    chains (coordconv /
#                                                     transanno liftover)
#
# Every per-haplotype step is its own rule so hap1 and hap2 are submitted as
# concurrent cluster jobs rather than looped inside one script; the chain step
# splits per reference as well, giving four independent alignments.
#
# Each rule is keyed on the assembly source sample (see annotation_src in
# commons.smk), so a tumor/normal pair sharing an assembly annotates once.
# Steps not covered here — sedef, cenSat, RepeatMasker, misassembly — remain
# config-supplied inputs.

import re

ANNOTATION_SCRIPTS = os.path.join(SCRIPTS_DIR, "annotation")

LIFTOFF_REF_DIR = MASKED_REF_DIR + "/liftoff"


wildcard_constraints:
    # Keeps annotation/{asmsrc}/... from also matching annotation/references/...
    asmsrc="|".join(re.escape(str(s)) for s in samples.index),
    hap="hap1|hap2",
    ref="GRCh38|chm13",


def _assembly_of(src, hap):
    return samples.loc[src, "assembly_" + hap]


def _masked_reference(ref):
    """The masked reference the configured sex selects."""
    suffix = "_noY" if config.get("sex", "female") == "female" else ""
    return "{}/{}.masked{}.fa".format(MASKED_REF_DIR, ref, suffix)


# ====================================================================
# dna-brnn: alpha-satellite / HSat regions, one job per haplotype
# ====================================================================

rule dna_brnn:
    input:
        fasta=lambda wc: _assembly_of(wc.asmsrc, wc.hap),
    output:
        bed=ANNOTATION_DIR + "/{asmsrc}/dna_nn/{asmsrc}.{hap}_dna-brnn.bed.gz",
        tbi=ANNOTATION_DIR + "/{asmsrc}/dna_nn/{asmsrc}.{hap}_dna-brnn.bed.gz.tbi",
    message:
        "--- dna-brnn satellite annotation for {wildcards.asmsrc} {wildcards.hap}"
    params:
        sample="{asmsrc}",
        hap="{hap}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/dna_nn",
    threads:
        get_threads("dna_brnn", 8)
    resources:
        mem_mb=get_mem_mb("dna_brnn", 8000)
    log:
        "logs/annotation/{asmsrc}_dna_brnn_{hap}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_dna_brnn_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("dna_nn", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/dna_brnn.sh \
            "{params.sample}" \
            "{params.hap}" \
            "{input.fasta}" \
            "{params.output_dir}" \
            "{threads}" &> {log}
        """


# ====================================================================
# liftoff: GRCh38 gene annotation → assembly coordinates
# ====================================================================
#
# The reference is staged once (chrY dropped for a female sample) and shared by
# both haplotype jobs, so the 3 GB filter is not repeated per haplotype.

rule liftoff_reference:
    input:
        grch38=config.get("grch38_fasta", "") or [],
        grch38_gtf=config.get("grch38_gtf", "") or [],
    output:
        fasta=LIFTOFF_REF_DIR + "/GRCh38.fa",
        gtf=LIFTOFF_REF_DIR + "/GRCh38.gtf",
    message:
        "--- Staging the GRCh38 reference for liftoff"
    params:
        sex=config.get("sex", "female"),
    threads:
        get_threads("liftoff_reference", 1)
    resources:
        mem_mb=get_mem_mb("liftoff_reference", 8000)
    log:
        "logs/annotation/liftoff_reference.log"
    benchmark:
        "benchmarks/annotation/liftoff_reference.tsv"
    singularity:
        config.get("singularity_images", {}).get("liftoff", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/liftoff_reference.sh \
            "{params.sex}" \
            "{input.grch38}" \
            "{input.grch38_gtf}" \
            "{output.fasta}" \
            "{output.gtf}" &> {log}
        """


rule liftoff_hap:
    input:
        fasta=lambda wc: _assembly_of(wc.asmsrc, wc.hap),
        grch38=LIFTOFF_REF_DIR + "/GRCh38.fa",
        grch38_gtf=LIFTOFF_REF_DIR + "/GRCh38.gtf",
    output:
        gff=ANNOTATION_DIR + "/{asmsrc}/liftoff/workspace/{asmsrc}.{hap}.liftoff.gff",
    message:
        "--- liftoff gene annotation for {wildcards.asmsrc} {wildcards.hap}"
    params:
        sample="{asmsrc}",
        hap="{hap}",
    threads:
        get_threads("liftoff", 8)
    resources:
        mem_mb=get_mem_mb("liftoff", 64000)
    log:
        "logs/annotation/{asmsrc}_liftoff_{hap}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_liftoff_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("liftoff", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/liftoff_hap.sh \
            "{params.sample}" \
            "{params.hap}" \
            "{input.fasta}" \
            "{input.grch38}" \
            "{input.grch38_gtf}" \
            "{output.gff}" \
            "{threads}" &> {log}
        """


rule liftoff_merge:
    input:
        hap1=ANNOTATION_DIR + "/{asmsrc}/liftoff/workspace/{asmsrc}.hap1.liftoff.gff",
        hap2=ANNOTATION_DIR + "/{asmsrc}/liftoff/workspace/{asmsrc}.hap2.liftoff.gff",
    output:
        gff=ANNOTATION_DIR + "/{asmsrc}/liftoff/{asmsrc}.liftoff.gff.gz",
        gff_tbi=ANNOTATION_DIR + "/{asmsrc}/liftoff/{asmsrc}.liftoff.gff.gz.tbi",
        gtf=ANNOTATION_DIR + "/{asmsrc}/liftoff/{asmsrc}.liftoff.gtf.gz",
    message:
        "--- Merging liftoff haplotype GFFs for {wildcards.asmsrc}"
    params:
        sample="{asmsrc}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/liftoff",
    threads:
        get_threads("liftoff_merge", 1)
    resources:
        mem_mb=get_mem_mb("liftoff_merge", 16000)
    log:
        "logs/annotation/{asmsrc}_liftoff_merge.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_liftoff_merge.tsv"
    singularity:
        config.get("singularity_images", {}).get("liftoff", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/liftoff_merge.sh \
            "{params.sample}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{params.output_dir}" &> {log}
        """


# ====================================================================
# chain files: masked references, then assembly → CHM13/GRCh38 chains
# ====================================================================
#
# The masked references depend only on the reference inputs, so they are built
# once per run and shared by every haplotype's alignment job.

rule prepare_mask_regions:
    input:
        chm13=config.get("chm13_fasta", "") or [],
        chm13_censat=config.get("chm13_censat", "") or [],
        grch38=config.get("grch38_fasta", "") or [],
        grch38_centromeres=config.get("grch38_centromeres", "") or [],
        grch38_exclusions=config.get("grch38_exclusions", "") or [],
    output:
        # Only the two _masked_reference() will ask for. Both variants used to be
        # built and the sex picked one pair at use time, leaving the other pair
        # -- 6 GB on a whole genome -- written and never read.
        chm13_masked=_masked_reference("chm13"),
        chm13_masked_fai=_masked_reference("chm13") + ".fai",
        grch38_masked=_masked_reference("GRCh38"),
        grch38_masked_fai=_masked_reference("GRCh38") + ".fai",
    message:
        "--- Masking satellite/centromere regions in CHM13 and GRCh38"
    params:
        output_dir=MASKED_REF_DIR,
        work_dir=MASKED_REF_DIR + "/workspace",
        drop_chry="true" if config.get("sex", "female") == "female" else "false",
    threads:
        get_threads("prepare_mask_regions", 1)
    resources:
        mem_mb=get_mem_mb("prepare_mask_regions", 10240)
    log:
        "logs/annotation/prepare_mask_regions.log"
    benchmark:
        "benchmarks/annotation/prepare_mask_regions.tsv"
    singularity:
        config.get("singularity_images", {}).get("chain_files", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/prepare_mask_regions.sh \
            "{input.chm13}" \
            "{input.chm13_censat}" \
            "{input.grch38}" \
            "{input.grch38_centromeres}" \
            "{input.grch38_exclusions}" \
            "{params.output_dir}" \
            "{params.work_dir}" \
            "{ANNOTATION_SCRIPTS}" \
            "{params.drop_chry}" &> {log}
        """


rule chain_files_hap:
    input:
        fasta=lambda wc: _assembly_of(wc.asmsrc, wc.hap),
        reference=lambda wc: _masked_reference(wc.ref),
    output:
        chain=(ANNOTATION_DIR
               + "/{asmsrc}/chain_files/workspace/{asmsrc}_{hap}.{ref}.inverted.chain"),
    message:
        "--- Chain {wildcards.asmsrc} {wildcards.hap} → {wildcards.ref}"
    params:
        sample="{asmsrc}",
        hap="{hap}",
        ref="{ref}",
    threads:
        get_threads("make_chain_files", 8)
    resources:
        mem_mb=get_mem_mb("make_chain_files", 48000)
    log:
        "logs/annotation/{asmsrc}_chain_{hap}_{ref}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_chain_{hap}_{ref}.tsv"
    singularity:
        config.get("singularity_images", {}).get("chain_files", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/chain_files_hap.sh \
            "{params.sample}" \
            "{params.hap}" \
            "{input.fasta}" \
            "{params.ref}" \
            "{input.reference}" \
            "{output.chain}" \
            "{threads}" &> {log}
        """


rule chain_files_merge:
    input:
        hap1=(ANNOTATION_DIR
              + "/{asmsrc}/chain_files/workspace/{asmsrc}_hap1.{ref}.inverted.chain"),
        hap2=(ANNOTATION_DIR
              + "/{asmsrc}/chain_files/workspace/{asmsrc}_hap2.{ref}.inverted.chain"),
    output:
        chain=ANNOTATION_DIR + "/{asmsrc}/chain_files/{asmsrc}_to_{ref}.chain",
    message:
        "--- Merging haplotype chains for {wildcards.asmsrc} → {wildcards.ref}"
    threads:
        get_threads("chain_files_merge", 1)
    resources:
        mem_mb=get_mem_mb("chain_files_merge", 8000)
    log:
        "logs/annotation/{asmsrc}_chain_merge_{ref}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_chain_merge_{ref}.tsv"
    singularity:
        config.get("singularity_images", {}).get("chain_files", "")
    shell:
        """
        cat {input.hap1} {input.hap2} > {output.chain} 2> {log}
        """


# ====================================================================
# LINE-1: full-length young elements, one job per haplotype
# ====================================================================
#
# nanomonsv insert_classify needs the source elements a transduction can come
# from, which is the full-length L1HS / L1PA2-L1PA5 set, not a general L1
# annotation. Running RepeatMasker for that costs ~240 core-hours per sample;
# aligning L1.3 and typing the hits with Dfam's subunit models reproduces the
# same list in about four minutes. See workspace/repeat_bench/README.md,
# Experiment 8, for the validation and the parameter choices.

# resource/line1 at the repo root, built by resource/scripts/
# build_line1_resources.sh; workflow.basedir is workflow/.
LINE1_RESOURCES = os.path.join(workflow.basedir, os.pardir, "resource", "line1")


rule line1_hap:
    input:
        fasta=lambda wc: _assembly_of(wc.asmsrc, wc.hap),
    output:
        bed=ANNOTATION_DIR + "/{asmsrc}/line1/{asmsrc}.{hap}.LINE1.bed",
    message:
        "--- LINE-1 elements for {wildcards.asmsrc} {wildcards.hap}"
    params:
        sample="{asmsrc}",
        hap="{hap}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/line1",
        resource_dir=LINE1_RESOURCES,
    threads:
        get_threads("line1", 8)
    resources:
        mem_mb=get_mem_mb("line1", 16000)
    log:
        "logs/annotation/{asmsrc}_line1_{hap}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_line1_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("nanomonsv", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/line1_hap.sh \
            "{params.sample}" \
            "{params.hap}" \
            "{input.fasta}" \
            "{params.output_dir}" \
            "{threads}" \
            "{SCRIPTS_DIR}" \
            "{params.resource_dir}" &> {log}
        """


rule line1_merge:
    input:
        hap1=ANNOTATION_DIR + "/{asmsrc}/line1/{asmsrc}.hap1.LINE1.bed",
        hap2=ANNOTATION_DIR + "/{asmsrc}/line1/{asmsrc}.hap2.LINE1.bed",
    output:
        bed=ANNOTATION_DIR + "/{asmsrc}/line1/{asmsrc}.LINE1.bed.gz",
        tbi=ANNOTATION_DIR + "/{asmsrc}/line1/{asmsrc}.LINE1.bed.gz.tbi",
    message:
        "--- Merging LINE-1 haplotype BEDs for {wildcards.asmsrc}"
    params:
        sample="{asmsrc}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/line1",
    threads:
        get_threads("line1_merge", 1)
    resources:
        mem_mb=get_mem_mb("line1_merge", 8000)
    log:
        "logs/annotation/{asmsrc}_line1_merge.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_line1_merge.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/line1_merge.sh \
            "{params.sample}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{params.output_dir}" &> {log}
        """


# ====================================================================
# Simple repeats: ULTRA with the period capped at 10, one job per haplotype
# ====================================================================
#
# nanomonsv get uses this to drop indel-like SVs whose breakpoints both sit in
# the same tandem repeat, so the period cap matches the MaxPeriod RepeatMasker
# gives TRF in its Simple_repeat stage rather than annotating satellite arrays.


rule simple_repeat_hap:
    input:
        fasta=lambda wc: _assembly_of(wc.asmsrc, wc.hap),
    output:
        bed=ANNOTATION_DIR + "/{asmsrc}/simple_repeat/{asmsrc}.{hap}.simple_repeats.bed",
    message:
        "--- Tandem repeat annotation for {wildcards.asmsrc} {wildcards.hap}"
    params:
        sample="{asmsrc}",
        hap="{hap}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/simple_repeat",
        min_copies=config.get("simple_repeat_min_copies", 4),
    threads:
        get_threads("simple_repeat", 8)
    resources:
        mem_mb=get_mem_mb("simple_repeat", 32000)
    log:
        "logs/annotation/{asmsrc}_simple_repeat_{hap}.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_simple_repeat_{hap}.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/simple_repeat_hap.sh \
            "{params.sample}" \
            "{params.hap}" \
            "{input.fasta}" \
            "{params.output_dir}" \
            "{threads}" \
            "{params.min_copies}" &> {log}
        """


rule simple_repeat_merge:
    input:
        hap1=ANNOTATION_DIR + "/{asmsrc}/simple_repeat/{asmsrc}.hap1.simple_repeats.bed",
        hap2=ANNOTATION_DIR + "/{asmsrc}/simple_repeat/{asmsrc}.hap2.simple_repeats.bed",
    output:
        bed=ANNOTATION_DIR + "/{asmsrc}/simple_repeat/{asmsrc}.simple_repeats.bed.gz",
        tbi=ANNOTATION_DIR + "/{asmsrc}/simple_repeat/{asmsrc}.simple_repeats.bed.gz.tbi",
    message:
        "--- Merging tandem repeat BEDs for {wildcards.asmsrc}"
    params:
        sample="{asmsrc}",
        output_dir=ANNOTATION_DIR + "/{asmsrc}/simple_repeat",
    threads:
        get_threads("simple_repeat_merge", 1)
    resources:
        mem_mb=get_mem_mb("simple_repeat_merge", 8000)
    log:
        "logs/annotation/{asmsrc}_simple_repeat_merge.log"
    benchmark:
        "benchmarks/annotation/{asmsrc}_simple_repeat_merge.tsv"
    singularity:
        config.get("singularity_images", {}).get("annotation", "")
    shell:
        """
        /bin/bash {ANNOTATION_SCRIPTS}/simple_repeat_merge.sh \
            "{params.sample}" \
            "{input.hap1}" \
            "{input.hap2}" \
            "{params.output_dir}" &> {log}
        """
