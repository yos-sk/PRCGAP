#!/bin/bash
# Regenerate the configs and runners for the HG008 chr20 test case, one per
# sequencing-type pattern. Run from test/HG008/:
#
#   bash test_configure.sh          # writes config/ and run_{hifi,ont,both}.sh
#   bash run_hifi.sh                # HiFi-only pair
#   bash run_ont.sh                 # ONT-only pair
#   bash run_both.sh                # HiFi + ONT
#
# Test order (project rule): this HG008 case must pass before the H2009
# full-pipeline case under test/H2009/ is started.
#
# Each pattern gets its own sample sheet, config, runner and output directory:
#   config/samples_<pattern>.tsv, config/config_<pattern>.yaml,
#   run_<pattern>.sh, ./<pattern>/
#
# Set MUTATION_CALLER to exercise the other point-mutation caller path:
#   MUTATION_CALLER=both bash test_configure.sh
#
# By default the pre-extracted chr20 annotation under resources/annotation/ is
# used and PRCGAP's own annotation steps are off, so this case tests the calling
# pipeline against a fixed annotation. Set RUN_ANNOTATION=1 to instead exercise
# dna-brnn / liftoff / chain-file generation inside the workflow:
#   RUN_ANNOTATION=1 bash test_configure.sh
# That variant additionally needs the GTF / centromeres / exclusions fetched by
# ../../resource/scripts/download_reference.sh.
#
# Prerequisites (produced by helpers under resources/scripts/):
#   - bash resources/scripts/download_data.sh
#       resources/reads/{ont,hifi}/HG008{T,N}.chr20.*.bam
#       resources/asm/HG008N.hap{1,2}.chr20.fa
#   - bash resources/scripts/extract_haplotypes.sh
#       resources/annotation/HG008N.*.{bed,gtf,gff}.gz
#       resources/annotation/HG008N.to_{grch38,chm13}.chain
#   - bash resources/scripts/download_annotation.sh
#       resources/annotation/cmrg_genes.list
#       resources/annotation/gencode.v46.basic.annotation.chr20.mane.transcript.bed.gz
#       resources/annotation/gnomad.v4.1.sv.sites.bed.gz{,.tbi}
#       resources/annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz{,.tbi}
#   - bash ../../resource/scripts/download_reference.sh   (writes ../../resource/reference/)
#       ../../resource/reference/chm13v2.0_maskedY_rCRS.fa{,.fai}
#       ../../resource/reference/GRCh38.d1.vd1.fa{,.fai}
#
# Manual prerequisites (no helper script):
#   - resources/annotation/cancer_gene_census.tsv  (requires COSMIC login; download manually)
#   - ../../images/*.sif singularity images pre-staged

set -euo pipefail

MUTATION_CALLER=${MUTATION_CALLER:-deepsomatic}
RUN_ANNOTATION=${RUN_ANNOTATION:-0}

if [ "${RUN_ANNOTATION}" = "1" ]; then
    # Build the annotation in-workflow. The satellite / liftoff / chain flags
    # below are still passed but ignored, since the generated files win.
    ANNOTATION_ARGS=(--run-dna-brnn --run-liftoff --run-chain-files
                     --run-line1 --run-simple-repeat
                     --grch38-gtf         ../../resource/reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf
                     --grch38-centromeres ../../resource/reference/centromeres.txt.gz
                     --grch38-exclusions  ../../resource/reference/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed)
else
    ANNOTATION_ARGS=(--no-run-dna-brnn --no-run-liftoff --no-run-chain-files
                     --no-run-line1 --no-run-simple-repeat)
fi

# $HOME may be a symlink (here /home/<user> -> /lustre10/home/<user>). Config
# paths keep the user-facing form, but rules that resolve symlinks (e.g.
# assembly_bwa_index's `readlink -f`) produce the target form, so both have to
# be visible inside the container.
SING_BIND="$HOME"
HOME_REAL=$(readlink -f "$HOME")
if [ "${HOME_REAL}" != "$HOME" ]; then
    SING_BIND="$HOME,${HOME_REAL}"
fi

# The sample sheet's read paths are resolved relative to the run's output
# directory (snakemake runs with cwd=output_dir), which is ./<pattern> here, so
# they are one level deeper than this script.
for PATTERN in hifi ont both; do
    case ${PATTERN} in
        hifi) SEQ_ARGS=(--tumor-hifi ../resources/reads/hifi/HG008T.chr20.hifi.bam
                        --normal-hifi ../resources/reads/hifi/HG008N.chr20.hifi.bam) ;;
        ont)  SEQ_ARGS=(--tumor-ont  ../resources/reads/ont/HG008T.chr20.ont.bam
                        --normal-ont  ../resources/reads/ont/HG008N.chr20.ont.bam) ;;
        both) SEQ_ARGS=(--tumor-hifi ../resources/reads/hifi/HG008T.chr20.hifi.bam
                        --normal-hifi ../resources/reads/hifi/HG008N.chr20.hifi.bam
                        --tumor-ont  ../resources/reads/ont/HG008T.chr20.ont.bam
                        --normal-ont  ../resources/reads/ont/HG008N.chr20.ont.bam) ;;
    esac

    echo "--- configuring HG008 pattern: ${PATTERN}"
    mkdir -p "./${PATTERN}"

    # --no-check-exists: configure does not require the multi-GB read/asm files
    # to be staged yet. --no-absolutize: keep the resource paths relative.
    python3 ../../set_sample_sheet.py \
        --tumor  HG008T \
        --normal HG008N \
        "${SEQ_ARGS[@]}" \
        --assembly-hap1 ../resources/asm/HG008N.hap1.chr20.fa \
        --assembly-hap2 ../resources/asm/HG008N.hap2.chr20.fa \
        --no-absolutize --no-check-exists \
        --output "config/samples_${PATTERN}.tsv" \
        --force

    python3 ../../setup_workflow.py \
        --samplesheet "config/samples_${PATTERN}.tsv" \
        --chm13-fasta ../../resource/reference/chm13v2.0_maskedY_rCRS.fa \
        --sex female \
        --mutation-caller "${MUTATION_CALLER}" \
        "${ANNOTATION_ARGS[@]}" \
        --hap1-satellite         resources/annotation/HG008N.hap1_dna-brnn.bed.gz \
        --hap2-satellite         resources/annotation/HG008N.hap2_dna-brnn.bed.gz \
        --simple-repeat          resources/annotation/HG008N.simple_repeats.bed.gz \
        --gtf-file               resources/annotation/HG008N.liftoff.gtf.gz \
        --gff-file               resources/annotation/HG008N.liftoff.gff.gz \
        --line1-bed              resources/annotation/HG008N.LINE1.bed.gz \
        --repeat-masker-bed      resources/annotation/HG008N.rmsk.bed.gz \
        --segdup-bed             resources/annotation/HG008N.segdup.bed.gz \
        --censat-bed             resources/annotation/HG008N.censat.bed.gz \
        --misassembly-hap1-bed   resources/annotation/HG008N.misassembly_hap1.bed.gz \
        --misassembly-hap2-bed   resources/annotation/HG008N.misassembly_hap2.bed.gz \
        --chain-to-grch38        resources/annotation/HG008N.to_grch38.chain \
        --chain-to-chm13         resources/annotation/HG008N.to_chm13.chain \
        --cancer-gene-census-tsv resources/annotation/cancer_gene_census.tsv \
        --cmrg-gene-list         resources/annotation/cmrg_genes.list \
        --gencode-transcript-bed resources/annotation/gencode.v46.basic.annotation.chr20.mane.transcript.bed.gz \
        --gnomad-bed             resources/annotation/gnomad.v4.1.sv.sites.bed.gz \
        --gnomad-vcf             resources/annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz \
        --grch38-fasta           ../../resource/reference/GRCh38.d1.vd1.fa \
        --copynumber-plot-sex-chrom false \
        --chm13-censat           ../../resource/reference/chm13v2.0_censat_v2.1.bed.gz \
        --pileup-no-baq true \
        --singularity-bind "${SING_BIND}" \
        --output-dir "./${PATTERN}" \
        --output "config/config_${PATTERN}.yaml" \
        --runner "run_${PATTERN}.sh" \
        --workflow-dir ../../workflow \
        --profile ../../profile/slurm \
        --images-dir ../../images \
        --jobs 24 \
        --nanomonsv-get-mem-mb 16000 \
        --deepsomatic-threads 4 \
        --deepsomatic-mem-mb 32000 \
        --clairs-threads 4 \
        --clairs-mem-mb 20000 \
        --bam-refiner-mem-mb 20000 \
        --bam-refiner-kmer-mem-mb 20000 \
        --methylation-mem-mb 24000 \
        --copynumber-mem-mb 24000 \
        --nanomonsv-parse-mem-mb 12000 \
        --nanomonsv-insert-classify-mem-mb 12000 \
        --nanomonsv-connect-mem-mb 8000 \
        --nanomonsv-merge-mem-mb 8000 \
        --nanomonsv-postprocess-mem-mb 8000 \
        --clairs-postprocess-realign-mem-mb 10000 \
        --deepsomatic-postprocess-realign-mem-mb 10000 \
        --clairs-postprocess-pileup-mem-mb 16000 \
        --deepsomatic-postprocess-pileup-mem-mb 16000 \
        --clairs-postprocess-haplotype-mem-mb 8000 \
        --deepsomatic-postprocess-haplotype-mem-mb 8000 \
        --assembly-bwa-index-mem-mb 16000 \
        --bam-refiner-kmer-threads 4 \
        --bam-refiner-threads 4 \
        --methylation-threads 4 \
        --copynumber-threads 4 \
        --nanomonsv-parse-threads 4 \
        --nanomonsv-get-threads 4 \
        --nanomonsv-insert-classify-threads 4 \
        --clairs-postprocess-realign-threads 4 \
        --deepsomatic-postprocess-realign-threads 4 \
        --force
done

echo
echo "configured: run_hifi.sh / run_ont.sh / run_both.sh (mutation_caller=${MUTATION_CALLER}, run_annotation=${RUN_ANNOTATION})"
