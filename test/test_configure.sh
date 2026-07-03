#!/bin/bash
# Regenerate test/config/test_config.yaml and test/run_test.sh for the
# HG008 chr20 test case. Run from the test/ directory:
#
#   bash test_configure.sh
#   bash run_test.sh
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
#   - bash resources/scripts/download_reference.sh
#       resources/reference/chm13v2.0_maskedY_rCRS.fa{,.fai}
#       resources/reference/GRCh38.d1.vd1.fa{,.fai}
#
# Manual prerequisites (no helper script):
#   - resources/annotation/cancer_gene_census.tsv  (requires COSMIC login; download manually)
#   - ../images/*.sif singularity images pre-staged

set -euo pipefail

# Generate the sample sheet (one tumor-normal pair sharing one assembly).
# --no-check-exists: configure does not require the multi-GB read/asm files to
# be staged yet. --no-absolutize: keep the resource paths relative to test/.
python3 ../set_sample_sheet.py \
    --tumor  HG008T \
    --tumor-ont  ../resources/reads/ont/HG008T.chr20.ont.bam \
    --tumor-hifi ../resources/reads/hifi/HG008T.chr20.hifi.bam \
    --normal HG008N \
    --normal-ont  ../resources/reads/ont/HG008N.chr20.ont.bam \
    --normal-hifi ../resources/reads/hifi/HG008N.chr20.hifi.bam \
    --assembly-hap1 ../resources/asm/HG008N.hap1.chr20.fa \
    --assembly-hap2 ../resources/asm/HG008N.hap2.chr20.fa \
    --no-absolutize --no-check-exists \
    --output config/test_samples.tsv \
    --force

python3 ../setup_workflow.py \
    --samplesheet config/test_samples.tsv \
    --chm13-fasta resources/reference/chm13v2.0_maskedY_rCRS.fa \
    --sex female \
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
    --grch38-fasta           resources/reference/GRCh38.d1.vd1.fa \
    --copynumber-plot-sex-chrom false \
    --chm13-censat           resources/reference/chm13v2.0_censat_v2.1.bed.gz \
    --pileup-no-baq true \
    --singularity-bind $HOME,/hshare1/ZETTAI_path_WA_slash_home_KARA/home/yosakam \
    --output-dir ./HG008 \
    --output config/test_config.yaml \
    --runner run_test.sh \
    --workflow-dir ../workflow \
    --images-dir ../images \
    --profile ../profile/sge \
    --jobs 8 \
    --nanomonsv-get-mem-mb 16000 \
    --deepsomatic-threads 4 \
    --deepsomatic-mem-mb 32000 \
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
    --clairs-threads 4 \
    --clairs-postprocess-realign-threads 4 \
    --deepsomatic-postprocess-realign-threads 4 \
    --force