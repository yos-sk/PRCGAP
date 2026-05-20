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
    --singularity-bind /lustre1/home/yosakam2,/home/yosakam2 \
    --bam-refiner-kmer-threads          8 --bam-refiner-kmer-mem-mb           16000 \
    --bam-refiner-threads               8 --bam-refiner-mem-mb                16000 \
    --assembly-bwa-index-threads        1 --assembly-bwa-index-mem-mb          8000 \
    --methylation-threads               8 --methylation-mem-mb                16000 \
    --copynumber-threads                8 --copynumber-mem-mb                 16000 \
    --nanomonsv-parse-threads           8 --nanomonsv-parse-mem-mb            16000 \
    --nanomonsv-get-threads             8 --nanomonsv-get-mem-mb              32000 \
    --nanomonsv-postprocess-threads     1 --nanomonsv-postprocess-mem-mb       8000 \
    --nanomonsv-insert-classify-threads 8 --nanomonsv-insert-classify-mem-mb  16000 \
    --nanomonsv-connect-threads         1 --nanomonsv-connect-mem-mb           8000 \
    --nanomonsv-merge-threads           1 --nanomonsv-merge-mem-mb             8000 \
    --clairs-threads                    8 --clairs-mem-mb                     16000 \
    --deepsomatic-threads               8 --deepsomatic-mem-mb                16000 \
    --clairs-postprocess-threads        1 --clairs-postprocess-mem-mb          8000 \
    --clairs-postprocess-split-threads  1 --clairs-postprocess-split-mem-mb    8000 \
    --clairs-postprocess-realign-threads        8 --clairs-postprocess-realign-mem-mb        16000 \
    --clairs-postprocess-pileup-threads         4 --clairs-postprocess-pileup-mem-mb         16000 \
    --clairs-postprocess-haplotype-threads      4 --clairs-postprocess-haplotype-mem-mb      16000 \
    --deepsomatic-postprocess-threads          1 --deepsomatic-postprocess-mem-mb            8000 \
    --deepsomatic-postprocess-split-threads    1 --deepsomatic-postprocess-split-mem-mb      8000 \
    --deepsomatic-postprocess-realign-threads  8 --deepsomatic-postprocess-realign-mem-mb   16000 \
    --deepsomatic-postprocess-pileup-threads   4 --deepsomatic-postprocess-pileup-mem-mb    16000 \
    --deepsomatic-postprocess-haplotype-threads 4 --deepsomatic-postprocess-haplotype-mem-mb 16000 \
    --prep-sv-threads        1 --prep-sv-mem-mb         4000 \
    --coordconv-sv-threads   1 --coordconv-sv-mem-mb    4000 \
    --annotate-sv-threads    1 --annotate-sv-mem-mb     8000 \
    --reclassify-sv-threads  1 --reclassify-sv-mem-mb   8000 \
    --prep-mut-threads       1 --prep-mut-mem-mb        4000 \
    --coordconv-mut-threads  1 --coordconv-mut-mem-mb   4000 \
    --annotate-mut-threads   1 --annotate-mut-mem-mb    8000 \
    --output-dir ./HG008 \
    --output config/test_config.yaml \
    --runner run_test.sh \
    --workflow-dir ../workflow \
    --images-dir ../images \
    --jobs 8 \
    --force
