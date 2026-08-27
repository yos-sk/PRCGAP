#!/bin/bash
# Configure the HG008 chr20 test case: one tumor/normal pair with both HiFi and
# ONT reads, against the pre-extracted chr20 annotation under resources/.
#
#   bash test_configure.sh    # writes config/ and run_hg008.sh
#   sbatch run_hg008.sh       # or: bash run_hg008.sh, depending on --profile
#
# This is the reference test case: it exercises every stage of the pipeline on a
# single chromosome, so a full run finishes in a couple of hours rather than
# days. Resource requests below are sized for chr20 and are deliberately far
# smaller than the whole-genome defaults in setup_workflow.py.
#
# The switchable annotation steps (dna-brnn, chain files, LINE-1, tandem
# repeats) are off here: the fixture supplies them, which keeps the run focused
# on the calling pipeline. liftoff has no switch -- the gene annotation is
# always built by the workflow -- so the GRCh38 FASTA and GTF are required.
#
# Set MUTATION_CALLER to exercise the other point-mutation caller:
#   MUTATION_CALLER=clairs bash test_configure.sh
#   MUTATION_CALLER=both   bash test_configure.sh
#
# Point-mutation calling runs on HiFi only (mutation_seqtypes defaults to
# `primary`). Pass --mutation-seqtypes all to call on ONT as well.
#
# Prerequisites. References are downloaded once for the whole repository and cut
# down here; everything else is specific to this test case.
#
#   1. bash ../../resource/scripts/download_reference.sh
#        ../../resource/reference/{chm13v2.0_maskedY_rCRS,GRCh38.d1.vd1}.fa{,.fai}
#        ../../resource/reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf
#        ../../resource/reference/{centromeres.txt.gz,chm13v2.0_censat_v2.1.bed.gz}
#        ../../resource/reference/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed
#        ../../resource/reference/MANE.GRCh38.v1.3.summary.txt.gz   (read directly)
#   2. bash resources/scripts/extract_chr20_reference.sh
#        resources/reference/*_chr20.*   -- the chr20 cut of the above; the _chr20
#        suffix keeps a cut from being mistaken for the genome-wide original
#   3. bash resources/scripts/download_data.sh
#        resources/reads/{hifi,ont}/HG008{T,N}.chr20.{hifi,ont}.bam
#        resources/asm/HG008N.hap{1,2}.chr20.fa
#   4. bash resources/scripts/download_annotation.sh
#        resources/annotation/cmrg_genes.list
#        resources/annotation/gnomad.v4.1.sv.sites.bed.gz{,.tbi}
#        resources/annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz{,.tbi}
#   5. ../../images/*.sif singularity images pre-staged (see images/pull_images.sh)
#
# One prerequisite has no helper script: resources/annotation/cancer_gene_census.tsv
# requires a COSMIC login. Drop --cancer-gene-census-tsv to run without it; the
# `cgc` column is then filled with the not-evaluated placeholder.

set -euo pipefail

MUTATION_CALLER=${MUTATION_CALLER:-deepsomatic}
PATTERN=HG008

# $HOME may be a symlink (e.g. /home/<user> -> /lustre/home/<user>). Config paths
# keep the user-facing form, but rules that resolve symlinks (assembly_bwa_index's
# `readlink -f`) produce the target form, so both have to be visible inside the
# container.
SING_BIND="$HOME"
HOME_REAL=$(readlink -f "$HOME")
if [ "${HOME_REAL}" != "$HOME" ]; then
    SING_BIND="$HOME,${HOME_REAL}"
fi

mkdir -p "./${PATTERN}"

# Read and annotation paths are resolved relative to the run's output directory
# (snakemake runs with cwd=output_dir), which is ./${PATTERN} here. The pair
# options write the sample sheet rather than reading one.
python3 ../../setup_workflow.py \
    --samplesheet "config/samples_${PATTERN}.tsv" \
    --tumor  HG008T \
    --normal HG008N \
    --tumor-hifi  resources/reads/hifi/HG008T.chr20.hifi.bam \
    --normal-hifi resources/reads/hifi/HG008N.chr20.hifi.bam \
    --tumor-ont   resources/reads/ont/HG008T.chr20.ont.bam \
    --normal-ont  resources/reads/ont/HG008N.chr20.ont.bam \
    --assembly-hap1 resources/asm/HG008N.hap1.chr20.fa \
    --assembly-hap2 resources/asm/HG008N.hap2.chr20.fa \
    --sex female \
    --mutation-caller "${MUTATION_CALLER}" \
    --no-run-dna-brnn --no-run-chain-files \
    --no-run-line1 --no-run-simple-repeat \
    --chm13-fasta            resources/reference/chm13v2.0_maskedY_rCRS_chr20.fa \
    --grch38-fasta           resources/reference/GRCh38.d1.vd1_chr20.fa \
    --grch38-gtf             resources/reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format_chr20.gtf \
    --chm13-censat           resources/reference/chm13v2.0_censat_v2.1_chr20.bed.gz \
    --mane-summary           ../../resource/reference/MANE.GRCh38.v1.3.summary.txt.gz \
    --hap1-satellite         resources/annotation/HG008N.hap1_dna-brnn.bed.gz \
    --hap2-satellite         resources/annotation/HG008N.hap2_dna-brnn.bed.gz \
    --simple-repeat          resources/annotation/HG008N.simple_repeats.bed.gz \
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
    --gnomad-bed             resources/annotation/gnomad.v4.1.sv.sites.bed.gz \
    --gnomad-vcf             resources/annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz \
    --copynumber-plot-sex-chrom false \
    --pileup-no-baq true \
    --singularity-bind "${SING_BIND}" \
    --output-dir "./${PATTERN}" \
    --output "config/config_${PATTERN}.yaml" \
    --runner "run_${PATTERN}.sh" \
    --workflow-dir ../../workflow \
    --profile ../../profile/slurm \
    --images-dir ../../images \
    --jobs 24 \
    --liftoff-threads 4 --liftoff-mem-mb 16000 \
    --bam-refiner-threads 4 --bam-refiner-mem-mb 20000 \
    --bam-refiner-kmer-threads 4 --bam-refiner-kmer-mem-mb 20000 \
    --assembly-bwa-index-mem-mb 16000 \
    --methylation-threads 4 --methylation-mem-mb 24000 \
    --copynumber-threads 4 --copynumber-mem-mb 24000 \
    --deepsomatic-threads 4 --deepsomatic-mem-mb 32000 \
    --clairs-threads 4 --clairs-mem-mb 20000 \
    --nanomonsv-parse-threads 4 --nanomonsv-parse-mem-mb 12000 \
    --nanomonsv-get-threads 4 --nanomonsv-get-mem-mb 16000 \
    --nanomonsv-insert-classify-threads 4 --nanomonsv-insert-classify-mem-mb 12000 \
    --nanomonsv-postprocess-mem-mb 8000 \
    --nanomonsv-connect-mem-mb 8000 \
    --nanomonsv-merge-mem-mb 8000 \
    --clairs-postprocess-realign-threads 4 --clairs-postprocess-realign-mem-mb 10000 \
    --deepsomatic-postprocess-realign-threads 4 --deepsomatic-postprocess-realign-mem-mb 10000 \
    --clairs-postprocess-pileup-mem-mb 16000 \
    --deepsomatic-postprocess-pileup-mem-mb 16000 \
    --clairs-postprocess-haplotype-mem-mb 8000 \
    --deepsomatic-postprocess-haplotype-mem-mb 8000 \
    --force

echo
echo "configured: run_${PATTERN}.sh (HiFi + ONT, mutation_caller=${MUTATION_CALLER})"
