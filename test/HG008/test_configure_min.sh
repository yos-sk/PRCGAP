#!/bin/bash
# Configure the HG008 chr20 case with the MINIMAL input set: reads, assembly and
# the references only. Every per-assembly annotation PRCGAP can build itself is
# built (dna-brnn, liftoff, chain files, LINE-1, tandem repeats) and no
# annotation path key is passed. The optional annotation databases
# (repeat-masker / segdup / misassembly / cenSat / COSMIC / gnomAD / gencode)
# are left out too -- they are _opt_path() in the rules, so the annotate steps
# simply skip those columns.
#
# The point of the run is to prove the pipeline needs nothing but data +
# references, and to record per-rule time and memory. test_configure_full.sh remains
# the fixed-annotation case.
#
#   bash test_configure_min.sh   # writes config/config_min.yaml, run_min.sh
#   bash run_min.sh
#
# Only the HiFi pattern is configured: the memory/time record is per rule, and
# adding ONT would duplicate every rule rather than cover a new one.

set -euo pipefail

MUTATION_CALLER=${MUTATION_CALLER:-deepsomatic}
PATTERN=min

SING_BIND="$HOME"
HOME_REAL=$(readlink -f "$HOME")
if [ "${HOME_REAL}" != "$HOME" ]; then
    SING_BIND="$HOME,${HOME_REAL}"
fi

mkdir -p "./${PATTERN}"

# References only. No --hap1-satellite / --gtf-file / --line1-bed /
# --simple-repeat / --chain-to-* : the run_* switches generate all of those.
# The pair options write --samplesheet rather than reading it.
python3 ../../setup_workflow.py \
    --samplesheet "config/samples_${PATTERN}.tsv" \
    --tumor  HG008T \
    --normal HG008N \
    --tumor-hifi  resources/reads/hifi/HG008T.chr20.hifi.bam \
    --normal-hifi resources/reads/hifi/HG008N.chr20.hifi.bam \
    --assembly-hap1 resources/asm/HG008N.hap1.chr20.fa \
    --assembly-hap2 resources/asm/HG008N.hap2.chr20.fa \
    --sex female \
    --mutation-caller "${MUTATION_CALLER}" \
    --run-dna-brnn --run-liftoff --run-chain-files \
    --run-line1 --run-simple-repeat \
    --chm13-fasta            resources/reference/chm13v2.0_maskedY_rCRS.fa \
    --grch38-fasta           resources/reference/GRCh38.d1.vd1.fa \
    --grch38-gtf             resources/reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf \
    --grch38-centromeres     resources/reference/centromeres.txt.gz \
    --grch38-exclusions      resources/reference/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed \
    --chm13-censat           resources/reference/chm13v2.0_censat_v2.1.bed.gz \
    --copynumber-plot-sex-chrom false \
    --singularity-bind "${SING_BIND}" \
    --output-dir "./${PATTERN}" \
    --output "config/config_${PATTERN}.yaml" \
    --runner "run_${PATTERN}.sh" \
    --workflow-dir ../../workflow \
    --profile ../../profile/slurm \
    --images-dir ../../images \
    --jobs 24 \
    --force

echo
echo "configured: run_${PATTERN}.sh (minimal inputs, mutation_caller=${MUTATION_CALLER})"
