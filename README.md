# PRCGAP - Personalized Reference genome-based Cancer Genome Analysis Pipeline

Snakemake-based workflow for comprehensive analysis of complete cancer genomes using long-read sequencing data (PacBio HiFi and Oxford Nanopore).

> 📖 [**Full documentation:**](https://yos-sk.github.io/prcgapdoc)
>
> This README is a short overview. Setup, configuration, sample sheet format, cluster execution, output layout, and container/image details are all covered in prcgapdoc.

## Analysis Steps

| Step | Description | Container |
|------|-------------|-----------|
| BAM Refiner | Align reads to phased de novo assemblies | bam_refiner |
| Methylation | Call methylation (HiFi: pb-CpG-tools, ONT: modkit) | methylation |
| Copy Number | Analyze copy number variations | copynumber |
| Nanomonsv Parse | Parse structural variants from BAM | nanomonsv |
| Nanomonsv Get | Extract SV calls (tumor-normal paired) | nanomonsv |
| Nanomonsv Postprocess | Post-process SV results | nanomonsv_postprocess |
| Nanomonsv Insert Classify | Classify insertion types | nanomonsv |
| Nanomonsv Connect | Connect fragmented SVs | nanomonsv_postprocess |
| Nanomonsv Merge | Merge HiFi and ONT SV results | nanomonsv_postprocess |
| ClairS | Call point mutations with ClairS | clairs |
| DeepSomatic | Call somatic variants with DeepSomatic | deepsomatic |
| ClairS Postprocess | Realignment, pileup, haplotyping for ClairS results | point_mutation_postprocess |
| DeepSomatic Postprocess | Realignment, pileup, haplotyping for DeepSomatic results | point_mutation_postprocess |

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) (currently tested with 7.32.4). Snakemake 8.x changes several CLI flags (e.g. `--sdm` / `--software-deployment-method` replaces `--use-singularity`, profile format differs), so the commands shown in prcgapdoc assume 7.x.
- [Apptainer/Singularity](https://apptainer.org/)
- Singularity images for each analysis tool (see prcgapdoc for the registry URLs / build instructions)

## Quick Start

See prcgapdoc for the full quick-start guide. The general flow is:

1. Prepare a sample sheet (TSV with sample, type, ONT/HiFi data paths, and per-haplotype assemblies).
2. Generate a `config.yaml` with `setup_workflow.py`.
3. Run with `snakemake --sdm apptainer ...`, optionally using the `profile/uge` or `profile/slurm` cluster profile.
