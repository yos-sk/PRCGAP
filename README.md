# PRCGAP - Personalized Reference genome-based Cancer Genome Analysis Pipeline

Snakemake-based workflow for comprehensive analysis of cancer genomes using long-read sequencing data (PacBio HiFi and Oxford Nanopore) based on personalized reference genome.

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
| SV Prep | Filter nanomonsv result, extract breakpoint BED | annotation |
| SV Coordconv | Lift breakpoints to GRCh38 / CHM13 (optional) | annotation |
| SV Annotation | Gene / RepeatMasker / centromere / segdup / kmer / liftover / gnomAD / misassembly | annotation |
| SV Reclassify | Reclassify inter-contig SVs using CHM13-normalized directions | annotation |
| Mut Prep (SNV/INDEL) | Split haplotyped.bed → SNV/INDEL TSV + coordconv BED | annotation |
| Mut Coordconv (SNV/INDEL) | Lift mutation positions to GRCh38 / CHM13 (optional) | annotation |
| SNV / INDEL Annotation | Lifted coords / Gene / RepeatMasker / centromere / segdup / misassembly / cross-tool check | annotation |

### Annotation resources (optional)

Annotation steps each consume an external resource file. **All keys are optional**; if the path is left empty (the default), the corresponding annotation column is skipped — Snakemake's DAG still resolves, the chain just emits fewer columns. Provide them under `setup_workflow.py` flags or directly in `config.yaml`:

| config key | CLI flag | Used by |
|---|---|---|
| `gtf_file` | `--gtf-file` | nanomonsv insert_classify (liftoff GTF) |
| `gff_file` | `--gff-file` | SV / SNV / INDEL gene annotation (tabix-indexed liftoff GFF) |
| `chain_to_grch38` | `--chain-to-grch38` | SV / SNV / INDEL liftover to GRCh38 |
| `chain_to_chm13` | `--chain-to-chm13` | SV / SNV / INDEL liftover to CHM13 |
| `repeat_masker_bed` | `--repeat-masker-bed` | SV / SNV / INDEL |
| `segdup_bed` | `--segdup-bed` | SV / SNV / INDEL |
| `censat_bed` | `--censat-bed` | SV / SNV / INDEL |
| `misassembly_hap{1,2}_bed` | `--misassembly-hap{1,2}-bed` | SV / SNV / INDEL (optional) |
| `cancer_gene_census_tsv` | `--cancer-gene-census-tsv` | SV / SNV / INDEL gene annotation |
| `cmrg_gene_tsv` | `--cmrg-gene-tsv` | SNV / INDEL gene annotation |
| `gencode_transcript_bed` | `--gencode-transcript-bed` | SNV / INDEL gene annotation |
| `gnomad_bed` | `--gnomad-bed` | SV gnomAD annotation (requires `chain_to_grch38`) |
| `gnomad_vcf` | `--gnomad-vcf` | SNV / INDEL gnomAD annotation (requires `chain_to_grch38`) |
| `grch38_fasta` | `--grch38-fasta` | INDEL transanno liftover `--query` (requires `chain_to_grch38`) |

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) (currently tested with 7.32.4). Snakemake 8.x changes several CLI flags (e.g. `--sdm` / `--software-deployment-method` replaces `--use-singularity`, profile format differs), so the commands shown in prcgapdoc assume 7.x.
- [Apptainer/Singularity](https://apptainer.org/)
- Singularity images for each analysis tool (see prcgapdoc for the registry URLs / build instructions)

## Quick Start

See prcgapdoc for the full quick-start guide. The general flow is:

1. Prepare a sample sheet (TSV with sample, type, ONT/HiFi data paths, and per-haplotype assemblies).
2. Generate a `config.yaml` with `setup_workflow.py`.
3. Run with `snakemake --sdm apptainer ...`, optionally using the `profile/uge` or `profile/slurm` cluster profile.
