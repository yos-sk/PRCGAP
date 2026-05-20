# test/resources/scripts

Helpers for materializing the chr20 / HG008 test fixtures consumed by the PRCGAP test workflow. All paths below are relative to this directory (`test/resources/scripts/`) — run the scripts from here.

## Dependency

Before running any of the scripts below, make sure the following tools are on `PATH`:

- `python3`
- [samtools](https://github.com/samtools/samtools) (used by `download_data.sh` and `download_reference.sh` for `view` / `faidx`)
- [bgzip and tabix](https://github.com/samtools/htslib) (used by `download_annotation.sh` and `extract_haplotypes.sh` for compression and indexing)

A simple `which samtools bgzip tabix wget python3` should print a path for each before you start.

## download_data.sh

Downloads HG008 tumor / normal long reads and the matched HG008N de novo assembly, slicing each to the chr20-equivalent contigs used by the test workflow.

Outputs (written one level up under `test/resources/`):

| Path | Source |
|---|---|
| `../reads/ont/HG008{T,N}.chr20.ont.bam` | GIAB Liss_lab ONT R10.4.1 BAMs (CHM13v2.0), chr20 slice via `samtools view -F 2308` |
| `../reads/hifi/HG008{T,N}.chr20.hifi.bam` | GIAB Liss_lab PacBio HiFi Revio BAMs (CHM13v2.0), chr20 slice via `samtools view -F 2308` |
| `../asm/HG008N.hap{1,2}.fa{,.fai}` | Verkko 2.2 HG008N haplotype assemblies (full) |
| `../asm/HG008N.hap{1,2}.chr20.fa` | `samtools faidx` extract of `haplotype1-0000020` / `haplotype2-0000110` (the chr20 contigs) |

Notes:
- Filter flag `-F 2308` drops unmapped, secondary, and supplementary alignments.
- The full-assembly FASTAs are kept alongside the chr20 extracts because `samtools faidx` needs the index built against the full assembly.

## download_reference.sh

Downloads the CHM13 and GRCh38 reference assemblies into `../reference/`.

| Path | Source |
|---|---|
| `../reference/chm13v2.0_maskedY_rCRS.fa{,.fai}` | T2T-CHM13 v2.0 (masked Y + rCRS), from the Human Pangenomics S3 bucket |
| `../reference/GRCh38.d1.vd1.fa{,.fai}` | GRCh38.d1.vd1 from GDC (used as transanno `--query` for INDEL liftover) |

## download_annotation.sh

Downloads the annotation resources referenced by the PRCGAP annotation steps and prepares the chr20 / MANE subsets used in tests. All outputs land under `../annotation/`.

| Path | Source / Step |
|---|---|
| `../annotation/cmrg_genes.list` | Challenging Medically Relevant Gene list (Wagner et al., Nat Biotechnol 2022 supplementary TSV) filtered by `extract_cmrg_gene.py` (autosomal, `Category == TRUE`) |
| `../annotation/gnomad.v4.1.sv.sites.bed.gz{,.tbi}` | gnomAD v4.1 SV sites (whole genome), re-sorted + bgzipped + tabix-indexed |
| `../annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz{,.tbi}` | gnomAD v4.1 genome SNV/indel sites, chr20 |
| `../annotation/gencode.v46.basic.annotation.chr20.mane.transcript.bed.gz` | chr20 transcript records whose transcript_id appears in MANE v1.3, produced by `proc_gencode_bed_mane_chr20.py` then `gzip`-ed |

The intermediate `gencode.v46.basic.annotation.gff3.gz`, `MANE.GRCh38.v1.3.summary.txt.gz`, and the CMRG supplementary TSV are removed at the end.

Python helpers `extract_cmrg_gene.py` and `proc_gencode_bed_mane_chr20.py` live in this directory.

Manual step:
- **Cancer Gene Census** must be downloaded after registering / logging in at [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v104/cancergenecensus) and placed at `../annotation/cancer_gene_census.tsv`; the script intentionally skips it.

## Companion scripts (run separately)

- `extract_haplotypes.sh` — filters externally-prepared full-assembly annotation files (dna-brnn, RepeatMasker, Liftoff, CenSat, SegDup, misassembly, liftover chains) down to the two test contigs and writes them under `../annotation/`. Edit the `*_SRC` paths at the top before running.
