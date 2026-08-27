# test/resources/scripts

Helpers for materializing the chr20 / HG008 test fixtures consumed by the PRCGAP test workflow. All paths below are relative to this directory (`test/HG008/resources/scripts/`) — run the scripts from here.

The reference genomes are not downloaded here. `resource/scripts/download_reference.sh` at the repository root fetches them once into `resource/reference/`, and every run in the repository reads them from there. `extract_chr20_reference.sh` below cuts that one copy down to chr20 for this test case; the scripts here download only what is specific to HG008.

## Dependency

Before running any of the scripts below, make sure the following tools are on `PATH`:

- `python3`
- [samtools](https://github.com/samtools/samtools) (used by `download_data.sh` for `view` / `faidx`)
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

## extract_chr20_reference.sh

Cuts the genome-wide references in `resource/reference/` down to chr20 and writes them under `../reference/` with a `_chr20` suffix, so a cut is never mistaken for the genome-wide original it came from. The assembly is chr20 only, so aligning it against 3 Gb of reference would spend most of its time building a minimap2 index of sequence that cannot match.

| Path | Cut from |
|---|---|
| `../reference/chm13v2.0_maskedY_rCRS_chr20.fa{,.fai}` | `resource/reference/chm13v2.0_maskedY_rCRS.fa` |
| `../reference/GRCh38.d1.vd1_chr20.fa{,.fai}` | `resource/reference/GRCh38.d1.vd1.fa` |
| `../reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format_chr20.gtf` | the matching genome-wide GTF |
| `../reference/centromeres_chr20.txt.gz` | UCSC hg38 centromeres table (chrom is field 2) |
| `../reference/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2_chr20.bed` | GRC exclusion regions |
| `../reference/chm13v2.0_censat_v2.1_chr20.bed.gz{,.tbi}` | CHM13 cenSat BED |

The GTF has to be cut as well: against a chr20-only FASTA, liftoff would still try to lift every other chromosome's genes and find no sequence for them.

`MANE.GRCh38.v1.3.summary.txt.gz` is read straight from `resource/reference/` and is not cut — the gene annotation matches on transcript ID, not coordinates.

## download_annotation.sh

Downloads the annotation resources referenced by the PRCGAP annotation steps and prepares the chr20 / MANE subsets used in tests. All outputs land under `../annotation/`.

| Path | Source / Step |
|---|---|
| `../annotation/cmrg_genes.list` | Challenging Medically Relevant Gene list (Wagner et al., Nat Biotechnol 2022 supplementary TSV) filtered by `extract_cmrg_gene.py` (autosomal, `Category == TRUE`) |
| `../annotation/gnomad.v4.1.sv.sites.bed.gz{,.tbi}` | gnomAD v4.1 SV sites (whole genome), re-sorted + bgzipped + tabix-indexed |
| `../annotation/gnomad.genomes.v4.1.sites.chr20.vcf.bgz{,.tbi}` | gnomAD v4.1 genome SNV/indel sites, chr20 |

The CMRG supplementary TSV is removed at the end. The Python helper `extract_cmrg_gene.py` lives in this directory.

The gene annotation now selects transcripts from the MANE summary directly, so the GENCODE-derived transcript BED this script used to build is gone.

Manual step:
- **Cancer Gene Census** must be downloaded after registering / logging in at [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v104/cancergenecensus) and placed at `../annotation/cancer_gene_census.tsv`; the script intentionally skips it.