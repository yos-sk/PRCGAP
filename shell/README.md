# PRCGAP shell implementation

This directory holds the shell-script implementation of PRCGAP (SLURM `sbatch`) that was
used for the analyses in the PRCGAP paper. It runs with the container images of this
release, the same images as the Snakemake workflow (`../workflow`), which implements the
same steps.

Annotation and reclassification of the SV / point-mutation calls are not part of this
implementation (in the Snakemake workflow they are the `annotate_*` rules).

## Steps

| # | Step | Scripts | Image |
|---|---|---|---|
| 1 | Read alignment to the personal diploid assembly and refinement | `scripts/bam_refiner` | `bam_refiner` |
| 2 | Methylation (HiFi: pb-CpG-tools, ONT: modkit) | `scripts/methylation` | `methylation` |
| 3 | Copy number (reference tables, depth, CBS, plot) | `scripts/copynumber` | `copynumber` |
| 4–7 | SV: nanomonsv parse / get, postprocess, insert classification | `scripts/nanomonsv` | `nanomonsv`, `nanomonsv_postprocess` |
| 8–9 | nanomonsv connect / HiFi+ONT merge (optional, off by default) | `scripts/nanomonsv` | `nanomonsv`, `nanomonsv_postprocess` |
| 10–11 | Point mutations: ClairS (optional, off by default), DeepSomatic | `scripts/clairs`, `scripts/deepsomatic` | `clairs`, `deepsomatic` |
| 12–13 | Point-mutation postprocess (realignment, pileup, haplotyping) | `scripts/mutation_postprocess` | `point_mutation_postprocess` |

## Container images

The scripts read the images from `./images`. Pull them into `../images` and link that
directory here:

```bash
(cd ../images && bash pull_images.sh)
ln -s ../images images
```

| Image | Docker tag |
|---|---|
| `bam_refiner.sif` | `yosakam2/bam_refiner:v0.3.6` |
| `methylation.sif` | `yosakam2/methylation:v0.1.0` |
| `copynumber.sif` | `yosakam2/copynumber:v0.3.0` |
| `nanomonsv.sif` | `yosakam2/nanomonsv:v0.8.0` |
| `nanomonsv_postprocess.sif` | `yosakam2/nanomonsv_postprocess:v0.2.5` |
| `clairs.sif` | `yosakam2/clairs:v0.4.0` |
| `deepsomatic.sif` | `yosakam2/deepsomatic:v1.8.0` |
| `point_mutation_postprocess.sif` | `yosakam2/mutation_postprocess:v0.1.2` |

## Usage

1. Copy `conf.sh` and fill in the sample, assembly, sequencing-data and reference paths.
   Steps are switched on and off with the `TRUE` / `FALSE` flags at the end.
2. If the data live outside `$HOME` and the current directory, list their directories in
   `BIND_PATHS` (comma-separated) so the containers can see them.
3. Run from this directory:

```bash
bash run.sh my_conf.sh
```

Results are written under `OUTPUT_PREFIX` (`bam_refiner/`, `methylation/`, `copynumber/`,
`nanomonsv/`, `deepsomatic/`, `deepsomatic_post/`, ...). Each step waits for the jobs of
the step it depends on (`sbatch --dependency`) when that step is switched on in the same run.

The `#SBATCH` lines (partition `-p`, cores, memory) are those of the original cluster;
adjust them to yours.

### Copy-number parameters

| Variable | Meaning |
|---|---|
| `CHM13_CENSAT` | CHM13 cenSat v2.1 BED(.gz) |
| `CENSAT_BED` | cenSat annotation of the assembly for the plot (empty: `HAP1_SATELLITE` + `HAP2_SATELLITE`) |
| `CN_BIN_WIDTH` | bin width of the depth-ratio histogram for the ploidy estimate (default 0.05) |
| `CN_PLOIDY_HAP1`, `CN_PLOIDY_HAP2` | `-` to estimate the tumor ploidy; set both to override |
| `CN_HAP1_LABEL`, `CN_HAP2_LABEL` | haplotype labels in the plot |
