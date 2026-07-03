#!/usr/bin/env python3
"""
Setup script for PRCGAP Snakemake workflow.
Generates config.yaml and a runner shell script from command-line arguments.
"""

import argparse
import os
import stat
import sys
import yaml
from pathlib import Path


def _image_path(explicit, images_dir, tool):
    """Return user-supplied image path or fall back to <images_dir>/<tool>.sif."""
    if explicit:
        return explicit
    return os.path.join(images_dir, f"{tool}.sif")


# (module_name, default_threads, default_mem_mb) for each rule whose
# resources can be overridden from the CLI.
_RESOURCE_DEFAULTS = [
    ("bam_refiner_kmer", 16, 128000),
    ("bam_refiner", 16, 128000),
    ("assembly_bwa_index", 1, 16000),
    ("methylation", 16, 128000),
    ("copynumber", 16, 128000),
    ("nanomonsv_parse", 16, 128000),
    ("nanomonsv_get", 16, 240000),
    ("nanomonsv_postprocess", 1, 30000),
    ("nanomonsv_insert_classify", 16, 128000),
    ("nanomonsv_connect", 1, 30000),
    ("nanomonsv_merge", 1, 30000),
    ("clairs", 16, 128000),
    ("deepsomatic", 16, 128000),
    ("clairs_postprocess", 1, 32000),
    ("clairs_postprocess_split", 1, 32000),
    ("clairs_postprocess_realign", 16, 128000),
    ("clairs_postprocess_pileup", 6, 240000),
    ("clairs_postprocess_haplotype", 4, 64000),
    ("deepsomatic_postprocess", 1, 32000),
    ("deepsomatic_postprocess_split", 1, 32000),
    ("deepsomatic_postprocess_realign", 16, 128000),
    ("deepsomatic_postprocess_pileup", 6, 240000),
    ("deepsomatic_postprocess_haplotype", 4, 64000),
    ("prep_sv", 1, 8000),
    ("coordconv_sv", 1, 8000),
    ("annotate_sv", 1, 16000),
    ("reclassify_sv", 1, 16000),
    ("gff_to_bed", 1, 8000),
    ("prep_mut", 1, 8000),
    ("coordconv_mut", 1, 8000),
    ("annotate_mut", 1, 16000),
    ("indel_hap_reference", 1, 8000),
    ("mut_vcf_index", 1, 8000),
    ("bed2vcf_mut", 1, 8000),
    ("liftvcf_mut", 1, 16000),
]


def _abs(path):
    """Absolutise a path WITHOUT following symlinks. Empty/None passes through.

    Snakemake runs with cwd=--directory (output_dir), so any relative path
    in config.yaml would resolve under output_dir and miss the user's input
    files. We absolutise at config-generation time.

    We use os.path.abspath rather than Path.resolve() so symlinks like
    /home/<user> -> /hshare1/.../home/<user> stay as the user-facing path.
    Following symlinks would yield paths that aren't bind-mounted into the
    singularity/apptainer container, breaking --pwd inside the container.
    """
    if not path:
        return path
    return os.path.abspath(os.path.expanduser(path))


# Profile config.yaml keys whose values are paths to scripts living inside
# the profile dir; snakemake resolves them relative to cwd (not the profile
# dir), so we rewrite them to absolute paths.
_PROFILE_PATH_KEYS = (
    "jobscript",
    "cluster",
    "cluster-status",
    "cluster-cancel",
    "cluster-generic-submit-cmd",
    "cluster-generic-status-cmd",
    "cluster-generic-cancel-cmd",
)


def _absolutise_profile_paths(profile_dir: Path):
    """Rewrite path-valued keys in <profile_dir>/config.yaml to absolute.

    Only rewrites values that resolve to an existing file inside the profile
    dir (so already-absolute or already-resolved paths are left alone).
    """
    cfg = profile_dir / "config.yaml"
    if not cfg.exists():
        return
    with open(cfg) as f:
        data = yaml.safe_load(f) or {}

    changed = False
    for key in _PROFILE_PATH_KEYS:
        val = data.get(key)
        if not isinstance(val, str) or not val:
            continue
        if Path(val).is_absolute():
            continue
        # Use abspath (no symlink follow) so the rewritten path matches what
        # gets bind-mounted into the container.
        candidate = Path(os.path.abspath(profile_dir / val))
        if candidate.exists():
            data[key] = str(candidate)
            changed = True

    if changed:
        with open(cfg, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def create_config(args):
    """Create configuration dictionary from arguments."""
    config = {
        "output_dir": _abs(args.output_dir),
        "samplesheet": _abs(args.samplesheet),
        "chm13_fasta": _abs(args.chm13_fasta),
        "singularity_images": {
            "bam_refiner": _abs(_image_path(args.bam_refiner_image, args.images_dir, "bam_refiner")),
            "methylation": _abs(_image_path(args.methylation_image, args.images_dir, "methylation")),
            "copynumber": _abs(_image_path(args.copynumber_image, args.images_dir, "copynumber")),
            "nanomonsv": _abs(_image_path(args.nanomonsv_image, args.images_dir, "nanomonsv")),
            "nanomonsv_postprocess": _abs(_image_path(args.nanomonsv_postprocess_image, args.images_dir, "nanomonsv_postprocess")),
            "clairs": _abs(_image_path(args.clairs_image, args.images_dir, "clairs")),
            "deepsomatic": _abs(_image_path(args.deepsomatic_image, args.images_dir, "deepsomatic")),
            "point_mutation_postprocess": _abs(_image_path(args.point_mutation_postprocess_image, args.images_dir, "point_mutation_postprocess")),
            "annotation": _abs(_image_path(args.annotation_image, args.images_dir, "annotation")),
        },
        "resources": {
            name: {
                "threads": getattr(args, f"{name}_threads"),
                "mem_mb": getattr(args, f"{name}_mem_mb"),
            }
            for name, _, _ in _RESOURCE_DEFAULTS
        },
        "hap1_satellite": _abs(args.hap1_satellite) or "",
        "hap2_satellite": _abs(args.hap2_satellite) or "",
        "sex": args.sex,
        "simple_repeat": _abs(args.simple_repeat) or "",
        "gtf_file": _abs(args.gtf_file) or "",
        "gff_file": _abs(args.gff_file) or "",
        "line1_bed": _abs(args.line1_bed) or "",
        # ---- copy number plot params ----
        "copynumber_binwidth": args.copynumber_binwidth,
        "copynumber_hap1_label": args.copynumber_hap1_label,
        "copynumber_hap2_label": args.copynumber_hap2_label,
        "copynumber_ploidy_hap1": args.copynumber_ploidy_hap1,
        "copynumber_ploidy_hap2": args.copynumber_ploidy_hap2,
        "copynumber_plot_sex_chrom": args.copynumber_plot_sex_chrom,
        "chm13_censat": _abs(args.chm13_censat) or "",
        # ---- pileup (mutation postprocess) params ----
        "pileup_no_baq": args.pileup_no_baq,
        # ---- annotation resources (optional) ----
        "chain_to_grch38": _abs(args.chain_to_grch38) or "",
        "chain_to_chm13": _abs(args.chain_to_chm13) or "",
        "repeat_masker_bed": _abs(args.repeat_masker_bed) or "",
        "segdup_bed": _abs(args.segdup_bed) or "",
        "censat_bed": _abs(args.censat_bed) or "",
        "misassembly_hap1_bed": _abs(args.misassembly_hap1_bed) or "",
        "misassembly_hap2_bed": _abs(args.misassembly_hap2_bed) or "",
        "cancer_gene_census_tsv": _abs(args.cancer_gene_census_tsv) or "",
        "cmrg_gene_tsv": _abs(args.cmrg_gene_list) or "",
        "gencode_transcript_bed": _abs(args.gencode_transcript_bed) or "",
        "gnomad_bed": _abs(args.gnomad_bed) or "",
        "gnomad_vcf": _abs(args.gnomad_vcf) or "",
        "grch38_fasta": _abs(args.grch38_fasta) or "",
    }
    return config


def write_runner(args, config_path: Path, runner_path: Path):
    """Emit a runner shell script invoking snakemake.

    Pattern (after cosigt/organize.py):
      - With --profile  -> snakemake --profile <profile> ...
      - Without         -> snakemake -j <threads> ...
      - Always uses --use-singularity (conda is not supported).
    """
    # Use _abs (os.path.abspath) so symlinks aren't followed; the resolved
    # path must match what's bind-mounted into the singularity container.
    workflow_dir = _abs(args.workflow_dir)
    snakefile = os.path.join(workflow_dir, "Snakefile")
    output_dir = _abs(args.output_dir)
    config_abs = _abs(str(config_path))

    cmd_lines = ["#!/bin/bash", "set -euo pipefail", ""]

    # Targets to build. When --targets is given, bake them in as defaults so the
    # runner builds only those (e.g. copynumber alone); targets passed on the
    # runner's command line still override them via "${@:-...}". Otherwise pass
    # "$@" straight through and let snakemake build the full pipeline by default.
    if args.targets:
        targets_array = " ".join(args.targets)
        cmd_lines.append(
            "# Default targets baked in by setup_workflow.py --targets;\n"
            "# arguments passed to this script override them.")
        cmd_lines.append(f"DEFAULT_TARGETS=({targets_array})")
        cmd_lines.append("")
        targets_expr = '"${@:-${DEFAULT_TARGETS[@]}}"'
    else:
        targets_expr = '"$@"'

    cmd_parts = [
        "snakemake",
        f"--snakefile {snakefile}",
        f"--configfile {config_abs}",
        f"--directory {output_dir}",
    ]

    if args.profile:
        # snakemake runs with cwd=--directory (output_dir), so a relative
        # profile path would miss. Absolutise here (no symlink follow).
        profile_path = Path(_abs(args.profile))
        # Snakemake interprets path-valued keys inside profile/config.yaml
        # (jobscript, cluster, cluster-status, cluster-cancel, and the v8
        # cluster-generic-* equivalents) relative to cwd, not to the profile
        # dir, so they break under --directory. Rewrite any such bare-name
        # values to absolute paths.
        _absolutise_profile_paths(profile_path)
        cmd_parts.append(f"--profile {profile_path}")
    # snakemake always requires -j: locally it caps parallel rules,
    # under a cluster profile it caps the number of jobs queued at once.
    cmd_parts.append(f"-j {args.jobs}")

    cmd_parts.append("--use-singularity")
    if args.singularity_bind:
        cmd_parts.append(f'--singularity-args "-B {args.singularity_bind} -e"')

    cmd_parts.extend([
        "--rerun-triggers mtime",
        "--rerun-incomplete",
        "--keep-going",
        targets_expr,
    ])

    cmd_lines.append(" \\\n    ".join(cmd_parts))
    cmd_lines.append("")

    runner_path.parent.mkdir(parents=True, exist_ok=True)
    runner_path.write_text("\n".join(cmd_lines))
    runner_path.chmod(runner_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def main():
    parser = argparse.ArgumentParser(
        description="Setup PRCGAP Snakemake workflow (config + runner script)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 1. Pre-stage images (once)
  bash Dockerfile/pull_images.sh   # populates ./images/<tool>.sif

  # 2. Local execution
  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --chm13-fasta chm13.fa \\
    --jobs 8

  # 3. SLURM cluster via profile
  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --chm13-fasta chm13.fa \\
    --profile profile/slurm

  # Override a single image path:
  python setup_workflow.py ... --bam-refiner-image /path/to/custom.sif
        """
    )

    # ---------- Required I/O ----------
    parser.add_argument("--samplesheet", required=True, help="sample sheet TSV")
    parser.add_argument("--chm13-fasta", required=True,
                        help="CHM13 reference FASTA (consumed by copynumber and by the "
                             "INDEL liftvcf_indel_chm13 rule as transanno --query)")

    # ---------- Singularity images ----------
    # Defaults assume Dockerfile/pull_images.sh has populated ./images/.
    # Pass an explicit --*-image to override.
    parser.add_argument("--images-dir", default="images",
                        help="directory containing prepared singularity images (default: images)")
    parser.add_argument("--bam-refiner-image", default=None, help="bam_refiner image (default: <images-dir>/bam_refiner.sif)")
    parser.add_argument("--methylation-image", default=None, help="methylation image (default: <images-dir>/methylation.sif)")
    parser.add_argument("--copynumber-image", default=None, help="copynumber image (default: <images-dir>/copynumber.sif)")
    parser.add_argument("--nanomonsv-image", default=None, help="nanomonsv image (default: <images-dir>/nanomonsv.sif)")
    parser.add_argument("--nanomonsv-postprocess-image", default=None, help="nanomonsv_postprocess image (default: <images-dir>/nanomonsv_postprocess.sif)")
    parser.add_argument("--clairs-image", default=None, help="ClairS image (default: <images-dir>/clairs.sif)")
    parser.add_argument("--deepsomatic-image", default=None, help="DeepSomatic image (default: <images-dir>/deepsomatic.sif)")
    parser.add_argument("--point-mutation-postprocess-image", default=None,
                        help="mutation postprocess image (default: <images-dir>/point_mutation_postprocess.sif)")
    parser.add_argument("--annotation-image", default=None,
                        help="annotation image bundling pysam, samtools, coordconv, "
                             "and transanno — used by every SV/SNV/INDEL annotation "
                             "rule (default: <images-dir>/annotation.sif)")

    # ---------- Resources (per-module overrides) ----------
    res_group = parser.add_argument_group(
        "Per-module resources",
        "Override threads / memory for individual rules. All optional; "
        "omitted values fall back to the defaults shown below.",
    )
    for name, default_threads, default_mem_mb in _RESOURCE_DEFAULTS:
        flag_base = name.replace("_", "-")
        res_group.add_argument(
            f"--{flag_base}-threads",
            type=int, default=default_threads,
            help=f"threads for {name} (default: {default_threads})",
        )
        res_group.add_argument(
            f"--{flag_base}-mem-mb",
            type=int, default=default_mem_mb,
            help=f"memory MB for {name} (default: {default_mem_mb})",
        )

    # ---------- Sample-level params ----------
    parser.add_argument("--sex", choices=["female", "male"], default="female")
    parser.add_argument("--hap1-satellite", default="")
    parser.add_argument("--hap2-satellite", default="")
    parser.add_argument("--simple-repeat", default="")
    parser.add_argument("--gtf-file", default="",
                        help="liftoff GTF (.gtf / .gtf.gz). Consumed by nanomonsv "
                             "insert_classify only.")
    parser.add_argument("--line1-bed", default="")

    # ---------- copy number plot params (all optional) ----------
    parser.add_argument("--copynumber-binwidth", type=float, default=0.05,
                        help="bin width for the depth-ratio mode/ploidy estimate "
                             "(default: 0.05)")
    parser.add_argument("--copynumber-hap1-label", default="Haplotype1",
                        help="haplotype 1 label on the copy number plot "
                             "(default: Haplotype1)")
    parser.add_argument("--copynumber-hap2-label", default="Haplotype2",
                        help="haplotype 2 label on the copy number plot "
                             "(default: Haplotype2)")
    parser.add_argument("--copynumber-ploidy-hap1", default="",
                        help="manual hap1 tumor ploidy override; leave empty to "
                             "auto-estimate")
    parser.add_argument("--copynumber-ploidy-hap2", default="",
                        help="manual hap2 tumor ploidy override; leave empty to "
                             "auto-estimate")
    parser.add_argument("--copynumber-plot-sex-chrom", default="true",
                        choices=["true", "false"],
                        help="force chrX/chrY onto the copy-number plot even when "
                             "absent from the data (default: true; set false to "
                             "show sex chromosomes only when present)")
    parser.add_argument("--chm13-censat", default="",
                        help="CHM13 cenSat v2.1 BED (.gz) for the copy-number plot; "
                             "fills assembly gaps with reference satellite (optional). "
                             "CHM13 chromosome lengths are derived at plot time from "
                             "the --chm13-fasta via chromosome_length.py.")
    parser.add_argument("--pileup-no-baq", default="false",
                        choices=["true", "false"],
                        help="pass --no-BAQ to samtools mpileup in the "
                             "clairs/deepsomatic pileup step (default: false). "
                             "true skips BAQ computation, reducing memory/CPU.")

    # ---------- annotation resources (all optional) ----------
    ann_group = parser.add_argument_group(
        "Annotation resources",
        "All optional; leave empty to skip the corresponding annotation step.",
    )
    ann_group.add_argument("--gff-file", default="",
                           help="tabix-indexed liftoff GFF (.gff.gz + .tbi). Used by SV / SNV / INDEL gene annotation.")
    ann_group.add_argument("--chain-to-grch38", default="",
                           help="chain file from personalized assembly to GRCh38 (SV liftover flag)")
    ann_group.add_argument("--chain-to-chm13", default="",
                           help="chain file from personalized assembly to CHM13 (SV liftover flag)")
    ann_group.add_argument("--repeat-masker-bed", default="",
                           help="tabix-indexed RepeatMasker BED.gz")
    ann_group.add_argument("--segdup-bed", default="",
                           help="tabix-indexed segdup BED.gz")
    ann_group.add_argument("--censat-bed", default="",
                           help="tabix-indexed centromere/satellite BED.gz")
    ann_group.add_argument("--misassembly-hap1-bed", default="",
                           help="hap1 misassembly BED")
    ann_group.add_argument("--misassembly-hap2-bed", default="",
                           help="hap2 misassembly BED")
    ann_group.add_argument("--cancer-gene-census-tsv", default="",
                           help="Cancer Gene Census TSV")
    ann_group.add_argument("--cmrg-gene-list", default="",
                           help="CMRG gene list (one HGNC gene symbol per line; SNV/INDEL gene annotation)")
    ann_group.add_argument("--gencode-transcript-bed", default="",
                           help="GENCODE transcript BED.gz (SNV/INDEL gene annotation)")
    ann_group.add_argument("--gnomad-bed", default="",
                           help="gnomAD SV BED.gz (requires --chain-to-grch38)")
    ann_group.add_argument("--gnomad-vcf", default="",
                           help="gnomAD SNV/INDEL VCF.gz tabix-indexed (requires --chain-to-grch38)")
    ann_group.add_argument("--grch38-fasta", default="",
                           help="GRCh38 reference FASTA (used as transanno --query for INDEL liftover)")

    # ---------- Output / runner ----------
    parser.add_argument("--output-dir", default="results",
                        help="snakemake working directory for results (default: results)")
    parser.add_argument("--output", "-o", default="config/config.yaml",
                        help="path for generated config.yaml (default: config/config.yaml)")
    parser.add_argument("--runner", default="run_workflow.sh",
                        help="path for generated runner shell script (default: run_workflow.sh)")
    parser.add_argument("--workflow-dir", default="workflow",
                        help="snakemake workflow directory containing Snakefile (default: workflow)")
    parser.add_argument("--targets", nargs="+", default=None,
                        help="snakemake targets to bake into the runner so it builds only "
                             "those by default (e.g. 'copynumber/HG008T/output'). Targets are "
                             "relative to --output-dir. Omit to build the full pipeline. "
                             "Targets passed to the runner on the command line still take "
                             "precedence over the baked-in ones.")
    parser.add_argument("--force", "-f", action="store_true", default=False,
                        help="Overwrite existing output files")

    # ---------- Executor / container backend ----------
    parser.add_argument("--profile", default=None,
                        help="snakemake profile path (e.g. profile/slurm). When set, "
                             "snakemake runs jobs through that profile instead of locally.")
    parser.add_argument("--jobs", "-j", type=int, default=8,
                        help="max concurrent jobs (locally: parallel rules; "
                             "with --profile: jobs queued in the cluster). default: 8")
    parser.add_argument("--singularity-bind", default="",
                        help='extra bind paths for singularity, e.g. "/data,/scratch" '
                             "(ignored when --profile is set)")

    args = parser.parse_args()

    config_path = Path(args.output)
    runner_path = Path(args.runner)

    for p in (config_path, runner_path):
        if p.exists() and not args.force:
            print(f"Error: {p} already exists. Use --force to overwrite.")
            sys.exit(1)

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config = create_config(args)

    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    print(f"✓ Configuration written to {config_path}")

    write_runner(args, config_path, runner_path)
    print(f"✓ Runner script written to {runner_path}")

    if args.profile:
        print(f"\nExecutor: snakemake profile {args.profile}")
    else:
        print(f"\nExecutor: local (-j {args.jobs})")
    print("Container backend: singularity")


if __name__ == "__main__":
    main()
