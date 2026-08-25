#!/usr/bin/env python3
"""
Setup script for PRCGAP Snakemake workflow.
Generates config.yaml and a runner shell script from command-line arguments.
"""

import argparse
import csv
import json
import os
import stat
import sys
import yaml
from pathlib import Path


# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------
# A run is one tumor-normal pair on one assembly: the per-assembly annotation in
# config is global, so a second pair needs its own config and output dir. The
# sheet is therefore written from --tumor/--normal here rather than by a separate
# command; --samplesheet alone still reads a sheet somebody wrote by hand.

COLUMNS = ["sample", "type", "ont", "hifi", "assembly_hap1", "assembly_hap2"]
SEQTYPE_FIELDS = ("ont", "hifi")
PATH_FIELDS = ("ont", "hifi", "assembly_hap1", "assembly_hap2")


def _split_paths(value):
    """Split a (possibly comma-separated) data field into individual paths."""
    return [p.strip() for p in str(value).split(",") if p.strip()]


def _join_paths(paths):
    return ",".join(paths)


def _collect_paths(values):
    """Flatten a repeatable + comma-separated path option into one path list."""
    paths = []
    for value in values or []:
        paths.extend(_split_paths(value))
    return paths


class SampleSheetError(Exception):
    """Raised for invalid sample-sheet input."""


def build_rows_from_pair(args):
    """Turn the --tumor/--normal option group into tumor + normal row dicts.

    The assembly is shared by both samples (one assembly per run), so it is
    given once via --assembly-hap1/--assembly-hap2 and copied into each row.
    """
    missing = [
        flag for flag, value in (
            ("--tumor", args.tumor),
            ("--normal", args.normal),
            ("--assembly-hap1", args.assembly_hap1),
            ("--assembly-hap2", args.assembly_hap2),
        ) if not value
    ]
    if missing:
        raise SampleSheetError(
            "option mode requires " + ", ".join(missing)
            + " (or use --input to validate an existing TSV)")

    specs = [
        ("tumor", args.tumor, args.tumor_ont, args.tumor_hifi),
        ("normal", args.normal, args.normal_ont, args.normal_hifi),
    ]
    rows = []
    for sample_type, name, ont_values, hifi_values in specs:
        ont = _collect_paths(ont_values)
        hifi = _collect_paths(hifi_values)
        # A sample may be HiFi-only or ONT-only; require at least one type.
        if not ont and not hifi:
            raise SampleSheetError(
                f"{sample_type} sample {name}: at least one of "
                f"--{sample_type}-ont / --{sample_type}-hifi is required")
        rows.append({
            "sample": name,
            "type": sample_type,
            "ont": _join_paths(ont),
            "hifi": _join_paths(hifi),
            "assembly_hap1": args.assembly_hap1,
            "assembly_hap2": args.assembly_hap2,
        })
    return rows


def read_input_tsv(path):
    """Read a draft sample sheet TSV into row dicts."""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SampleSheetError(
                f"input TSV {path} is missing required columns: {', '.join(missing)}")
        for raw in reader:
            rows.append({c: (raw.get(c) or "").strip() for c in COLUMNS})
    return rows


def normalize_and_validate(rows, check_exists=True, absolutize=True):
    """Validate rows, absolutise paths, and enforce the one-pair-per-run spec."""
    if not rows:
        raise SampleSheetError("no samples provided (use --tumor/--normal and/or --input)")

    seen = set()
    for row in rows:
        name = row["sample"]
        if not name:
            raise SampleSheetError("a sample has an empty name")
        if name in seen:
            raise SampleSheetError(f"duplicate sample name: {name}")
        seen.add(name)

        if row["type"] not in ("tumor", "normal"):
            raise SampleSheetError(
                f"sample {name}: type must be 'tumor' or 'normal', got '{row['type']}'")

        for field in PATH_FIELDS:
            paths = _split_paths(row[field])
            if not paths:
                # ont / hifi are optional (a sample may be HiFi-only or
                # ONT-only); the "at least one" check below enforces coverage.
                if field in SEQTYPE_FIELDS:
                    row[field] = ""
                    continue
                raise SampleSheetError(f"sample {name}: empty '{field}' field")
            if absolutize:
                paths = [_abs(p) for p in paths]
            if check_exists:
                for p in paths:
                    if not os.path.exists(p):
                        raise SampleSheetError(
                            f"sample {name}: {field} path does not exist: {p} "
                            "(use --no-check-exists if the files live elsewhere)")
            # assembly_* are single files; ont/hifi may be comma-separated.
            if field not in SEQTYPE_FIELDS and len(paths) > 1:
                raise SampleSheetError(
                    f"sample {name}: {field} must be a single file, got {len(paths)}")
            row[field] = _join_paths(paths)

        # Enforce that every sample carries at least one sequencing type.
        if not any(_split_paths(row[f]) for f in SEQTYPE_FIELDS):
            raise SampleSheetError(
                f"sample {name}: at least one of {', '.join(SEQTYPE_FIELDS)} "
                "must be provided")

    normals = [r["sample"] for r in rows if r["type"] == "normal"]
    tumors = [r["sample"] for r in rows if r["type"] == "tumor"]

    # One run = one normal assembly (the per-assembly annotations are global).
    if len(normals) == 0:
        raise SampleSheetError("no normal sample provided (exactly one is required)")
    if len(normals) > 1:
        raise SampleSheetError(
            f"{len(normals)} normal samples ({', '.join(sorted(normals))}); "
            "a run is tied to one normal assembly. Run distinct TN pairs "
            "separately (own config + sample sheet + output dir).")
    if not tumors:
        print("Warning: no tumor sample provided.", file=sys.stderr)

    assemblies = {(r["assembly_hap1"], r["assembly_hap2"]) for r in rows}
    if len(assemblies) > 1:
        raise SampleSheetError(
            "samples reference more than one assembly pair; the workflow's "
            "annotation resources are global (one assembly per run). All samples "
            "must share the same assembly_hap1/assembly_hap2. Distinct assemblies "
            "must be run separately.")

    return rows


def write_sample_sheet(rows, output):
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([row.get(col, "") for col in COLUMNS])


# Directory containing this script (the PRCGAP repo root). Used to default
# repo-internal paths (workflow/, images/) so they resolve correctly no matter
# which working directory setup_workflow.py is invoked from — relative defaults
# would otherwise be resolved against the CWD by _abs().
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _image_path(explicit, images_dir, tool):
    """Return user-supplied image path or fall back to <images_dir>/<tool>.sif."""
    if explicit:
        return explicit
    return os.path.join(images_dir, f"{tool}.sif")


# (module_name, default_threads, default_mem_mb) for each rule whose
# resources can be overridden from the CLI.
_RESOURCE_DEFAULTS = [
    # In-workflow annotation (opt-in via --run-dna-brnn / --run-liftoff /
    # --run-chain-files).
    # In-workflow annotation. Every rule is capped at 8 threads. liftoff keeps
    # the 128 GB assembly_workflow declares for it (unmeasured here, and it did
    # need that much); make_chain_files runs the same minimap2 asm5 as
    # copynumber_ref_table, measured at 34.6 GB. dna_brnn is 8 GB against a
    # measured 4.4 GB -- it only annotates satellite now that ref.table runs
    # unmasked, so the plot's cen/sat track is its one consumer, and censat_bed
    # supersedes it there (69 min per haplotype).
    ("dna_brnn", 8, 8000),   # measured 4.4 GB / 69 min per haplotype
    ("liftoff_reference", 1, 8000),
    ("liftoff", 8, 96000),
    ("liftoff_merge", 1, 16000),
    ("prepare_mask_regions", 1, 10240),
    ("make_chain_files", 8, 48000),             # same minimap2 asm5 shape as copynumber_ref_table, 34.6 GB there
    ("chain_files_merge", 1, 8000),
    # LINE-1 and tandem repeats: both are blastn/ULTRA over one haplotype, far
    # lighter than the RepeatMasker route they replace.
    ("line1", 8, 16000),
    ("line1_merge", 1, 8000),
    ("simple_repeat", 8, 32000),
    ("simple_repeat_merge", 1, 8000),
    ("bam_refiner_kmer", 8, 32000),  # measured 15.2 GB (H2009); HG008 2.3 GB
    ("bam_refiner", 8, 64000),                  # measured 37.8 GB
    ("assembly_bwa_index", 1, 16000),
    ("methylation", 8, 32000),                  # measured 13.5 GB
    ("copynumber", 8, 48000),                   # legacy fallback; matches copynumber_ref_table, the heaviest split rule
    ("copynumber_reference", 1, 4000),    # measured 0.07 GB
    ("copynumber_ref_table", 8, 48000),   # minimap2 unmasked: 34.1 GB at 8 threads
    ("copynumber_ref_table_final", 1, 8000),
    ("copynumber_depth", 4, 8000),        # mosdepth: measured 4.3 GB, flat in threads
    ("copynumber_segment", 1, 4000),      # measured 188 MB (cbs.R + split_gaps.py)
    ("copynumber_plot", 1, 4000),         # measured 0.23 GB
    ("nanomonsv_parse", 8, 16000),              # measured 3.1 GB
    ("nanomonsv_get", 8, 64000),                # measured 2.2 GB; --max_memory_minimap2 16 caps the aligner
    ("nanomonsv_postprocess", 1, 8000),         # measured 0.1 GB
    ("nanomonsv_insert_classify", 8, 64000),    # measured 54.9 GB -- the tightest of these
    ("nanomonsv_connect", 1, 30000),
    ("nanomonsv_merge", 1, 30000),
    # Callers are scattered over contig chunks: these are per-chunk resources.
    # 8 threads / 32000 MB is the benchmarked setting
    # (plan/mutation_calling_performance.md 6.6).
    #
    # --num_shards is the rule's thread count, and each shard is a separate
    # TensorFlow process. Measured on an unbound 192-core host, where every
    # library sizes itself from the host core count: 8 shards reached >58.5 GB
    # of address space against a 32 GB s_vmem grant and the process group was
    # killed with SIGXCPU ~10 s in; 2 shards completed at 31.09 GB, i.e. 3 %
    # under the same grant. RSS was never the problem -- 18 MB at kill time.
    # 8 is kept here because profile/sge/cluster.yaml now core-binds every job,
    # which makes nproc report the allocation rather than 192. Re-measure
    # maxvmem before trusting this on a scheduler without such binding.
    ("clairs", 8, 32000),
    ("clairs_scatter_setup", 8, 16000),
    ("clairs_chunks", 1, 4000),
    ("clairs_merge", 1, 16000),
    ("deepsomatic", 8, 32000),
    ("deepsomatic_scatter_setup", 8, 16000),
    ("deepsomatic_chunks", 1, 4000),
    ("deepsomatic_merge", 1, 16000),
    ("clairs_postprocess", 1, 32000),
    ("clairs_postprocess_split", 1, 32000),
    ("clairs_postprocess_realign", 8, 16000),  # measured 0.1 GB (bwa on the candidate-flank FASTA)
    ("clairs_postprocess_pileup", 8, 32000),
    ("clairs_postprocess_haplotype", 4, 64000),
    ("deepsomatic_postprocess", 1, 32000),
    ("deepsomatic_postprocess_split", 1, 32000),
    ("deepsomatic_postprocess_realign", 8, 16000),  # measured 0.1 GB (bwa on the candidate-flank FASTA)
    ("deepsomatic_postprocess_pileup", 8, 32000),
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


# Rules split out of a formerly monolithic one inherit the parent rule's CLI
# resources unless their own flag was passed, so an existing invocation that only
# sets --copynumber-threads / --copynumber-mem-mb still sizes every split job.
_RESOURCE_INHERIT = {
    "copynumber_ref_table": "copynumber",
    "copynumber_depth": "copynumber",
    "copynumber_segment": "copynumber",
    "copynumber_plot": "copynumber",
}


def _apply_resource_inheritance(args):
    """Copy a parent rule's resources onto children whose flags were defaulted."""
    declared = {name: (threads, mem_mb) for name, threads, mem_mb in _RESOURCE_DEFAULTS}
    for child, parent in _RESOURCE_INHERIT.items():
        if child not in declared or parent not in declared:
            continue
        for idx, attr in ((0, "threads"), (1, "mem_mb")):
            parent_val = getattr(args, f"{parent}_{attr}")
            if parent_val == declared[parent][idx]:
                continue  # parent left at its default; nothing to propagate
            if getattr(args, f"{child}_{attr}") != declared[child][idx]:
                continue  # child was set explicitly; keep it
            setattr(args, f"{child}_{attr}", parent_val)


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
            "dna_nn": _abs(_image_path(args.dna_nn_image, args.images_dir, "dna_nn")),
            "liftoff": _abs(_image_path(args.liftoff_image, args.images_dir, "liftoff")),
            "chain_files": _abs(_image_path(args.chain_files_image, args.images_dir, "chain_files")),
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
        "pileup_max_depth": args.pileup_max_depth,
        "pileup_num_chunks": args.pileup_num_chunks,
        # ---- optional nanomonsv steps (not used in the paper) ----
        "run_nanomonsv_connect": args.run_nanomonsv_connect,
        "run_nanomonsv_merge": args.run_nanomonsv_merge,
        # ---- point-mutation caller ----
        "mutation_caller": args.mutation_caller,
        "caller_solo_contig_min_bp": args.caller_solo_contig_min_bp,
        "deepsomatic_postprocess_variants_cpus": args.deepsomatic_postprocess_variants_cpus,
        # ---- ClairS platform model ----
        "clairs_model": args.clairs_model,
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
        # ---- in-workflow annotation (opt-in per step) ----
        "run_dna_brnn": args.run_dna_brnn,
        "run_liftoff": args.run_liftoff,
        "run_chain_files": args.run_chain_files,
        "run_line1": args.run_line1,
        "run_simple_repeat": args.run_simple_repeat,
        "grch38_gtf": _abs(args.grch38_gtf) or "",
        "grch38_centromeres": _abs(args.grch38_centromeres) or "",
        "grch38_exclusions": _abs(args.grch38_exclusions) or "",
    }
    return config

def _is_slurm_profile(profile: str) -> bool:
    """True when the profile path looks like a SLURM profile (name contains
    'slurm'). Used to decide whether the runner needs #SBATCH headers."""
    return bool(profile) and "slurm" in os.path.basename(os.path.normpath(profile)).lower()


def _profile_partition(profile_dir: Path):
    """Extract the partition from <profile>/settings.json SBATCH_DEFAULTS.

    Returns the partition string (e.g. 'mjobs,rjobs') or None when it cannot be
    determined (missing settings.json, absent/empty partition key). The runner
    omits the '#SBATCH -p' line in that case, deferring to SLURM's default
    partition.
    """
    settings = profile_dir / "settings.json"
    if not settings.is_file():
        return None
    try:
        defaults = json.loads(settings.read_text()).get("SBATCH_DEFAULTS", "")
    except (json.JSONDecodeError, OSError):
        return None
    for tok in defaults.split():
        if tok.startswith("partition="):
            part = tok.split("=", 1)[1]
            return part or None
    return None


def _is_sge_profile(profile: str) -> bool:
    """True when the profile path looks like an SGE/UGE profile (name contains
    'sge' or 'uge'). Used to decide whether the runner needs #$ headers."""
    if not profile:
        return False
    name = os.path.basename(os.path.normpath(profile)).lower()
    return "sge" in name or "uge" in name


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

    # `set -euo pipefail` is emitted below the header block, not here: both
    # sbatch and qsub stop scanning for directives at the first command, so any
    # header after a command is silently ignored.
    cmd_lines = ["#!/bin/bash", ""]

    # When the runner drives a SLURM profile, make the driver itself a batch job
    # by emitting #SBATCH headers directly below the shebang. The partition is
    # read from the profile's settings.json so it stays defined in one place; if
    # it can't be determined, the -p line is omitted and SLURM's default
    # partition applies.
    is_cluster_driver = False
    if _is_slurm_profile(args.profile):
        is_cluster_driver = True
        cmd_lines += [
            "#SBATCH -c 1",
            f"#SBATCH --mem-per-cpu={args.driver_mem}G",
            "#SBATCH -J PRCGAP",
            # The driver only orchestrates, but it has to outlive the whole
            # workflow. Without --time it inherits the partition's DefaultTime
            # (3 days here), which killed a full-size run at 72%.
            f"#SBATCH --time={args.driver_time}",
        ]
        partition = _profile_partition(Path(_abs(args.profile)))
        if partition:
            cmd_lines.append(f"#SBATCH -p {partition}")
        cmd_lines.append("#SBATCH -o log/%x_%j.log -e log/%x_%j.log")

    # Same idea for SGE/UGE. There is no walltime line to match the SLURM
    # branch's --time: s_rt/h_rt are INFINITY on the queues this targets, and a
    # guessed limit would kill the driver mid-workflow rather than protect it.
    elif _is_sge_profile(args.profile):
        is_cluster_driver = True
        cmd_lines += [
            # posix_compliant queues run the script with /bin/sh unless told
            # otherwise, and the snakemake command below is bash.
            "#$ -S /bin/bash",
            "#$ -cwd",
            "#$ -N PRCGAP",
            "#$ -o log/ -e log/",
            f"#$ -l s_vmem={args.driver_mem}G",
        ]

    if is_cluster_driver:
        # Both header blocks send the driver's own stdout/stderr into log/, and
        # the scheduler will not create it: qsub/sbatch fail the job instead.
        (runner_path.parent / "log").mkdir(parents=True, exist_ok=True)

    # Site-specific setup (module load, PATH for snakemake, ...). Emitted before
    # `set -u` on purpose: the environment-modules init scripts reference unset
    # variables and would abort the runner under it.
    if args.runner_preamble:
        cmd_lines.append("")
        cmd_lines.append("# Site setup from setup_workflow.py --runner-preamble.")
        cmd_lines += list(args.runner_preamble)
        cmd_lines.append("")

    cmd_lines += ["set -euo pipefail", ""]

    if _is_sge_profile(args.profile):
        # The driver only orchestrates, but numpy/OpenBLAS (pulled in via
        # pandas) sizes its thread-local buffers from the host's core count. On
        # a 192-core exec host that overran the driver's s_vmem grant before
        # snakemake had built the DAG ("OpenBLAS error: Memory allocation still
        # failed after 10 retries"). SGE's s_vmem is a virtual-memory limit, so
        # the reservation counts even though the driver's RSS stays tiny.
        cmd_lines += [
            "# One BLAS thread is plenty for a process that only submits jobs.",
            "export OPENBLAS_NUM_THREADS=1",
            "export OMP_NUM_THREADS=1",
            "export MKL_NUM_THREADS=1",
            "export NUMEXPR_NUM_THREADS=1",
            "",
        ]

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
  # 1. One-time setup
  bash images/pull_images.sh                    # populates ./images/<tool>.sif
  bash resource/scripts/download_reference.sh   # references + LINE-1 models

  # 2. A tumor-normal pair. The --tumor/--normal options WRITE --samplesheet;
  #    a run is one pair on one assembly, so a second pair needs its own config
  #    and --output-dir.
  python setup_workflow.py \\
    --samplesheet config/samples.tsv \\
    --tumor  T --tumor-hifi  t.hifi.bam --tumor-ont  t.ont.bam \\
    --normal N --normal-hifi n.hifi.bam --normal-ont n.ont.bam \\
    --assembly-hap1 asm/hap1.fa --assembly-hap2 asm/hap2.fa \\
    --chm13-fasta resource/reference/chm13v2.0_maskedY_rCRS.fa \\
    --jobs 8

  # HiFi-only or ONT-only: omit the other flag. Repeat a flag (or comma-separate)
  # for several files per type:
  #   --tumor-hifi t.1.bam --tumor-hifi t.2.bam

  # 3. Build the annotation instead of importing it, so nothing but reads, the
  #    assembly and the references is needed:
  python setup_workflow.py \\
    --samplesheet config/samples.tsv \\
    --tumor T --tumor-hifi t.bam --normal N --normal-hifi n.bam \\
    --assembly-hap1 asm/hap1.fa --assembly-hap2 asm/hap2.fa \\
    --run-dna-brnn --run-liftoff --run-chain-files \\
    --run-line1 --run-simple-repeat \\
    --chm13-fasta        resource/reference/chm13v2.0_maskedY_rCRS.fa \\
    --grch38-fasta       resource/reference/GRCh38.d1.vd1.fa \\
    --grch38-gtf         resource/reference/Homo_sapiens.GRCh38.Ensembl.112.chr.format.gtf \\
    --grch38-centromeres resource/reference/centromeres.txt.gz \\
    --grch38-exclusions  resource/reference/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed \\
    --chm13-censat       resource/reference/chm13v2.0_censat_v2.1.bed.gz \\
    --profile profile/slurm

  # 4. Reuse a sheet somebody wrote by hand: pass --samplesheet with no pair
  #    options and it is read as-is.

  # Override a single image path:
  python setup_workflow.py ... --bam-refiner-image /path/to/custom.sif
        """
    )

    # ---------- Required I/O ----------
    parser.add_argument("--samplesheet", required=True,
                        help="sample sheet TSV. Written from the pair options below "
                             "when any of them is given, otherwise read as-is.")

    pair = parser.add_argument_group(
        "sample sheet (writes --samplesheet instead of reading it)",
        "Name the tumor and normal and point each at its ONT and/or HiFi data. "
        "Each is optional but a sample needs at least one, so HiFi-only and "
        "ONT-only pairs work. Repeatable, and comma-separated lists are accepted.")
    pair.add_argument("--tumor", metavar="NAME", help="tumor sample name")
    pair.add_argument("--tumor-ont", action="append", metavar="PATH",
                      help="tumor ONT data (fastq.gz or bam). Repeatable.")
    pair.add_argument("--tumor-hifi", action="append", metavar="PATH",
                      help="tumor HiFi data (fastq.gz or bam). Repeatable.")
    pair.add_argument("--normal", metavar="NAME", help="normal sample name")
    pair.add_argument("--normal-ont", action="append", metavar="PATH",
                      help="normal ONT data (fastq.gz or bam). Repeatable.")
    pair.add_argument("--normal-hifi", action="append", metavar="PATH",
                      help="normal HiFi data (fastq.gz or bam). Repeatable.")
    pair.add_argument("--assembly-hap1", metavar="FASTA",
                      help="assembly hap1 fasta, shared by both samples")
    pair.add_argument("--assembly-hap2", metavar="FASTA",
                      help="assembly hap2 fasta, shared by both samples")
    pair.add_argument("--no-check-exists", dest="check_exists",
                      action="store_false", default=True,
                      help="do not check that the data files exist yet")
    pair.add_argument("--no-absolutize", dest="absolutize",
                      action="store_false", default=True,
                      help="keep the data paths as given instead of absolutising them")
    parser.add_argument("--chm13-fasta", required=True,
                        help="CHM13 reference FASTA (consumed by copynumber and by the "
                             "INDEL liftvcf_indel_chm13 rule as transanno --query)")

    # ---------- Singularity images ----------
    # Defaults assume Dockerfile/pull_images.sh has populated ./images/.
    # Pass an explicit --*-image to override.
    parser.add_argument("--images-dir", default=os.path.join(_SCRIPT_DIR, "images"),
                        help="directory containing prepared singularity images "
                             "(default: the images/ dir next to setup_workflow.py, so it "
                             "works regardless of the current working directory)")
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
    parser.add_argument("--dna-nn-image", default=None,
                        help="dna-brnn image (default: <images-dir>/dna_nn.sif); "
                             "only used with --run-dna-brnn")
    parser.add_argument("--liftoff-image", default=None,
                        help="liftoff + minimap2 + gffread image (default: "
                             "<images-dir>/liftoff.sif); only used with --run-liftoff")
    parser.add_argument("--chain-files-image", default=None,
                        help="minimap2 + transanno + chaintools + rustybam + "
                             "paf2chain image (default: <images-dir>/chain_files.sif); "
                             "only used with --run-chain-files")

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
    parser.add_argument("--pileup-no-baq", default="true",
                        choices=["true", "false"],
                        help="pass --no-BAQ to samtools mpileup in the "
                             "clairs/deepsomatic pileup step (default: true). "
                             "true skips BAQ computation, cutting peak RSS ~26x "
                             "for long-read/high-depth data; set false to keep BAQ.")
    parser.add_argument("--pileup-max-depth", type=int, default=0,
                        help="samtools mpileup -d/--max-depth for the pileup step "
                             "(default 0 = samtools default 8000). Set >0 to cap "
                             "per-file depth in very high-depth regions.")
    parser.add_argument("--pileup-num-chunks", type=int, default=16,
                        help="target maximum number of pileup chunks produced by "
                             "split_bed.sh (default 16). Contigs are packed into "
                             "balanced chunks pileup'd with `samtools mpileup -r`.")
    parser.add_argument("--run-nanomonsv-connect", action="store_true", default=False,
                        help="run the optional nanomonsv connect step "
                             "(not used in the paper; default off).")
    parser.add_argument("--run-nanomonsv-merge", action="store_true", default=False,
                        help="run the optional nanomonsv merge step combining "
                             "HiFi+ONT results (not used in the paper; default "
                             "off; only runs when both seqtypes are present).")
    parser.add_argument("--mutation-caller", default="deepsomatic",
                        choices=["deepsomatic", "clairs", "both"],
                        help="point-mutation caller to run (default: "
                             "deepsomatic). 'both' runs ClairS and DeepSomatic "
                             "together, as earlier versions always did.")
    parser.add_argument("--driver-time", default="14-00:00:00",
                        help="#SBATCH --time for the snakemake driver job when a "
                             "SLURM profile is used (default 14-00:00:00). The "
                             "driver must outlive the whole workflow; without it "
                             "the partition's DefaultTime applies.")
    parser.add_argument("--driver-mem", type=int, default=8,
                        help="memory in GB for the snakemake driver job when a "
                             "SLURM or SGE profile is used (default 8). Becomes "
                             "'#SBATCH --mem-per-cpu' or '#$ -l s_vmem'. The "
                             "driver only submits jobs, so this is about the "
                             "python process, not the rules.")
    parser.add_argument("--runner-preamble", action="append", default=[],
                        metavar="LINE",
                        help="shell line to emit in the runner before the "
                             "snakemake call; repeatable. For site setup the "
                             "generated script cannot know, e.g. "
                             "--runner-preamble 'module load apptainer' "
                             "--runner-preamble 'export PATH=~/miniconda3/envs/"
                             "snakemake7/bin:$PATH'. Emitted above "
                             "'set -euo pipefail' so module init scripts that "
                             "read unset variables do not abort the runner.")
    parser.add_argument("--caller-solo-contig-min-bp", type=int, default=1000000,
                        help="contigs at least this long get their own caller "
                             "chunk/job; shorter ones are bundled into one "
                             "(default 1000000).")
    parser.add_argument("--deepsomatic-postprocess-variants-cpus", type=int, default=1,
                        help="worker count for DeepSomatic's postprocess_variants "
                             "(default 1). Its own default is the node's core "
                             "count and ignores the scheduler allocation.")
    parser.add_argument("--clairs-model", default="hifi_sequel2",
                        choices=["hifi_sequel2", "hifi_revio",
                                 "ont_r10_dorado_sup_5khz_ssrs",
                                 "ont_r10_dorado_sup_4khz"],
                        help="ClairS platform model (default: hifi_sequel2). "
                             "Choose to match your data; applies to every "
                             "seqtype ClairS runs on.")

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
                           help="GRCh38 reference FASTA. transanno --query for INDEL "
                                "liftover, and the source for --run-liftoff / "
                                "--run-chain-files.")

    # ---------- in-workflow annotation generation (opt-in per step) ----------
    gen_group = parser.add_argument_group(
        "In-workflow annotation generation",
        "Build the minimum annotation set inside PRCGAP instead of importing it "
        "from assembly_workflow. Each switch supersedes the corresponding path "
        "flag above, which can then be omitted. Run download_reference.sh to "
        "fetch the reference inputs.",
    )
    # Satellite and chain generation are on by default, so a run given only the
    # assembly and the CHM13/GRCh38 references produces those itself; pass
    # --no-run-<step> to import them from assembly_workflow instead. liftoff is
    # off by default and has to be asked for.
    gen_group.add_argument("--run-dna-brnn", action=argparse.BooleanOptionalAction,
                           default=True,
                           help="build the per-haplotype dna-brnn satellite BEDs "
                                "(supersedes --hap1-satellite / --hap2-satellite); "
                                "default on, --no-run-dna-brnn to disable")
    gen_group.add_argument("--run-liftoff", action=argparse.BooleanOptionalAction,
                           default=False,
                           help="build the liftoff gene GFF + GTF (supersedes "
                                "--gff-file / --gtf-file); requires --grch38-fasta "
                                "and --grch38-gtf; default off (heaviest step: "
                                "8 threads x 128 GB per haplotype), pass "
                                "--run-liftoff to enable")
    gen_group.add_argument("--run-chain-files", action=argparse.BooleanOptionalAction,
                           default=True,
                           help="build the assembly → CHM13/GRCh38 chain files "
                                "(supersedes --chain-to-chm13 / --chain-to-grch38); "
                                "requires --grch38-fasta, --chm13-censat, "
                                "--grch38-centromeres and --grch38-exclusions; "
                                "default on, --no-run-chain-files to disable")
    gen_group.add_argument("--run-line1", action=argparse.BooleanOptionalAction,
                           default=True,
                           help="build the full-length young LINE-1 BED that "
                                "nanomonsv insert_classify uses (supersedes "
                                "--line1-bed); needs no extra reference, the "
                                "L1.3 query and Dfam subunit models ship with "
                                "the workflow; default on, --no-run-line1 to "
                                "disable")
    gen_group.add_argument("--run-simple-repeat", action=argparse.BooleanOptionalAction,
                           default=True,
                           help="build the tandem repeat BED that nanomonsv get "
                                "uses to filter indel-like SVs (supersedes "
                                "--simple-repeat); default on, "
                                "--no-run-simple-repeat to disable")
    gen_group.add_argument("--grch38-gtf", default="",
                           help="GRCh38 GTF with chr* contig names (plain .gtf) for "
                                "--run-liftoff")
    gen_group.add_argument("--grch38-centromeres", default="",
                           help="UCSC hg38 centromeres.txt(.gz) for --run-chain-files")
    gen_group.add_argument("--grch38-exclusions", default="",
                           help="GRC exclusion regions BED for --run-chain-files")

    # ---------- Output / runner ----------
    parser.add_argument("--output-dir", default="results",
                        help="snakemake working directory for results (default: results)")
    parser.add_argument("--output", "-o", default="config/config.yaml",
                        help="path for generated config.yaml (default: config/config.yaml)")
    parser.add_argument("--runner", default="run_workflow.sh",
                        help="path for generated runner shell script (default: run_workflow.sh)")
    parser.add_argument("--workflow-dir", default=os.path.join(_SCRIPT_DIR, "workflow"),
                        help="snakemake workflow directory containing Snakefile "
                             "(default: the workflow/ dir next to setup_workflow.py, so "
                             "it works regardless of the current working directory)")
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
    _apply_resource_inheritance(args)

    config_path = Path(args.output)
    runner_path = Path(args.runner)
    sheet_path = Path(args.samplesheet)

    # Any pair option means the sheet is an output, not an input.
    writes_sheet = any((args.tumor, args.normal, args.assembly_hap1, args.assembly_hap2,
                        args.tumor_ont, args.tumor_hifi, args.normal_ont, args.normal_hifi))

    outputs = [config_path, runner_path] + ([sheet_path] if writes_sheet else [])
    for p in outputs:
        if p.exists() and not args.force:
            print(f"Error: {p} already exists. Use --force to overwrite.")
            sys.exit(1)

    if writes_sheet:
        try:
            rows = normalize_and_validate(build_rows_from_pair(args),
                                          check_exists=args.check_exists,
                                          absolutize=args.absolutize)
        except (SampleSheetError, FileNotFoundError) as exc:
            parser.error(str(exc))
        sheet_path.parent.mkdir(parents=True, exist_ok=True)
        write_sample_sheet(rows, str(sheet_path))
        n_tumor = sum(1 for r in rows if r["type"] == "tumor")
        print(f"✓ Sample sheet written to {sheet_path}: {len(rows)} samples "
              f"({n_tumor} tumor, {len(rows) - n_tumor} normal).")
    elif not sheet_path.exists():
        parser.error(f"{sheet_path} does not exist; give the --tumor/--normal "
                     "options to write it")

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
