# import basic packages
import pandas as pd
import os
import sys
from snakemake.utils import validate

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# The ont / hifi columns are optional per sample: a sample may carry only one
# sequencing type (HiFi-only or ONT-only). Empty cells are read as NaN by
# pandas; normalise them to empty strings so the schema (type: string)
# validates and the seqtype-availability helpers below can test them uniformly.
for _seqcol in ("ont", "hifi"):
    if _seqcol in samples.columns:
        samples[_seqcol] = samples[_seqcol].fillna("").astype(str)

tumors = samples.loc[samples["type"] == "tumor", "sample"].tolist()
normals = samples.loc[samples["type"] == "normal", "sample"].tolist()

# validate sample sheet and config file
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

# Helper functions for getting sample metadata
def get_sample_hifi(sample):
    """Return comma-separated HiFi file paths for a sample."""
    return samples.loc[sample, "hifi"]

def get_sample_ont(sample):
    """Return comma-separated ONT file paths for a sample."""
    return samples.loc[sample, "ont"]

def get_sample_files_list(sample, seqtype):
    """Return list of individual file paths for a sample and seqtype.
    Handles comma-separated paths in the sample sheet.
    """
    raw = samples.loc[sample, seqtype]
    return [f.strip() for f in raw.split(",")]

def has_multiple_inputs(sample, seqtype):
    """Check if a sample has multiple input files for a given seqtype."""
    return len(get_sample_files_list(sample, seqtype)) > 1

def get_first_input_extension(sample, seqtype):
    """Get the file extension of the first input file."""
    files = get_sample_files_list(sample, seqtype)
    return files[0].split(".")[-1].lower()

# ---------------------------------------------------------------------------
# Sequencing-type availability
# ---------------------------------------------------------------------------
# A sample may provide HiFi-only, ONT-only, or both. These helpers drive which
# per-seqtype rules fire and how the analysis targets are enumerated so the
# workflow no longer requires both sequencing types to be present.

def has_seqtype(sample, seqtype):
    """True if the sample sheet lists non-empty data for this seqtype."""
    val = samples.loc[sample, seqtype]
    if not isinstance(val, str):
        return False
    return val.strip() not in ("", "nan", "NA", "NaN")

def sample_seqtypes(sample):
    """Seqtypes present for a sample, in preference order (hifi first)."""
    return [s for s in ("hifi", "ont") if has_seqtype(sample, s)]

def primary_seqtype(sample):
    """Preferred single seqtype for a sample (hifi if present, else ont)."""
    seqtypes = sample_seqtypes(sample)
    if not seqtypes:
        raise ValueError(
            "sample " + str(sample) + " has neither hifi nor ont data")
    return seqtypes[0]

def paired_seqtypes(tumor):
    """Seqtypes available for BOTH a tumor and its paired normal.

    Tumor/normal callers (nanomonsv get, clairs, deepsomatic, copynumber)
    consume matched tumor+normal BAMs of the same seqtype, so a seqtype is only
    analysable for a tumor when its paired normal also has that seqtype.
    """
    normal = get_paired_normal(tumor)
    tset = set(sample_seqtypes(tumor))
    nset = set(sample_seqtypes(normal))
    return [s for s in ("hifi", "ont") if s in tset and s in nset]

def primary_paired_seqtype(tumor):
    """Preferred single seqtype available for both tumor and its normal."""
    seqtypes = paired_seqtypes(tumor)
    if not seqtypes:
        raise ValueError(
            "tumor " + str(tumor) + " and its paired normal share no common "
            "seqtype (hifi/ont)")
    return seqtypes[0]

def get_sample_assembly_hap1(sample):
    return samples.loc[sample, "assembly_hap1"]

def get_sample_assembly_hap2(sample):
    return samples.loc[sample, "assembly_hap2"]

def get_kmer_source(sample):
    """Return the canonical sample name whose kmer outputs `sample` should use.

    Tumor and normal frequently share an assembly fasta in this workflow, so
    haplotype-specific k-mer extraction can be deduplicated. Among all samples
    that share the same (assembly_hap1, assembly_hap2) pair, the lexicographically
    smallest sample name is the canonical kmer source. Other samples reference
    that sample's bam_refiner/{kmer_src}/kmer/ outputs, and Snakemake's DAG
    only fires bam_refiner_kmer for the canonical samples.
    """
    hap1 = samples.loc[sample, "assembly_hap1"]
    hap2 = samples.loc[sample, "assembly_hap2"]
    candidates = samples[
        (samples["assembly_hap1"] == hap1) & (samples["assembly_hap2"] == hap2)
    ].index.tolist()
    return min(candidates)


def get_paired_normal(tumor):
    """Get the paired normal sample for a tumor sample.

    If there's only one normal, return it.
    Otherwise, look for a 'paired_normal' column in the sample sheet.
    """
    if len(normals) == 1:
        return normals[0]
    elif "paired_normal" in samples.columns:
        return samples.loc[tumor, "paired_normal"]
    else:
        raise ValueError("Cannot determine paired normal for tumor " + tumor + ". "
                        "Either provide only one normal sample or add 'paired_normal' column to sample sheet.")

# Default resource configurations for each tool
# These can be overridden in config.yaml under 'resources' key
def get_threads(tool, default=8):
    """Get thread count for a tool from config or use default."""
    return config.get("resources", {}).get(tool, {}).get("threads", default)

def mutation_callers():
    """Point-mutation callers to run: DeepSomatic always, ClairS when run_clairs."""
    if str(config.get("run_clairs", False)).strip().lower() == "true":
        return ["clairs", "deepsomatic"]
    return ["deepsomatic"]

MUTATION_SEQTYPES_CHOICES = ("primary", "hifi", "ont", "all")

def mutation_seqtypes(tumor):
    """Seqtypes the point-mutation callers (DeepSomatic, and ClairS when
    run_clairs) run on for a tumor, per config["mutation_seqtypes"]
    (same key and values as the devel branch):

      primary (default) -> one seqtype, HiFi preferred (primary_paired_seqtype)
      hifi / ont        -> that seqtype only
      all               -> every paired seqtype

    Only seqtypes present for both the tumor and its paired normal are
    eligible; a requested seqtype the pair lacks is skipped with a warning.
    Both callers share the list, so the point_mutation_other cross-check
    always finds the other caller's VCF for the same seqtype.
    """
    mode = str(config.get("mutation_seqtypes", "primary")).strip().lower()
    if mode not in MUTATION_SEQTYPES_CHOICES:
        raise ValueError(
            "mutation_seqtypes must be one of " + ", ".join(MUTATION_SEQTYPES_CHOICES)
            + " (got " + repr(mode) + ")")
    paired = paired_seqtypes(tumor)
    if mode == "primary":
        return paired[:1]
    if mode == "all":
        return paired
    if mode not in paired:
        sys.stderr.write(
            "WARNING: mutation_seqtypes=" + mode + " but tumor " + str(tumor)
            + " and its paired normal have no " + mode
            + " data; skipping point-mutation calling for this tumor.\n")
        return []
    return [mode]

def clairs_platform(seqtype):
    """ClairS --platform for a seqtype (config clairs_platform_hifi / _ont).
    There is no default: the right model depends on the instrument (Sequel II
    vs Revio) and ONT chemistry/basecaller (R9.4.1 vs R10.4.1, guppy vs dorado,
    4 vs 5 kHz), so it must be set for every seqtype ClairS runs on."""
    key = "clairs_platform_" + seqtype
    value = str(config.get(key, "") or "").strip()
    if not value:
        raise ValueError(
            "run_clairs is on and ClairS calls on " + seqtype + " data, but "
            + key + " is not set (setup_workflow.py --clairs-platform-" + seqtype
            + "). Choose the ClairS platform matching the " + seqtype + " data.")
    if not value.startswith(seqtype + "_"):
        raise ValueError(
            key + "=" + value + " is not a " + seqtype + " ClairS platform")
    return value

def check_clairs_platforms(tumors):
    """Fail at DAG build if a seqtype ClairS will run on has no platform set."""
    if "clairs" not in mutation_callers():
        return
    if config.get("clairs_model"):
        raise ValueError(
            "clairs_model was replaced by clairs_platform_hifi / "
            "clairs_platform_ont (one ClairS platform per seqtype); "
            "re-run setup_workflow.py or edit the config.")
    for seqtype in sorted({s for t in tumors for s in mutation_seqtypes(t)}):
        clairs_platform(seqtype)

def pileup_no_baq():
    """'true'/'false' for pileup.sh: pass --no-BAQ to samtools mpileup (default off)."""
    return "true" if str(config.get("pileup_no_baq", False)).strip().lower() == "true" else "false"

def get_mem_mb(tool, default=32000):
    """Get memory (MB) for a tool from config or use default."""
    return config.get("resources", {}).get(tool, {}).get("mem_mb", default)


# Every sample must provide at least one sequencing type.
for _sample in samples.index:
    if not sample_seqtypes(_sample):
        raise ValueError(
            "sample " + str(_sample) + " has neither hifi nor ont data in the "
            "sample sheet; at least one is required")
