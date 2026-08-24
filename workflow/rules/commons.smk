# import basic packages
import pandas as pd
import os
import re
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

def get_mem_mb(tool, default=32000):
    """Get memory (MB) for a tool from config or use default."""
    return config.get("resources", {}).get(tool, {}).get("mem_mb", default)

def mutation_callers():
    """Point-mutation callers to run: deepsomatic (default), clairs, or both."""
    choice = str(config.get("mutation_caller", "deepsomatic")).strip().lower()
    if choice == "both":
        return ["clairs", "deepsomatic"]
    if choice not in ("clairs", "deepsomatic"):
        raise ValueError(
            "mutation_caller must be 'deepsomatic', 'clairs' or 'both'; got "
            + repr(config.get("mutation_caller")))
    return [choice]


# Every sample must provide at least one sequencing type.
for _sample in samples.index:
    if not sample_seqtypes(_sample):
        raise ValueError(
            "sample " + str(_sample) + " has neither hifi nor ont data in the "
            "sample sheet; at least one is required")


# ---------------------------------------------------------------------------
# In-workflow annotation (dna-brnn / liftoff / chain files)
# ---------------------------------------------------------------------------
# PRCGAP can either consume annotation built elsewhere (the assembly_workflow
# repo) through the config path keys, or build the minimum set itself from the
# assembly plus the CHM13/GRCh38 references. The run_* switches pick per step;
# when a step is on, its generated path wins over the corresponding config key,
# which can then be left empty.
#
# Annotation depends only on the assembly pair, which a tumor and its normal
# usually share, so outputs are keyed on the same canonical sample as k-mer
# extraction (get_kmer_source) and computed once per assembly pair.

ANNOTATION_DIR = "annotation"
MASKED_REF_DIR = "annotation/references"


# dna-brnn and the chain files default on, so a run given only the assembly and
# the CHM13/GRCh38 references produces those itself. liftoff is opt-in: it is the
# heaviest step (8 threads x 128 GB per haplotype) and gene annotation is not
# required by every downstream rule. Set a switch false, or leave run_liftoff
# off, to import that annotation from assembly_workflow through the path keys.
#
# dna-brnn no longer feeds the reference table -- that runs unmasked -- so the
# satellite BEDs exist for the copy-number plot's cen/sat track alone, where a
# contig-coordinate cenSat BED (config censat_bed) supersedes them when given.
def run_dna_brnn():
    """A contig-coordinate cenSat BED supersedes the satellite BEDs wherever
    they are used, so censat_bed turns the step off however the switch is set.
    dna-brnn costs 69 min per haplotype and its output would go unread."""
    if (config.get("censat_bed", "") or ""):
        return False
    return bool(config.get("run_dna_brnn", True))


def run_liftoff():
    return bool(config.get("run_liftoff", False))


def run_chain_files():
    return bool(config.get("run_chain_files", True))


# Both feed nanomonsv and are on by default: neither needs a reference beyond
# the assembly, and together they cost a few minutes per haplotype.
def run_line1():
    return bool(config.get("run_line1", True))


def run_simple_repeat():
    return bool(config.get("run_simple_repeat", True))


def annotation_src(sample):
    """Sample whose generated annotation `sample` should use."""
    return get_kmer_source(sample)


def satellite_bed_src(src, hap):
    """dna-brnn satellite BED.gz for an assembly source sample. hap: hap1|hap2."""
    if run_dna_brnn():
        return "{}/{}/dna_nn/{}.{}_dna-brnn.bed.gz".format(
            ANNOTATION_DIR, src, src, hap)
    return config.get(hap + "_satellite", "") or ""


def satellite_bed(sample, hap):
    return satellite_bed_src(annotation_src(sample), hap)


def liftoff_gff_src(src):
    """Tabix-indexed liftoff GFF for an assembly source sample."""
    if run_liftoff():
        return "{}/{}/liftoff/{}.liftoff.gff.gz".format(ANNOTATION_DIR, src, src)
    return config.get("gff_file", "") or ""


def liftoff_gff(sample):
    return liftoff_gff_src(annotation_src(sample))


def liftoff_gtf_src(src):
    """liftoff GTF.gz for an assembly source sample (nanomonsv insert_classify)."""
    if run_liftoff():
        return "{}/{}/liftoff/{}.liftoff.gtf.gz".format(ANNOTATION_DIR, src, src)
    return config.get("gtf_file", "") or ""


def liftoff_gtf(sample):
    return liftoff_gtf_src(annotation_src(sample))


def line1_bed_src(src):
    """Full-length young LINE-1 BED (nanomonsv insert_classify)."""
    if run_line1():
        return "{}/{}/line1/{}.LINE1.bed.gz".format(ANNOTATION_DIR, src, src)
    return config.get("line1_bed", "") or ""


def line1_bed(sample):
    return line1_bed_src(annotation_src(sample))


def simple_repeat_bed_src(src):
    """Tandem repeat BED (nanomonsv get indel filtering)."""
    if run_simple_repeat():
        return "{}/{}/simple_repeat/{}.simple_repeats.bed.gz".format(
            ANNOTATION_DIR, src, src)
    return config.get("simple_repeat", "") or ""


def simple_repeat_bed(sample):
    return simple_repeat_bed_src(annotation_src(sample))


def chain_file_src(src, ref):
    """assembly → reference chain for an assembly source sample. ref: GRCh38|chm13."""
    if run_chain_files():
        return "{}/{}/chain_files/{}_to_{}.chain".format(
            ANNOTATION_DIR, src, src, ref)
    key = "chain_to_grch38" if ref == "GRCh38" else "chain_to_chm13"
    return config.get(key, "") or ""


def chain_file(sample, ref):
    return chain_file_src(annotation_src(sample), ref)


def as_input(*paths):
    """Input-list form of the resolvers: drops the empty (step disabled) ones.

    Listing a resolved path as a rule input is correct whether it is a
    workflow-generated file (Snakemake builds it) or a user-supplied one
    (Snakemake just requires it to exist).
    """
    return [p for p in paths if p]


# Every run_* step that needs extra reference inputs is checked up front, so a
# missing reference fails at DAG construction rather than hours into the run.
_ANNOTATION_REQUIREMENTS = [
    (run_liftoff, "run_liftoff", ["grch38_fasta", "grch38_gtf"]),
    (run_chain_files, "run_chain_files",
     ["chm13_fasta", "grch38_fasta", "chm13_censat",
      "grch38_centromeres", "grch38_exclusions"]),
]

for _enabled, _flag, _needed in _ANNOTATION_REQUIREMENTS:
    if not _enabled():
        continue
    _missing = [k for k in _needed if not (config.get(k, "") or "")]
    if _missing:
        raise ValueError(
            _flag + " is enabled but these config keys are empty: "
            + ", ".join(_missing)
            + ". Run download_reference.sh and pass the corresponding "
              "setup_workflow.py flags.")
