# import basic packages
import pandas as pd
import os
from snakemake.utils import validate

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

tumors = samples.loc[samples["type"] == "tumor", "sample"].tolist()
normals = samples.loc[samples["type"] == "normal", "sample"].tolist()

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")

# Define analysis steps based on config flags
# Default flags if not specified
steps = {
    "bam_refiner": config.get("steps", {}).get("bam_refiner", True),
    "methylation": config.get("steps", {}).get("methylation", False),
    "copynumber": config.get("steps", {}).get("copynumber", False),
    "nanomonsv_parse": config.get("steps", {}).get("nanomonsv_parse", False),
    "nanomonsv_get": config.get("steps", {}).get("nanomonsv_get", False),
    "nanomonsv_postprocess": config.get("steps", {}).get("nanomonsv_postprocess", False),
    "nanomonsv_insert_classify": config.get("steps", {}).get("nanomonsv_insert_classify", False),
    "nanomonsv_connect": config.get("steps", {}).get("nanomonsv_connect", False),
    "nanomonsv_merge": config.get("steps", {}).get("nanomonsv_merge", False),
    "clairs": config.get("steps", {}).get("clairs", False),
    "deepsomatic": config.get("steps", {}).get("deepsomatic", False),
    "clairs_postprocess": config.get("steps", {}).get("clairs_postprocess", False),
    "deepsomatic_postprocess": config.get("steps", {}).get("deepsomatic_postprocess", False),
}

# Helper functions for getting sample metadata
def get_sample_hifi(sample):
    return samples.loc[sample, "hifi"]

def get_sample_ont(sample):
    return samples.loc[sample, "ont"]

def get_sample_assembly_hap1(sample):
    return samples.loc[sample, "assembly_hap1"]

def get_sample_assembly_hap2(sample):
    return samples.loc[sample, "assembly_hap2"]

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
        raise ValueError(f"Cannot determine paired normal for tumor {tumor}. "
                        "Either provide only one normal sample or add 'paired_normal' column to sample sheet.")

# Default resource configurations for each tool
# These can be overridden in config.yaml under 'resources' key
def get_threads(tool, default=8):
    """Get thread count for a tool from config or use default."""
    return config.get("resources", {}).get(tool, {}).get("threads", default)

def get_mem_mb(tool, default=32000):
    """Get memory (MB) for a tool from config or use default."""
    return config.get("resources", {}).get(tool, {}).get("mem_mb", default)
