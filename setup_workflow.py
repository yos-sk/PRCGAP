#!/usr/bin/env python3
"""
Setup script for PRCGAP Snakemake workflow configuration.
This script creates config.yaml from command-line arguments.
"""

import argparse
import sys
import yaml
from pathlib import Path


def create_config(args):
    """Create configuration dictionary from arguments."""
    config = {
        "samplesheet": args.samplesheet,
        "reference": args.reference,
        "singularity_images": {
            "bam_refiner": args.bam_refiner_image,
            "methylation": args.methylation_image,
            "copynumber": args.copynumber_image,
            "nanomonsv": args.nanomonsv_image,
            "nanomonsv_postprocess": args.nanomonsv_postprocess_image,
            "clairs": args.clairs_image,
            "deepsomatic": args.deepsomatic_image,
            "point_mutation_postprocess": args.point_mutation_postprocess_image,
        },
        "resources": {
            "bam_refiner": {"threads": args.threads, "mem_mb": args.mem_mb},
            "methylation": {"threads": 8, "mem_mb": 32000},
            "copynumber": {"threads": 8, "mem_mb": 32000},
            "nanomonsv_parse": {"threads": 4, "mem_mb": 16000},
            "nanomonsv_get": {"threads": 8, "mem_mb": 32000},
            "nanomonsv_postprocess": {"threads": 4, "mem_mb": 16000},
            "nanomonsv_insert_classify": {"threads": 4, "mem_mb": 16000},
            "nanomonsv_connect": {"threads": 4, "mem_mb": 16000},
            "nanomonsv_merge": {"threads": 2, "mem_mb": 8000},
            "clairs": {"threads": 16, "mem_mb": 64000},
            "deepsomatic": {"threads": 16, "mem_mb": 64000},
            "clairs_postprocess": {"threads": 4, "mem_mb": 16000},
            "clairs_postprocess_realign": {"threads": 16, "mem_mb": 32000},
            "clairs_postprocess_pileup": {"threads": 16, "mem_mb": 32000},
            "deepsomatic_postprocess": {"threads": 4, "mem_mb": 16000},
            "deepsomatic_postprocess_realign": {"threads": 16, "mem_mb": 32000},
            "deepsomatic_postprocess_pileup": {"threads": 16, "mem_mb": 32000},
        },
        "steps": {
            "bam_refiner": args.enable_bam_refiner,
            "methylation": args.enable_methylation,
            "copynumber": args.enable_copynumber,
            "nanomonsv_parse": args.enable_nanomonsv_parse,
            "nanomonsv_get": args.enable_nanomonsv_get,
            "nanomonsv_postprocess": args.enable_nanomonsv_postprocess,
            "nanomonsv_insert_classify": args.enable_nanomonsv_insert_classify,
            "nanomonsv_connect": args.enable_nanomonsv_connect,
            "nanomonsv_merge": args.enable_nanomonsv_merge,
            "clairs": args.enable_clairs,
            "deepsomatic": args.enable_deepsomatic,
            "clairs_postprocess": args.enable_clairs_postprocess,
            "deepsomatic_postprocess": args.enable_deepsomatic_postprocess,
        },
        "hap1_satellite": args.hap1_satellite or "",
        "hap2_satellite": args.hap2_satellite or "",
        "sex": args.sex,
        "simple_repeat": args.simple_repeat or "",
        "gtf_file": args.gtf_file or "",
        "line1_bed": args.line1_bed or "",
    }
    return config


def main():
    """Main function to parse arguments and create config."""
    parser = argparse.ArgumentParser(
        description="Setup PRCGAP Snakemake workflow configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic setup with BAM refiner only
  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --reference ref.fa \\
    --bam-refiner-image bam_refiner.sif

  # Full pipeline setup (using --full-pipeline shortcut)
  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --reference ref.fa \\
    --bam-refiner-image bam_refiner.sif \\
    --methylation-image methylation.sif \\
    --copynumber-image copynumber.sif \\
    --nanomonsv-image nanomonsv.sif \\
    --nanomonsv-postprocess-image nanomonsv_postprocess.sif \\
    --clairs-image clairs.sif \\
    --deepsomatic-image deepsomatic.sif \\
    --point-mutation-postprocess-image mutation_postprocess.sif \\
    --full-pipeline \\
    --hap1-satellite hap1_sat.bed \\
    --hap2-satellite hap2_sat.bed \\
    --sex female \\
    --gtf-file annotation.gtf \\
    --simple-repeat simple_repeat.bed \\
    --line1-bed line1.bed \\
    --output config/config.yaml

  # With copy number analysis
  python setup_workflow.py \\
    --samplesheet samples.tsv \\
    --reference ref.fa \\
    --bam-refiner-image bam_refiner.sif \\
    --enable-copynumber \\
    --hap1-satellite hap1_sat.bed \\
    --hap2-satellite hap2_sat.bed \\
    --sex female \\
    --output config/config.yaml
        """
    )

    # Required arguments
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="Path to sample sheet TSV file (required)"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference genome FASTA file (required)"
    )

    # Singularity images
    parser.add_argument(
        "--bam-refiner-image",
        required=True,
        help="Path to BAM refiner Singularity image (required)"
    )
    parser.add_argument(
        "--methylation-image",
        default="",
        help="Path to methylation calling Singularity image"
    )
    parser.add_argument(
        "--copynumber-image",
        default="",
        help="Path to copy number calling Singularity image"
    )
    parser.add_argument(
        "--nanomonsv-image",
        default="",
        help="Path to NanoMonSV Singularity image (for parse, get, and other steps)"
    )
    parser.add_argument(
        "--nanomonsv-postprocess-image",
        default="",
        help="Path to NanoMonSV postprocess Singularity image"
    )
    parser.add_argument(
        "--clairs-image",
        default="",
        help="Path to ClairS Singularity image"
    )
    parser.add_argument(
        "--deepsomatic-image",
        default="",
        help="Path to DeepSomatic Singularity image"
    )
    parser.add_argument(
        "--point-mutation-image",
        default="",
        help="Path to point mutation (DeepSomatic/ClairS) Singularity image (deprecated, use --clairs-image and --deepsomatic-image)"
    )

    # Full pipeline shortcut
    parser.add_argument(
        "--full-pipeline",
        action="store_true",
        default=False,
        help="Enable all analysis steps (equivalent to enabling all --enable-* options)"
    )

    # Analysis steps (enable/disable)
    parser.add_argument(
        "--enable-bam-refiner",
        action="store_true",
        default=True,
        help="Enable BAM refiner (default: True, this is typically the base step)"
    )
    parser.add_argument(
        "--disable-bam-refiner",
        action="store_false",
        dest="enable_bam_refiner",
        help="Disable BAM refiner"
    )
    parser.add_argument(
        "--enable-methylation",
        action="store_true",
        default=False,
        help="Enable methylation analysis"
    )
    parser.add_argument(
        "--enable-copynumber",
        action="store_true",
        default=False,
        help="Enable copy number analysis"
    )
    parser.add_argument(
        "--enable-nanomonsv-parse",
        action="store_true",
        default=False,
        help="Enable NanoMonSV parse step"
    )
    parser.add_argument(
        "--enable-nanomonsv-get",
        action="store_true",
        default=False,
        help="Enable NanoMonSV get step (requires parse step)"
    )
    parser.add_argument(
        "--enable-nanomonsv-postprocess",
        action="store_true",
        default=False,
        help="Enable NanoMonSV postprocess step (requires get step)"
    )
    parser.add_argument(
        "--enable-nanomonsv-insert-classify",
        action="store_true",
        default=False,
        help="Enable NanoMonSV insert classification (requires get step)"
    )
    parser.add_argument(
        "--enable-nanomonsv-connect",
        action="store_true",
        default=False,
        help="Enable NanoMonSV connect step (requires postprocess step)"
    )
    parser.add_argument(
        "--enable-nanomonsv-merge",
        action="store_true",
        default=False,
        help="Enable NanoMonSV merge step (requires postprocess step)"
    )
    parser.add_argument(
        "--enable-clairs",
        action="store_true",
        default=False,
        help="Enable ClairS variant calling"
    )
    parser.add_argument(
        "--enable-deepsomatic",
        action="store_true",
        default=False,
        help="Enable DeepSomatic variant calling"
    )
    parser.add_argument(
        "--enable-clairs-postprocess",
        action="store_true",
        default=False,
        help="Enable ClairS postprocessing (requires clairs step)"
    )
    parser.add_argument(
        "--enable-deepsomatic-postprocess",
        action="store_true",
        default=False,
        help="Enable DeepSomatic postprocessing (requires deepsomatic step)"
    )

    # Resource parameters (defaults apply to bam_refiner, other tools use their own defaults)
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Default number of threads for bam_refiner (default: 16). Other tools use their own defaults."
    )
    parser.add_argument(
        "--mem-mb",
        type=int,
        default=64000,
        help="Default memory in MB for bam_refiner (default: 64000). Other tools use their own defaults."
    )
    parser.add_argument(
        "--sex",
        choices=["female", "male"],
        default="female",
        help="Sex of the sample for copy number analysis (default: female)"
    )
    parser.add_argument(
        "--hap1-satellite",
        default="",
        help="Path to satellite regions for haplotype 1"
    )
    parser.add_argument(
        "--hap2-satellite",
        default="",
        help="Path to satellite regions for haplotype 2"
    )
    parser.add_argument(
        "--simple-repeat",
        default="",
        help="Path to simple repeat regions file"
    )
    parser.add_argument(
        "--gtf-file",
        default="",
        help="Path to GTF annotation file"
    )
    parser.add_argument(
        "--line1-bed",
        default="",
        help="Path to LINE1 element BED file"
    )
    parser.add_argument(
        "--point-mutation-postprocess-image",
        default="",
        help="Path to point mutation postprocessing Singularity image"
    )

    # Output
    parser.add_argument(
        "--output",
        "-o",
        default="config/config.yaml",
        help="Output path for config.yaml (default: config/config.yaml)"
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        default=False,
        help="Overwrite existing config file"
    )

    args = parser.parse_args()

    # If --full-pipeline is specified, enable all analysis steps
    if args.full_pipeline:
        args.enable_bam_refiner = True
        args.enable_methylation = True
        args.enable_copynumber = True
        args.enable_nanomonsv_parse = True
        args.enable_nanomonsv_get = True
        args.enable_nanomonsv_postprocess = True
        args.enable_nanomonsv_insert_classify = True
        args.enable_nanomonsv_connect = True
        args.enable_nanomonsv_merge = True
        args.enable_clairs = True
        args.enable_deepsomatic = True
        args.enable_clairs_postprocess = True
        args.enable_deepsomatic_postprocess = True

    # Check if output file exists
    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        print(f"Error: {args.output} already exists. Use --force to overwrite.")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create configuration
    config = create_config(args)

    # Write configuration to YAML
    try:
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        print(f"✓ Configuration written to {output_path}")
        print("\nConfiguration summary:")
        print(f"  Samplesheet: {args.samplesheet}")
        print(f"  Reference: {args.reference}")
        print("\nEnabled steps:")
        for step, enabled in config["steps"].items():
            if enabled:
                print(f"  ✓ {step}")
        print("\nResource defaults (can be overridden in config.yaml):")
        print(f"  bam_refiner: threads={args.threads}, mem_mb={args.mem_mb}")
        print("  (See config.yaml for all tool resource settings)")
        print("\nSingularity images:")
        for img, path in config["singularity_images"].items():
            if path:
                print(f"  {img}: {path}")
    except Exception as e:
        print(f"Error writing configuration: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
