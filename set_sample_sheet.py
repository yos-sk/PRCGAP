#!/usr/bin/env python3
"""Generate a PRCGAP sample sheet TSV for one tumor-normal pair.

A PRCGAP run is built around a single tumor-normal pair sharing one
personalized (diploid) assembly: the per-assembly annotation resources
(satellite, simple_repeat, gtf/gff, chains, repeat-masker, segdup, censat,
misassembly, ...) are GLOBAL in config -- one set per run -- so a run is tied
to a single normal assembly, and tumor<->normal pairing is resolved
automatically. Distinct pairs (or distinct assemblies) must be run separately
(own config + sample sheet + output dir).

This script does NOT auto-discover files. You name the tumor and normal samples
and point each at its ONT and HiFi data with dedicated options; the shared
assembly is given once. The script validates the entries, absolutises the data
paths, and writes the canonical sample sheet consumed by the workflow
(workflow/rules/commons.smk).

Columns written: sample, type, ont, hifi, assembly_hap1, assembly_hap2.

Examples
--------
  # One tumor + one normal sharing one assembly (the common case):
  python3 set_sample_sheet.py \
      --tumor  HG008T --tumor-ont  reads/HG008T.ont.bam  --tumor-hifi  reads/HG008T.hifi.bam \
      --normal HG008N --normal-ont reads/HG008N.ont.bam  --normal-hifi reads/HG008N.hifi.bam \
      --assembly-hap1 asm/hap1.fa --assembly-hap2 asm/hap2.fa \
      -o config/samples.tsv

  # Multiple files per data type (repeat the flag, or comma-separate):
  python3 set_sample_sheet.py \
      --tumor  T --tumor-ont  t.r1.bam --tumor-ont  t.r2.bam --tumor-hifi  t.hifi.bam \
      --normal N --normal-ont n.bam                          --normal-hifi n.hifi.bam \
      --assembly-hap1 asm/hap1.fa --assembly-hap2 asm/hap2.fa \
      -o config/samples.tsv

  # Validate/format an existing draft TSV (absolutise paths, run the checks):
  python3 set_sample_sheet.py --input draft.tsv -o config/samples.tsv
"""

import argparse
import csv
import os
import sys

COLUMNS = ["sample", "type", "ont", "hifi", "assembly_hap1", "assembly_hap2"]
SEQTYPE_FIELDS = ("ont", "hifi")
PATH_FIELDS = ("ont", "hifi", "assembly_hap1", "assembly_hap2")


def _abs(path):
    """Absolutise a path WITHOUT following symlinks. Empty/None passes through."""
    if not path:
        return path
    return os.path.abspath(os.path.expanduser(path))


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
        if not ont:
            raise SampleSheetError(f"{sample_type} sample {name}: no --{sample_type}-ont given")
        if not hifi:
            raise SampleSheetError(f"{sample_type} sample {name}: no --{sample_type}-hifi given")
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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    pair = parser.add_argument_group(
        "tumor-normal pair (option mode)",
        "Name the two samples and point each at its ONT and HiFi data. "
        "ONT/HiFi options are repeatable and also accept comma-separated lists.")
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
                      help="shared assembly hap1 fasta (used by both samples)")
    pair.add_argument("--assembly-hap2", metavar="FASTA",
                      help="shared assembly hap2 fasta (used by both samples)")

    parser.add_argument(
        "--input", "-i", default=None,
        help="draft sample sheet TSV to validate/format (same columns as the "
             "output) instead of building from --tumor/--normal options.")
    parser.add_argument(
        "--output", "-o", default="config/samples.tsv",
        help="output sample sheet path (default: config/samples.tsv)")
    parser.add_argument(
        "--no-check-exists", dest="check_exists", action="store_false", default=True,
        help="do not check that data files exist (e.g. when they live on a cluster)")
    parser.add_argument(
        "--no-absolutize", dest="absolutize", action="store_false", default=True,
        help="keep paths as given instead of absolutising them")
    parser.add_argument(
        "--force", "-f", action="store_true", default=False,
        help="overwrite the output if it already exists")

    args = parser.parse_args()

    option_mode = any((args.tumor, args.normal, args.assembly_hap1, args.assembly_hap2,
                       args.tumor_ont, args.tumor_hifi, args.normal_ont, args.normal_hifi))
    if not option_mode and not args.input:
        parser.error("provide the --tumor/--normal options or an --input TSV")

    if os.path.exists(args.output) and not args.force:
        parser.error(f"{args.output} already exists. Use --force to overwrite.")

    try:
        rows = []
        if args.input:
            rows.extend(read_input_tsv(args.input))
        if option_mode:
            rows.extend(build_rows_from_pair(args))
        rows = normalize_and_validate(
            rows, check_exists=args.check_exists, absolutize=args.absolutize)
    except (SampleSheetError, FileNotFoundError) as exc:
        parser.error(str(exc))

    write_sample_sheet(rows, args.output)

    n_tumor = sum(1 for r in rows if r["type"] == "tumor")
    n_normal = sum(1 for r in rows if r["type"] == "normal")
    print(f"✓ Wrote {args.output}: {len(rows)} samples "
          f"({n_tumor} tumor, {n_normal} normal).")


if __name__ == "__main__":
    main()
