#!/usr/bin/env python3
"""Annotate a prepped mutation table with GRCh38/CHM13 lifted coordinates.

This script replaces the *coordinate annotation* portion of the reference
`compare_current_reference.py` scripts (point_mutation_analysis/scripts/
annotation/{snv,indel}/compare_current_reference.py). The reference
scripts bundle two responsibilities:

  (a) add GRCh38/CHM13 lifted-coordinate columns to each variant, and
  (b) flag whether the same variant was called against GRCh38/CHM13
      reference VCFs.

(b) is the comparison step that proposal step 3 explicitly defers, so
this script keeps (a) only. The reference's `GRCh38_flag`/`chm13_flag`
columns (which would hold the result of (b)) are omitted entirely rather
than emitted as `-` placeholders; re-add them here (and re-add the two
index positions in check_homozygous_{snv,indel}.py) if (b) is implemented.

Two modes:
  --mode snv    Reads `coordconv` output BEDs (one for GRCh38, one for
                CHM13). The original BED fed to coordconv is the 14-col
                haplotyped.bed filtered for SNVs. coordconv preserves
                the input columns past the inserted dst/status columns,
                so the GRCh38/CHM13 outputs each carry every input
                column too. We use those output rows as the source of
                record.

  --mode indel  Reads `transanno liftvcf` output VCFs (gzipped). The
                INFO field's first five `;`-separated key=value pairs
                are `haplotype=...;contig=...;pos=...;ref=...;alt=...`
                (set by bed2vcf.py before lifting) — we use that
                composite key to join GRCh38 ↔ CHM13 lifts back onto
                the original BED rows.

Output header is fixed for each mode and aligns with what the reference
emits at this point in the pipeline; downstream `add_annotation_mut.py`
reads with csv.DictReader so column ORDER is irrelevant — column NAMES
are what matter.
"""

import argparse
import gzip
import sys


# ----------------------------- common helpers -----------------------------

# Haplotyped.bed (14 columns, no header) -> downstream-friendly field names.
# Index → column name; matches the rest of the workflow's expectations.
BED_COLUMNS = [
    "Contig", "Start", "End", "Ref", "Alt",
    "Score", "Filter",
    "VAF", "Depth", "Normal_VAF", "Normal_depth",
    "ID", "Haplotype", "Kmer_ratio",
]

# Columns we drop from the final output (matches reference: Score/Filter
# are not carried forward past compare_current_reference).
EMIT_COLUMNS_BASE = [
    "Contig", "Start", "End", "Ref", "Alt",
]
EMIT_COLUMNS_TAIL = [
    "VAF", "Depth", "Normal_VAF", "Normal_depth",
    "ID", "Haplotype", "Kmer_ratio",
]


def _parse_bed_row(items):
    row = dict(zip(BED_COLUMNS, items[: len(BED_COLUMNS)]))
    return row


# ----------------------------- SNV mode -----------------------------------
#
# coordconv output (one row per input BED row) has the shape:
#   src_contig src_pos dst_contig dst_pos status <input cols 3..N>
# When fed a 14-col haplotyped.bed, the trailing `<input cols 3..N>` is
# Ref, Alt, Score, Filter, VAF, Depth, Normal_VAF, Normal_depth, ID,
# Haplotype, Kmer_ratio. We key by ID since coordconv pass-throughs it
# unchanged, allowing GRCh38 and CHM13 outputs to be joined back.

def _read_coordconv(path):
    """Return {ID: (dst_contig, dst_pos, status)}.

    coordconv on a 14-column haplotyped.bed produces:
      Match row (16 cols): src_contig src_pos dst_contig dst_pos status
                           Ref Alt Score Filter VAF Depth Normal_VAF
                           Normal_depth ID Haplotype Kmer_ratio
      Gap row   (17 cols): like Match but with an extra dst_end inserted
                           between dst_pos and status.
    The pass-through block (Ref..Kmer_ratio) is always the last 11
    columns; status sits immediately before it (items[-12]); ID is the
    third-from-last (items[-3]); dst_contig/dst_pos are always at
    items[2]/items[3]. Empty `path` (chain not configured) returns {}.
    """
    out = {}
    if not path:
        return out
    with open(path, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            if len(items) < 16:
                continue
            dst_contig = items[2]
            dst_pos = items[3]
            status = items[-12]
            mid = items[-3]
            if status in ("Match", "Gap"):
                out[mid] = (dst_contig, dst_pos, status)
            else:
                out[mid] = ("-", "-", status)
    return out


def run_snv(args):
    grch38 = _read_coordconv(args.grch38)
    chm13 = _read_coordconv(args.chm13)

    header = (
        EMIT_COLUMNS_BASE
        + [
            "GRCh38_contig", "GRCh38_pos", "GRCh38_status",
            "chm13_contig", "chm13_pos", "chm13_status",
        ]
        + EMIT_COLUMNS_TAIL
    )

    with open(args.input_bed, "r") as f, open(args.output, "w") as w:
        w.write("\t".join(header) + "\n")
        for line in f:
            items = line.rstrip("\n").split("\t")
            if len(items) < len(BED_COLUMNS):
                continue
            row = _parse_bed_row(items)
            mid = row["ID"]
            g = grch38.get(mid, ("-", "-", "-"))
            c = chm13.get(mid, ("-", "-", "-"))
            out = (
                [row[k] for k in EMIT_COLUMNS_BASE]
                + list(g) + list(c)
                + [row[k] for k in EMIT_COLUMNS_TAIL]
            )
            w.write("\t".join(out) + "\n")


# ----------------------------- INDEL mode ---------------------------------
#
# transanno liftvcf preserves the INFO=`haplotype=...;contig=...;pos=...;
# ref=...;alt=...` set by bed2vcf.py and prepends additional info on
# allele transformations when applicable. The key for joining is the
# first five `;` fields of the INFO column.

def _read_liftvcf(path):
    """Return {info_key: (dst_contig, dst_pos, dst_ref, dst_alt)}.

    Empty path returns {}.
    """
    out = {}
    if not path:
        return out
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            items = line.rstrip("\n").split("\t")
            if len(items) < 8:
                continue
            chrom, pos, _, ref, alt = items[0], items[1], items[2], items[3], items[4]
            info_fields = items[7].split(";")
            # Reference uses items[7].split(";")[:5]; that assumes the
            # five tagged fields come first. transanno may prepend
            # `ReverseComplementedAlleles=...`; pick the five tagged
            # fields by name instead so it works either way.
            wanted = {"haplotype", "contig", "pos", "ref", "alt"}
            tagged = [kv for kv in info_fields if kv.split("=", 1)[0] in wanted]
            if len(tagged) < 5:
                continue
            # Re-order tagged so the key is canonical.
            d = dict(kv.split("=", 1) for kv in tagged)
            key = (
                f"haplotype={d['haplotype']};contig={d['contig']};"
                f"pos={d['pos']};ref={d['ref']};alt={d['alt']}"
            )
            out[key] = (chrom, pos, ref, alt)
    return out


def run_indel(args):
    grch38 = _read_liftvcf(args.grch38)
    chm13 = _read_liftvcf(args.chm13)

    header = (
        EMIT_COLUMNS_BASE
        + [
            "GRCh38_contig", "GRCh38_pos", "GRCh38_ref", "GRCh38_alt", "GRCh38_status",
            "chm13_contig", "chm13_pos", "chm13_ref", "chm13_alt", "chm13_status",
        ]
        + EMIT_COLUMNS_TAIL
    )

    with open(args.input_bed, "r") as f, open(args.output, "w") as w:
        w.write("\t".join(header) + "\n")
        for line in f:
            items = line.rstrip("\n").split("\t")
            if len(items) < len(BED_COLUMNS):
                continue
            row = _parse_bed_row(items)
            key = (
                f"haplotype={row['Haplotype']};contig={row['Contig']};"
                f"pos={row['End']};ref={row['Ref']};alt={row['Alt']}"
            )
            if key in grch38:
                g_contig, g_pos, g_ref, g_alt = grch38[key]
                g_status = "Match"
            else:
                g_contig = g_pos = g_ref = g_alt = "-"
                g_status = "Rejected"

            if key in chm13:
                c_contig, c_pos, c_ref, c_alt = chm13[key]
                c_status = "Match"
            else:
                c_contig = c_pos = c_ref = c_alt = "-"
                c_status = "Rejected"

            out = (
                [row[k] for k in EMIT_COLUMNS_BASE]
                + [g_contig, g_pos, g_ref, g_alt, g_status]
                + [c_contig, c_pos, c_ref, c_alt, c_status]
                + [row[k] for k in EMIT_COLUMNS_TAIL]
            )
            w.write("\t".join(out) + "\n")


# ----------------------------- entrypoint ---------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="add_lift_coords.py",
        description=(
            "Add GRCh38/CHM13 lifted-coordinate columns to a prepped "
            "mutation table. Mirrors the coordinate-annotation portion "
            "of compare_current_reference.py (comparison part skipped)."
        ),
    )
    parser.add_argument("--mode", required=True, choices=["snv", "indel"])
    parser.add_argument("--input_bed", "-i", required=True,
                        help="14-column haplotyped.bed filtered for SNV or INDEL")
    parser.add_argument("--output", "-o", required=True,
                        help="Output TSV with header")
    parser.add_argument("--grch38", "-g", default="",
                        help="GRCh38 lift output (coordconv BED for snv, transanno VCF for indel)")
    parser.add_argument("--chm13", "-c", default="",
                        help="CHM13 lift output (coordconv BED for snv, transanno VCF for indel)")

    args = parser.parse_args()
    if args.mode == "snv":
        run_snv(args)
    else:
        run_indel(args)


if __name__ == "__main__":
    main()
