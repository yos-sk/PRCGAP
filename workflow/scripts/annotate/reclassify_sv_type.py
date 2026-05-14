#!/usr/bin/env python3

"""
Reclassify SV types for inter-contig events using CHM13-normalized
breakpoint directions.

Reclassification is applied only when the two breakpoints fall on
different personalized contigs (Chr_1 != Chr_2). Intra-contig events
(Chr_1 == Chr_2) retain the original SV_Type.

For inter-contig events:
    - Different CHM13 chromosomes                       -> Translocation
    - Same CHM13 chromosome, different haplotypes
      (inter-haplotype)                                 -> Translocation
    - Same CHM13 chromosome, same haplotype             -> reclassify
      using CHM13-normalized directions following the nanomonsv rule
      set (https://github.com/friend1ws/nanomonsv/blob/master/misc/sv_type.py):
          Dir_1 == Dir_2                 -> Inversion
          Dir_1 == '-' and Dir_2 == '+'  -> Duplication
          Dir_1 == '+' and Dir_2 == '-'  -> Deletion / Insertion
                                           (Insertion if the inserted
                                           sequence is longer than the
                                           approximate CHM13 span)

Haplotype membership is inferred from which reference table file the
contig came from; supply one --ref-table per haplotype.
"""

import argparse
import csv
import sys


NEW_COLUMNS = [
    "CHM13_Chr_1",
    "CHM13_Chr_2",
    "Contig_Strand_1",
    "Contig_Strand_2",
    "Dir_1_chm13",
    "Dir_2_chm13",
    "SV_Type_chm13",
]


def load_ref_table(paths: list) -> dict:
    """Load reference table(s) (output of make_reference_table.py).

    Columns (TSV, no header):
        0: personalized contig (e.g., haplotype1-0000023)
        1: contig start
        2: contig end
        3: strand ('+' or '-'; orientation of the personalized contig relative to CHM13)
        4: CHM13 chromosome (e.g., chr1)
        5: CHM13 start
        6: CHM13 end

    Multiple table files can be provided (e.g., one per haplotype); their
    entries are merged into a single lookup. The zero-based index of the
    source file is retained so that inter-haplotype SVs can be detected.

    Returns:
        dict mapping personalized contig ->
            (chm13_chr, strand, chm13_start, chm13_end, hap_idx).
        If a contig appears multiple times, the entry with the longest
        CHM13 span is retained.
    """
    table = {}
    for idx, path in enumerate(paths):
        with open(path, "r") as f:
            for line in f:
                items = line.rstrip("\n").split("\t")
                if len(items) < 7:
                    continue
                contig = items[0]
                strand = items[3]
                chm13_chr = items[4]
                chm13_start = int(items[5])
                chm13_end = int(items[6])

                existing = table.get(contig)
                new_span = chm13_end - chm13_start
                if existing is None or new_span > (existing[3] - existing[2]):
                    table[contig] = (chm13_chr, strand, chm13_start, chm13_end, idx)
    return table


def flip_dir(d: str) -> str:
    return "-" if d == "+" else "+"


def normalize_dir(d: str, strand: str) -> str:
    return d if strand == "+" else flip_dir(d)


def classify_sv(row: dict, ref_table: dict) -> dict:
    chr_1, chr_2 = row["Chr_1"], row["Chr_2"]
    dir_1, dir_2 = row["Dir_1"], row["Dir_2"]
    inseq = "" if row["Inserted_Seq"] == "---" else row["Inserted_Seq"]
    original_sv_type = row["SV_Type"]

    result = {k: "NA" for k in NEW_COLUMNS}
    result["SV_Type_chm13"] = original_sv_type

    entry_1 = ref_table.get(chr_1)
    entry_2 = ref_table.get(chr_2)

    if entry_1 is not None:
        result["CHM13_Chr_1"] = entry_1[0]
        result["Contig_Strand_1"] = entry_1[1]
        result["Dir_1_chm13"] = normalize_dir(dir_1, entry_1[1])
    if entry_2 is not None:
        result["CHM13_Chr_2"] = entry_2[0]
        result["Contig_Strand_2"] = entry_2[1]
        result["Dir_2_chm13"] = normalize_dir(dir_2, entry_2[1])

    # Intra-contig events keep the original SV_Type.
    if chr_1 == chr_2:
        return result

    # Inter-contig events: reclassify.
    # If either contig is missing from the ref table, we cannot determine
    # the CHM13 chromosome, so fall back to the original SV_Type.
    if entry_1 is None or entry_2 is None:
        return result

    chm13_chr_1, _, ts_1, _, hap_idx_1 = entry_1
    chm13_chr_2, _, ts_2, _, hap_idx_2 = entry_2

    if chm13_chr_1 != chm13_chr_2:
        result["SV_Type_chm13"] = "Translocation"
        return result

    if hap_idx_1 != hap_idx_2:
        # Inter-haplotype SV mapping to the same CHM13 chromosome.
        result["SV_Type_chm13"] = "Translocation"
        return result

    # Same haplotype, different contigs, same CHM13 chromosome.
    # Sort the two breakpoints by CHM13 coordinate using each contig's
    # chm13_start as a proxy. This is an approximation: strictly, the
    # order should be determined from each breakpoint's actual CHM13
    # position (contig_start + (pos - contig_start) / chm13_end - (pos -
    # contig_start), depending on strand). The proxy fails only when the
    # two contigs' CHM13 ranges overlap AND both breakpoints fall inside
    # that overlap, which is rare in practice; switch to per-breakpoint
    # coordinate conversion if stricter behavior is needed.
    dir_1_chm13 = result["Dir_1_chm13"]
    dir_2_chm13 = result["Dir_2_chm13"]
    span = abs(ts_2 - ts_1) + 1
    if ts_1 > ts_2:
        dir_1_chm13, dir_2_chm13 = dir_2_chm13, dir_1_chm13

    if dir_1_chm13 == dir_2_chm13:
        sv_type = "Inversion"
    elif dir_1_chm13 == "-" and dir_2_chm13 == "+":
        sv_type = "Duplication"
    elif dir_1_chm13 == "+" and dir_2_chm13 == "-":
        if len(inseq) > span:
            sv_type = "Insertion"
        else:
            sv_type = "Deletion"
    else:
        sv_type = "Unresolved"

    result["SV_Type_chm13"] = sv_type
    return result


def reclassify_main(args) -> None:
    ref_table = load_ref_table(args.ref_table)

    with open(args.input, "r") as hin, open(args.output, "w") as hout:
        dreader = csv.DictReader(hin, delimiter="\t")
        header = list(dreader.fieldnames) + NEW_COLUMNS
        print("\t".join(header), file=hout)
        for row in dreader:
            annot = classify_sv(row, ref_table)
            out_values = [row[k] for k in dreader.fieldnames] + [annot[k] for k in NEW_COLUMNS]
            print("\t".join(out_values), file=hout)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog="reclassify_sv_type.py",
        description=(
            "Classify SV types using CHM13-normalized breakpoint directions "
            "(nanomonsv rule set applied on CHM13-normalized coordinates)."
        ),
    )
    parser.add_argument(
        "--input", "-i", type=str, required=True,
        help="Input nanomonsv result TSV with Chr_1, Pos_1, Dir_1, Chr_2, Pos_2, Dir_2, Inserted_Seq columns.",
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True,
        help="Output TSV file.",
    )
    parser.add_argument(
        "--ref-table", "-r", type=str, nargs="+", required=True,
        help="Reference table(s) output by make_reference_table.py "
             "(7-column TSV: contig, contig_start, contig_end, strand, chm13_chr, chm13_start, chm13_end). "
             "Provide one path per haplotype for diploid references.",
    )
    parser.set_defaults(func=reclassify_main)

    return parser


def main():
    parser = arg_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
