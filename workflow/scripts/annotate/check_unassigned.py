#!/usr/bin/env python3
"""Drop unassigned duplicates of a variant that a haplotype-assigned record
already covers at the same reference coordinate.

Both haplotype contigs can carry the same somatic variant, and the two records
lift to one reference position. When only one of them got a haplotype, the other
is the same event seen without phasing information and adds nothing; when
neither did, one of the two is enough. Records that both carry a haplotype are
both kept: that is the haplotype-resolved output this pipeline exists to
produce, and collapsing them would discard the assignment.

Rows are matched on the reference contig and position plus the alleles, and are
emitted in input order with every field untouched -- the only effect is that
some rows are absent.

  --mode snv    alleles compared as Ref/Alt, or their reverse complement (the
                two contigs can align to the reference in opposite
                orientations). This mode carries MNVs as well as single-base
                substitutions, so the alleles are treated as strings.
  --mode indel  alleles compared as the lifted <ref>_ref / <ref>_alt, which
                transanno resolves per reference

Columns are located by header name, so the annotated table's column order is
not baked in here.
"""

import argparse
import sys

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def _revcomp(allele):
    """Reverse complement, leaving anything that is not a base alone. An allele
    is a string, not a base: `snv` mode covers MNVs (prep_mut.sh sends every
    equal-length ref/alt pair there), so reversing matters and a per-base lookup
    that assumed length 1 would fail on them."""
    return "".join(COMPLEMENT.get(b.upper(), b) for b in reversed(allele))


def _allele_keys(row, mode, ref_prefix):
    """Allele forms under which two rows count as the same variant."""
    if mode == "indel":
        return {(row[f"{ref_prefix}_ref"], row[f"{ref_prefix}_alt"])}
    ref, alt = row["Ref"], row["Alt"]
    # The two haplotype contigs can align to the reference in opposite
    # orientations, so the same event appears reverse-complemented on one of
    # them. For a single base this is just the complement.
    return {(ref, alt), (_revcomp(ref), _revcomp(alt))}


def _is_assigned(row):
    return row["Haplotype"].startswith("haplotype")


def _float(value):
    try:
        return float(value)
    except ValueError:
        return 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--mode", choices=["snv", "indel"], required=True)
    parser.add_argument("--reference", choices=["GRCh38", "chm13"], required=True)
    args = parser.parse_args()

    with open(args.input) as hin:
        header = hin.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in hin]

    contig_col = f"{args.reference}_contig"
    pos_col = f"{args.reference}_pos"
    needed = [contig_col, pos_col, "Ref", "Alt", "Haplotype", "VAF", "Depth"]
    if args.mode == "indel":
        needed += [f"{args.reference}_ref", f"{args.reference}_alt"]
    missing = [c for c in needed if c not in header]
    if missing:
        sys.exit(f"{args.input}: missing column(s) {', '.join(missing)}")

    # One bucket per (reference position, allele form). A row joins every bucket
    # its allele forms produce, so a complement match lands in the same bucket
    # as its counterpart.
    buckets = {}
    for i, row in enumerate(rows):
        pos = row[pos_col]
        # An unlifted row has no reference coordinate to be a duplicate at.
        if not pos or "-" in pos or not row[contig_col]:
            continue
        for allele in _allele_keys(row, args.mode, args.reference):
            buckets.setdefault((row[contig_col], pos, allele), []).append(i)

    dropped = set()
    for members in buckets.values():
        members = [i for i in members if i not in dropped]
        if len(members) < 2:
            continue
        assigned = [i for i in members if _is_assigned(rows[i])]
        if assigned:
            # Keep every haplotype-assigned record; the unassigned ones are the
            # same event without phasing.
            dropped.update(i for i in members if i not in assigned)
        else:
            # Nothing to prefer on haplotype: keep the best-supported one.
            keep = max(members,
                       key=lambda i: (_float(rows[i]["VAF"]),
                                      _float(rows[i]["Depth"]), -i))
            dropped.update(i for i in members if i != keep)

    with open(args.output, "w") as hout:
        print("\t".join(header), file=hout)
        for i, row in enumerate(rows):
            if i in dropped:
                continue
            print("\t".join(row.get(c, "") for c in header), file=hout)

    print(f"[check_unassigned] {args.mode} / {args.reference}: "
          f"{len(rows)} rows in, {len(rows) - len(dropped)} out, "
          f"{len(dropped)} unassigned duplicates dropped", file=sys.stderr)


if __name__ == "__main__":
    main()
