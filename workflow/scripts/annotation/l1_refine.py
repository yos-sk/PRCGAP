#!/usr/bin/env python3
"""Set each element's boundary from the 5'end/3'end models, type it, emit the LINE1 BED.

The whole-length query cannot reach the 5'UTR of the older subfamilies -- it is
the most diverged part of the element, so blastn stops short and the chain loses
13-16% of the consensus. Dfam's 61 subfamily 5'end models cover exactly that
region, so the 5' boundary comes from them and the 3' boundary from the 3'end
model that also names the subfamily. The length window is applied to that
boundary, not to the extent of the chain.

Two rejections ride along. Scoring all 67 3'end models rather than only the five
wanted families lets an older element win on its own model and be dropped, which
is what keeps L1PA6/L1PA7 out. And a locus whose top families sit within
--tie-margin of each other is not evidence for one of them, so every family
inside the margin is reported, joined by "/" (L1PA2/L1PA3).

Loci are extracted in element orientation (minus-strand ones reverse
complemented), so inside a locus the 5'end model lands near the start and the
3'end model near the end whatever the contig strand.

Output is the BED6 nanomonsv insert_classify reads: contig, start, end,
"contig,start,end,strand,family", 0, strand.
"""

import argparse
import collections
import sys

YOUNG = ("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5")


def read_hits(path):
    """locus -> [(bits, model, lo, hi)] sorted by score.

    outfmt 6: qseqid sseqid pident length qlen qstart qend sstart send bitscore
    """
    per = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        ss, se = int(f[7]), int(f[8])
        lo, hi = (ss, se) if ss <= se else (se, ss)
        # bedtools getfasta -s -nameOnly names a region "<name>(<strand>)"
        locus = f[1][:-3] if f[1].endswith(("(+)", "(-)")) else f[1]
        per[locus].append((float(f[9]), f[0], lo, hi))
    for k in per:
        per[k].sort(reverse=True)
    return per


def main():
    ap = argparse.ArgumentParser(prog="l1_refine.py")
    ap.add_argument("--loci", required=True, help="padded candidate BED")
    ap.add_argument("--three", required=True)
    ap.add_argument("--five", required=True)
    ap.add_argument("--min-span", type=int, default=5800)
    ap.add_argument("--max-span", type=int, default=6300)
    ap.add_argument("--tie-margin", type=float, default=3.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    three, five = read_hits(args.three), read_hits(args.five)
    counts = collections.Counter()
    rows = []

    for line in open(args.loci):
        # padded candidate BED6: contig, start, end, name, score, strand
        g = line.rstrip("\n").split("\t")
        contig, ps, pe, strand = g[0], int(g[1]), int(g[2]), g[5]
        key = g[3]
        t = three.get(key)
        if not t:
            counts["no_3end_hit"] += 1
            continue

        per_family = {}
        for bits, model, lo, hi in t:
            fam = model.replace("_3end", "")
            if fam not in per_family:
                per_family[fam] = (bits, hi)
        ranked = sorted(per_family.items(), key=lambda kv: -kv[1][0])
        call, (bits, t_hi) = ranked[0]
        if call not in YOUNG:
            counts["rejected_decoy"] += 1
            continue
        tied = [fa for fa, (b, _) in ranked if fa in YOUNG and bits - b <= args.tie_margin]

        upstream = [h for h in five.get(key, ()) if h[2] < t_hi]
        if not upstream:
            counts["no_5end_hit"] += 1
            continue
        start_in_locus, end_in_locus = upstream[0][2] - 1, t_hi

        span = end_in_locus - start_in_locus
        if not (args.min_span <= span <= args.max_span):
            counts["outside_length_window"] += 1
            continue

        if strand == "+":
            cs, ce = ps + start_in_locus, ps + end_in_locus
        else:
            cs, ce = pe - end_in_locus, pe - start_in_locus
        counts["kept"] += 1
        rows.append((contig, cs, ce, "/".join(tied), strand))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, "w") as out:
        for contig, cs, ce, fam, strand in rows:
            name = "{},{},{},{},{}".format(contig, cs, ce, strand, fam)
            out.write("{}\t{}\t{}\t{}\t0\t{}\n".format(contig, cs, ce, name, strand))

    print("[l1_refine] candidates {:,} -> kept {:,}".format(
        sum(counts.values()), counts["kept"]), file=sys.stderr)
    for k, v in counts.most_common():
        if k != "kept":
            print("    {:<22s} {:>6,}".format(k, v), file=sys.stderr)


if __name__ == "__main__":
    main()
