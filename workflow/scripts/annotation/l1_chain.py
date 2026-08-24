#!/usr/bin/env python3
"""blastn HSPs of the L1 query against an assembly -> candidate element regions.

A full-length L1 never comes back as one HSP: blastn breaks it at indels, at
diverged stretches, and at whatever has inserted itself into the element. Blocks
are chained back together, but a target-gap rule alone also chains the *next*
copy, because tandem L1 sit a few hundred bp apart. Blocks of one copy are
collinear -- walking along the contig, the consensus interval advances (recedes
on the minus strand) -- and the next copy restarts the consensus, so a break in
collinearity is the copy boundary.

Only candidate regions are produced here. The element boundary is set later from
the 5'end/3'end models (l1_refine.py), so the filters are deliberately loose.
"""

import argparse
import collections
import sys


def read_blast(path):
    """outfmt 6: qseqid sseqid pident length qstart qend sstart send evalue bitscore"""
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        qs, qe = int(f[4]) - 1, int(f[5])
        ss, se = int(f[6]), int(f[7])
        strand = "+" if ss <= se else "-"
        ts, te = (ss - 1, se) if strand == "+" else (se - 1, ss)
        alnlen = int(f[3])
        yield dict(contig=f[1], tstart=ts, tend=te, strand=strand,
                   qstart=qs, qend=qe, alnlen=alnlen,
                   nmatch=int(round(alnlen * float(f[2]) / 100.0)))


def union_len(intervals):
    out, end = 0, -1
    for s, e in sorted(intervals):
        if s > end:
            out += e - s
            end = e
        elif e > end:
            out += e - end
            end = e
    return out


def main():
    ap = argparse.ArgumentParser(prog="l1_chain.py")
    ap.add_argument("--input", required=True)
    ap.add_argument("--query-len", type=int, required=True)
    ap.add_argument("--max-gap", type=int, default=300,
                    help="contig gap allowed between blocks of one element")
    ap.add_argument("--q-tolerance", type=int, default=200,
                    help="consensus bp a block may step back and still count as "
                         "collinear; adjacent HSPs of one element overlap a "
                         "little where blastn extended both sides of an indel")
    ap.add_argument("--min-cov", type=float, default=0.55,
                    help="fraction of the query the chain must cover")
    ap.add_argument("--min-span", type=int, default=3000)
    ap.add_argument("--max-span", type=int, default=6300)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    blocks = collections.defaultdict(list)
    n_in = 0
    for b in read_blast(args.input):
        n_in += 1
        blocks[(b["contig"], b["strand"])].append(b)

    loci = []
    tol = args.q_tolerance
    for (contig, strand), bl in blocks.items():
        bl.sort(key=lambda b: b["tstart"])
        chain, q_hi, q_lo = [bl[0]], bl[0]["qend"], bl[0]["qstart"]
        for b in bl[1:]:
            collinear = (b["qstart"] >= q_hi - tol if strand == "+"
                         else b["qend"] <= q_lo + tol)
            # an element cannot exceed max_span, so a block that would push the
            # chain past it belongs to something else; without this a short
            # low-identity fragment a kilobase away inflates an otherwise
            # correct locus out of the length window
            fits = b["tend"] - chain[0]["tstart"] <= args.max_span
            if (b["tstart"] - chain[-1]["tend"] <= args.max_gap
                    and collinear and fits):
                chain.append(b)
                q_hi, q_lo = max(q_hi, b["qend"]), min(q_lo, b["qstart"])
            else:
                loci.append((contig, strand, chain))
                chain, q_hi, q_lo = [b], b["qend"], b["qstart"]
        loci.append((contig, strand, chain))

    kept = []
    for contig, strand, chain in loci:
        ts = min(b["tstart"] for b in chain)
        te = max(b["tend"] for b in chain)
        frac = union_len([(b["qstart"], b["qend"]) for b in chain]) / args.query_len
        if frac < args.min_cov:
            continue
        if not (args.min_span <= te - ts <= args.max_span):
            continue
        alnlen = sum(b["alnlen"] for b in chain) or 1
        kept.append((contig, ts, te, strand, round(frac, 4), len(chain),
                     round(sum(b["nmatch"] for b in chain) / alnlen, 4)))

    kept.sort(key=lambda x: (x[0], x[1]))
    with open(args.out, "w") as out:
        for r in kept:
            out.write("\t".join(str(v) for v in r) + "\n")

    print("[l1_chain] {:,} HSP -> {:,} chain -> {:,} candidate loci".format(
        n_in, len(loci), len(kept)), file=sys.stderr)


if __name__ == "__main__":
    main()
