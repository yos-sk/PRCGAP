#!/usr/bin/env python3
"""Extract chr20 MANE-registered transcripts from a GENCODE GFF, with a
longest-exon fallback for genes not covered by MANE.

Adapted from gcat_database_scripts/scripts/GENCODE_v46/proc_gencode_bed.py.
NCC Select is dropped (test scope), but the longest-exon fallback is kept
so genes like readthroughs (e.g. PEDS1-UBE2V1), recently annotated
protein-coding genes (e.g. GCNT7), and Ensembl-only loci (e.g.
ENSG00000288000) — none of which appear in MANE v1.3 — still emit a
representative transcript.

Output: <prefix>.transcript.bed  (1-based GFF -> 0-based BED start)
"""

import sys
import gzip


def _open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _parse_attrs(attr_field):
    attrs = {}
    for kv in attr_field.split(";"):
        if "=" not in kv: continue
        k, v = kv.split("=", 1)
        attrs[k] = v
    return attrs


def get_mane_info(mane_summary_file):
    mane_info = {}
    with _open_text(mane_summary_file) as hin:
        for line in hin:
            if line.startswith("#"): continue
            F = line.rstrip("\n").split("\t")
            transcript_id = F[7]
            mane_status = F[9]
            key = transcript_id.split(".")[0]
            mane_info.setdefault(key, []).append(mane_status)
    return mane_info


def get_longest_exon_fallback(in_gff_file, mane_info):
    """Pick one representative transcript per chr20 gene that has no
    transcript registered in MANE. The representative is the transcript
    with the largest total exon length (sum of exon widths).
    """
    # gene_name is empty for some Ensembl-only loci; fall back to gene_id
    # to avoid collapsing all such entries onto a single empty-string key.
    transcript_gene = {}        # transcript_id (with version) -> gene_key
    transcript_exome_len = {}   # transcript_id (with version) -> int
    gene_has_mane = set()       # gene_keys with at least one MANE transcript

    with gzip.open(in_gff_file, "rt") as hin:
        for line in hin:
            if line.startswith("#"): continue
            F = line.rstrip("\n").split("\t")
            if F[0] != "chr20": continue

            if F[2] == "transcript":
                a = _parse_attrs(F[8])
                transcript_id = a.get("transcript_id", "")
                if not transcript_id: continue
                gene_key = a.get("gene_name") or a.get("gene_id", "")
                transcript_gene[transcript_id] = gene_key
                transcript_exome_len.setdefault(transcript_id, 0)
                if transcript_id.split(".")[0] in mane_info:
                    gene_has_mane.add(gene_key)

            elif F[2] == "exon":
                a = _parse_attrs(F[8])
                parent = a.get("transcript_id", "")
                if not parent: continue
                transcript_exome_len[parent] = (
                    transcript_exome_len.get(parent, 0) + int(F[4]) - int(F[3]) + 1
                )

    longest_per_gene = {}  # gene_key -> (transcript_id, length)
    for tid, gene_key in transcript_gene.items():
        if gene_key in gene_has_mane: continue
        length = transcript_exome_len.get(tid, 0)
        cur = longest_per_gene.get(gene_key)
        if cur is None or length >= cur[1]:
            longest_per_gene[gene_key] = (tid, length)

    return {tid for tid, _ in longest_per_gene.values()}


def write_chr20_bed(in_gff_file, mane_info, longest_fallback, output_file_prefix):
    with gzip.open(in_gff_file, "rt") as hin, \
         open(f"{output_file_prefix}.transcript.bed", "w") as hout:
        for line in hin:
            if line.startswith("#"): continue
            F = line.rstrip("\n").split("\t")
            if F[0] != "chr20": continue
            if F[2] != "transcript": continue

            a = _parse_attrs(F[8])
            gene_name = a.get("gene_name", "")
            transcript_id = a.get("transcript_id", "")
            exon_num = a.get("exon_number", "0")
            gene_id = a.get("gene_id", "")
            hgnc_id = a.get("hgnc_id", "")

            tid_without_ver = transcript_id.split(".")[0]
            in_mane = tid_without_ver in mane_info
            in_fallback = transcript_id in longest_fallback
            if not in_mane and not in_fallback: continue

            select_flag = ";".join(mane_info[tid_without_ver]) if in_mane else ""
            print(
                f"{F[0]}\t{int(F[3]) - 1}\t{F[4]}\t{transcript_id}\t{exon_num}\t{F[6]}\t{gene_name}\t{gene_id}\t{hgnc_id}\t{select_flag}",
                file=hout,
            )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: proc_gencode_bed_mane_chr20.py <gencode.gff3.gz> <MANE.summary.txt[.gz]> <output_prefix>",
            file=sys.stderr,
        )
        sys.exit(1)

    gff_file = sys.argv[1]
    mane_summary_file = sys.argv[2]
    output_file_prefix = sys.argv[3]

    mane_info = get_mane_info(mane_summary_file)
    longest_fallback = get_longest_exon_fallback(gff_file, mane_info)
    write_chr20_bed(gff_file, mane_info, longest_fallback, output_file_prefix)
