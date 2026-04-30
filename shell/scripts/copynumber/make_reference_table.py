#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse


class paf_parser:
    def __init__(self):
        self.target_contig = ''
        self.target_length = 0
        self.target_start = 0
        self.target_end = 0
        self.strand = ''
        self.query_length = 0
        self.query_start = 0
        self.query_end = 0
        self.sequence_match = 0
        self.alignment_length = 0
        self.mapping_quality = 0

    def add_info(self, items: list):
        self.target_contig = items[0]
        self.target_length = int(items[1])
        self.target_start = int(items[2])
        self.target_end = int(items[3])
        self.strand = items[4]
        self.query_length = int(items[6])
        self.query_start = int(items[7])
        self.query_end = int(items[8])
        self.sequence_match = int(items[9])
        self.alignment_length = int(items[10])
        self.mapping_quality = int(items[11])


def load_paf(input_paf: str) -> dict:
    out = dict()
    with open(input_paf, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            q_c = items[5]
            p = paf_parser()
            p.add_info(items)
            if q_c in out:
                out[q_c].append(p)
            else:
                out[q_c] = [p]
    return out


def classify_contig(paf_d: dict, mapq_threshold: int, bin_size: int, coverage_threshold: float,
                    min_span: int = 1000000, min_total_aligned: int = 5000000) -> None:
    results = []

    for contig in paf_d:
        # Filter alignments by quality and length
        filtered = [p for p in paf_d[contig]
                     if p.mapping_quality >= mapq_threshold and p.sequence_match >= 100000]
        if not filtered:
            continue

        # Group filtered alignments by reference chromosome
        chrom_alns = {}
        for p in filtered:
            chrom_alns.setdefault(p.target_contig, []).append(p)

        # For each chromosome, use bin-based coverage to decide assignment
        for chrom, alns in chrom_alns.items():
            chrom_length = alns[0].target_length
            n_bins = (chrom_length + bin_size - 1) // bin_size

            # Calculate alignment coverage per bin on the reference
            bin_coverage = [0] * n_bins
            for p in alns:
                start_bin = p.target_start // bin_size
                end_bin = min((p.target_end - 1) // bin_size, n_bins - 1)
                for b in range(start_bin, end_bin + 1):
                    bin_start = b * bin_size
                    bin_end = min((b + 1) * bin_size, chrom_length)
                    overlap_start = max(p.target_start, bin_start)
                    overlap_end = min(p.target_end, bin_end)
                    if overlap_end > overlap_start:
                        bin_coverage[b] += overlap_end - overlap_start

            # Identify bins with coverage exceeding the threshold
            covered_bins = []
            for i in range(n_bins):
                bin_len = min(bin_size, chrom_length - i * bin_size)
                if bin_coverage[i] >= coverage_threshold * bin_len:
                    covered_bins.append(i)

            if not covered_bins:
                continue

            # Merge contiguous covered bins into segments
            segments = []
            seg_start = covered_bins[0]
            seg_end = covered_bins[0]
            for b in covered_bins[1:]:
                if b == seg_end + 1:
                    seg_end = b
                else:
                    segments.append((seg_start, seg_end))
                    seg_start = b
                    seg_end = b
            segments.append((seg_start, seg_end))

            # For each segment, compute precise coordinates from actual alignments
            for seg_s, seg_e in segments:
                seg_ref_start = seg_s * bin_size
                seg_ref_end = min((seg_e + 1) * bin_size, chrom_length)

                # Collect alignments overlapping this segment
                seg_alns = [p for p in alns
                            if p.target_end > seg_ref_start and p.target_start < seg_ref_end]
                if not seg_alns:
                    continue

                # Precise coordinates from actual alignment boundaries
                query_start = min(p.query_start for p in seg_alns)
                query_end = max(p.query_end for p in seg_alns)
                ref_start = min(p.target_start for p in seg_alns)
                ref_end = max(p.target_end for p in seg_alns)

                # Strand by majority vote
                f_cnt = sum(1 for p in seg_alns if p.strand == '+')
                r_cnt = sum(1 for p in seg_alns if p.strand == '-')
                strand = '+' if f_cnt > r_cnt else '-'

                results.append((contig, query_start, query_end, strand, chrom, ref_start, ref_end))

    # Filter by minimum reference span before clustering
    results = [r for r in results if r[6] - r[5] >= min_span]

    # Fallback for sparse alignments: if a (contig, chrom) pair has sufficient
    # total aligned bases but the surviving records inadequately represent it,
    # replace with a full-range record
    surviving_span = {}
    for r in results:
        key = (r[0], r[4])
        surviving_span[key] = surviving_span.get(key, 0) + (r[6] - r[5])

    for contig in paf_d:
        filtered = [p for p in paf_d[contig]
                     if p.mapping_quality >= mapq_threshold and p.sequence_match >= 100000]
        if not filtered:
            continue
        chrom_alns = {}
        for p in filtered:
            chrom_alns.setdefault(p.target_contig, []).append(p)
        for chrom, alns in chrom_alns.items():
            total_aligned = sum(p.target_end - p.target_start for p in alns)
            if total_aligned < min_total_aligned:
                continue
            key = (contig, chrom)
            current_span = surviving_span.get(key, 0)
            if current_span >= total_aligned * 0.5:
                continue
            # Remove existing records for this pair and add full-range fallback
            results = [r for r in results if not (r[0] == contig and r[4] == chrom)]
            query_start = min(p.query_start for p in alns)
            query_end = max(p.query_end for p in alns)
            ref_start = min(p.target_start for p in alns)
            ref_end = max(p.target_end for p in alns)
            f_cnt = sum(1 for p in alns if p.strand == '+')
            r_cnt = sum(1 for p in alns if p.strand == '-')
            strand = '+' if f_cnt > r_cnt else '-'
            results.append((contig, query_start, query_end, strand, chrom, ref_start, ref_end))

    # Cluster records by (contig, chrom, strand) and merge coordinates
    clustered = {}
    for r in results:
        key = (r[0], r[4], r[3])  # (contig, chrom, strand)
        if key in clustered:
            clustered[key] = (
                key[0],
                min(clustered[key][1], r[1]),
                max(clustered[key][2], r[2]),
                key[2],
                key[1],
                min(clustered[key][5], r[5]),
                max(clustered[key][6], r[6]),
            )
        else:
            clustered[key] = r

    filtered = list(clustered.values())

    # Sort by (ref_chrom, ref_start) and remove contained records
    merged = sorted(filtered, key=lambda x: (x[4], x[5]))
    prev_chrom = ''
    prev_start = 0
    prev_end = 0
    for i, items in enumerate(merged):
        if i == 0:
            prev_chrom = items[4]
            prev_start = items[5]
            prev_end = items[6]
            print('\t'.join([str(k) for k in items]))
        else:
            if prev_chrom == items[4] and prev_start <= items[5] and items[6] <= prev_end:
                continue
            print('\t'.join([str(k) for k in items]))
            prev_chrom = items[4]
            prev_start = items[5]
            prev_end = items[6]


def create_parser():
    parser = argparse.ArgumentParser(prog="make_reference_table_v2.py")
    parser.add_argument('-i', '--input', required=True, help="Path to the PAF file.")
    parser.add_argument('-q', '--mapq', type=int, default=20, help="Mapping quality threshold (default: 20).")
    parser.add_argument('-b', '--bin_size', type=int, default=1000000,
                        help="Bin size in bp for coverage calculation (default: 1000000).")
    parser.add_argument('-c', '--coverage_threshold', type=float, default=0.8,
                        help="Minimum coverage fraction per bin (default: 0.8).")
    parser.add_argument('-s', '--min_span', type=int, default=1000000,
                        help="Minimum reference span to output after clustering (default: 1000000).")
    parser.add_argument('-t', '--min_total_aligned', type=int, default=5000000,
                        help="Minimum total aligned bases for fallback assignment (default: 5000000).")
    return parser.parse_args()


def main():
    args = create_parser()
    paf_info = load_paf(args.input)
    classify_contig(paf_info, args.mapq, args.bin_size, args.coverage_threshold, args.min_span,
                    args.min_total_aligned)


if __name__ == '__main__':
    main()
