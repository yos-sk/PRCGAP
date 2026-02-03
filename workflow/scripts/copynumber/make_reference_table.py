#! /usr/bin/env/ python
# -*- coding: utf-8 -*-

import sys
import gzip
import argparse

class paf_parser:
    def __init__(self):
        self.target_contig = ''
        self.target_length = 0
        self.target_start = 0
        self.target_end = 0
        self.strand = ''
        #self.query_contig = ''
        self.query_length = 0
        self.query_start = 0
        self.query_end = 0
        self.sequence_match = 0
        self.alignment_length = 0
        self.mapping_qutlity = 0
    
    def add_info(self, items: list):
        self.target_contig = items[0]
        self.target_length = int(items[1])
        self.target_start = int(items[2])
        self.target_end = int(items[3])
        self.strand = items[4]
        #self.query_contig = q_c
        self.query_length = int(items[6])
        self.query_start = int(items[7])
        self.query_end = int(items[8])
        self.sequence_match = int(items[9])
        self.alignment_length = int(items[10])
        self.mapping_quality = int(items[11])
    
    def identity(self):
        return self.sequence_match / self.alignment_length
    
    

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

def classify_contig(paf_d: dict, mapq_threshold: int) -> dict:

    unsort_out = list()
    for contig in paf_d:
        tmp_dict = dict()
        prev_contig = ''
        cnt = 0
        for p in sorted(paf_d[contig], key = lambda x: x.query_start):
            if p.mapping_quality < mapq_threshold: continue
            if p.sequence_match < 100000: continue
            
            if prev_contig != '' and prev_contig != p.target_contig:
                cnt += 1
            
            if p.target_contig not in tmp_dict:
                tmp_dict[p.target_contig] = [p]
            else:
                tmp_dict[p.target_contig].append(p)
            
            prev_contig = p.target_contig
        
        #if cnt > 2: continue
        #if len(tmp_dict) > 2: continue
        if len(tmp_dict) == 0: continue
        
        length = 0
        max_contig = ''
        for c in tmp_dict:
            if sum([p.query_end - p.query_start + 1 for p in tmp_dict[c]]) > length:
                length = sum([p.query_end - p.query_start + 1 for p in tmp_dict[c]])
                max_contig = c
        

        f_cnt, r_cnt = 0, 0
        for k in tmp_dict[max_contig]:
            if k.strand == '+': f_cnt += 1
            if k.strand == '-': r_cnt += 1
        if f_cnt > r_cnt:
            unsort_out.append((contig, tmp_dict[max_contig][0].query_start, tmp_dict[max_contig][-1].query_end, '+', max_contig, tmp_dict[max_contig][0].target_start, tmp_dict[max_contig][-1].target_end))
        else:
            if tmp_dict[max_contig][-1].target_start > tmp_dict[max_contig][0].target_end:
                continue
            unsort_out.append((contig, tmp_dict[max_contig][0].query_start, tmp_dict[max_contig][-1].query_end, '-', max_contig, tmp_dict[max_contig][-1].target_start, tmp_dict[max_contig][0].target_end))
        
    prev_chrom = ''
    prev_start = 0
    prev_end = 0
    for i, items in enumerate(sorted(unsort_out, key=lambda x: (x[4], x[5]))):
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
    parser = argparse.ArgumentParser(prog = "make_reference_table.py")
    parser.add_argument('-i', '--input', help="Enter the path of the PAF file.")
    parser.add_argument('-q', '--mapq', type=int, default=20, help="Enter the threshold of mapping quality.")
    
    args = parser.parse_args()
    
    return args
    
def main():
    args = create_parser()
    
    paf_info = load_paf(args.input)
    #print(paf_info['h1tg000001l'], file=sys.stderr)
    classify_contig(paf_info, args.mapq)
    
if __name__ == '__main__':
    main()
    
                    
