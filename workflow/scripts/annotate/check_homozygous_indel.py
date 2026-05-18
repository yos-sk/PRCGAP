#!/usr/bin/env python3

import sys

class SNVinfo:
    def __init__(self):
        self.contig = ""
        self.start = 0
        self.end = 0
        self.ref = ""
        self.alt = ""
        self.GRCh38_contig = ""
        self.GRCh38_pos = ""
        self.GRCh38_ref = ""
        self.GRCh38_alt = ""
        self.GRCh38_status = ""
        self.chm13_contig = ""
        self.chm13_pos = ""
        self.chm13_ref = ""
        self.chm13_alt = ""
        self.chm13_status = ""
        self.GRCh38 = ""
        self.chm13 = ""
        self.vaf = 0.0
        self.depth = 0
        self.normal_vaf = 0.0
        self.normal_depth = 0
        self.mutation_id = ""
        self.haplotype = ""
        self.kmer_ratio = 0.0
        self.liftoff_gene = ""
        self.liftoff_info = ""
        self.cgc = ""
        self.cmrg = ""
        self.rmsk = ""
        self.rmsk_type = ""
        self.contig_size = 0
        self.misassembly = ""
        self.centromere = ""
        self.segdup = ""
        self.segdup_similarity = ""
        self.point_mutation_other = ""
        
    def insert(self, info):
        self.contig = info[0]
        self.start = int(info[1])
        self.end = int(info[2])
        self.ref = info[3]
        self.alt = info[4]
        self.GRCh38_contig = info[5]
        self.GRCh38_pos = info[6]
        self.GRCh38_ref = info[7]
        self.GRCh38_alt = info[8]
        self.GRCh38_status = info[9]
        self.chm13_contig = info[10]
        self.chm13_pos = info[11]
        self.chm13_ref = info[12]
        self.chm13_alt = info[13]
        self.chm13_status = info[14]
        self.GRCh38 = info[15]
        self.chm13 = info[16]
        self.vaf = float(info[17])
        self.depth = int(info[18])
        self.normal_vaf = float(info[19])
        self.normal_depth = int(info[20])
        self.mutation_id = info[21]
        self.haplotype = info[22]
        self.kmer_ratio = float(info[23])
        self.liftoff_gene = info[24]
        self.liftoff_info = info[25]
        self.cgc = info[26]
        self.cmrg = info[27]
        self.rmsk = info[28]
        self.rmsk_type = info[29]
        self.contig_size = int(info[30])
        self.misassembly = info[31]
        self.centromere = info[32]
        self.segdup = info[33]
        self.segdup_similarity = info[34]
        self.other = info[35]
    
    
    def output(self):
        print(self.contig, self.start, self.end, self.ref, self.alt, self.GRCh38_contig, self.GRCh38_pos, self.GRCh38_ref, self.GRCh38_alt, self.GRCh38_status,
              self.chm13_contig, self.chm13_pos, self.chm13_ref, self.chm13_alt, self.chm13_status, self.GRCh38, self.chm13, self.vaf, self.depth, self.normal_vaf, self.normal_depth, self.mutation_id, self.haplotype,
              self.kmer_ratio, self.liftoff_gene, self.liftoff_info, self.cgc, self.cmrg, self.rmsk, self.rmsk_type, self.contig_size, self.misassembly, self.centromere, self.segdup, self.segdup_similarity, self.other, sep="\t")
    
def get_haplotype(contig, haplotype):
    if haplotype.startswith("haplotype"):
        return haplotype
    elif haplotype == "homozygous":
        if contig.startswith("haplotype1") or contig.startswith("h1"):
            return "haplotype1"
        else:
            return "haplotype2"
    else:
        return haplotype

def main():
    input_file = sys.argv[1]
    reference = sys.argv[2]
    mode = sys.argv[3]
    
    info_list = list()
    header = ""
    with open(input_file, "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                header = line.rstrip("\n")
                continue
            info = SNVinfo()
            items = line.rstrip("\n").split("\t")
            info.insert(items)
            info_list.append(info)
    
    print(header)
    prev_info = SNVinfo()
    if reference == "GRCh38":
        for i, info in enumerate(sorted(info_list, key=lambda x: (x.GRCh38_contig, x.GRCh38_pos))):
            if len(prev_info.contig) == 0:
                prev_info = info
                continue

            if prev_info.GRCh38_contig != info.GRCh38_contig:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif prev_info.GRCh38_pos != info.GRCh38_pos:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif "-" in info.GRCh38_pos:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif prev_info.GRCh38_ref == info.GRCh38_ref and prev_info.GRCh38_alt == info.GRCh38_alt:
                if prev_info.haplotype.startswith("haplotype") or prev_info.haplotype == "homozygous":
                    if info.haplotype.startswith("haplotype") or info.haplotype == "homozygous":
                        prev_info.haplotype = "homozygous"
                        info.haplotype = "homozygous"
                        prev_info.output()
                        info.output()
                        prev_info = SNVinfo()
                    else:
                        if mode == "update":
                            prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                        prev_info.output()
                        prev_info = SNVinfo()
                else:
                    if info.haplotype.startswith("haplotype") or info.haplotype == "homozygous":
                        if mode == "update":
                            info.haplotype = get_haplotype(info.contig, info.haplotype)
                        info.output()
                        prev_info = SNVinfo()
                    else:
                        if prev_info.vaf > info.vaf:
                            prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                            prev_info.output()
                        else:
                            info.haplotype = get_haplotype(info.contig, info.haplotype)
                            info.output()
                        prev_info = SNVinfo()
                continue

            prev_info = info 

        if prev_info.contig != "":
            prev_info.output()

    else:
        for i, info in enumerate(sorted(info_list, key=lambda x: (x.chm13_contig, x.chm13_pos))):
            if len(prev_info.contig) == 0:
                prev_info = info
                continue

            if prev_info.chm13_contig != info.chm13_contig:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif prev_info.chm13_pos != info.chm13_pos:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif "-" in info.chm13_pos:
                if mode == "update":
                    prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                prev_info.output()

            elif prev_info.chm13_ref == info.chm13_ref and prev_info.chm13_alt == info.chm13_alt:
                if prev_info.haplotype.startswith("haplotype") or prev_info.haplotype == "homozygous":
                    if info.haplotype.startswith("haplotype") or info.haplotype == "homozygous":
                        prev_info.haplotype = "homozygous"
                        info.haplotype = "homozygous"
                        prev_info.output()
                        info.output()
                        prev_info = SNVinfo()
                    else:
                        if mode == "update":
                            prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                        prev_info.output()
                        prev_info = SNVinfo()
                else:
                    if info.haplotype.startswith("haplotype") or info.haplotype == "homozygous":
                        if mode == "update":
                            info.haplotype = get_haplotype(info.contig, info.haplotype)
                        info.output()
                        prev_info = SNVinfo()
                    else:
                        if prev_info.vaf > info.vaf:
                            prev_info.haplotype = get_haplotype(prev_info.contig, prev_info.haplotype)
                            prev_info.output()
                        else:
                            info.haplotype = get_haplotype(info.contig, info.haplotype)
                            info.output()
                        prev_info = SNVinfo()
                continue

            prev_info = info 

        if prev_info.contig != "":
            prev_info.output()

if __name__ == "__main__":
    main()
