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
        self.GRCh38_status = ""
        self.chm13_contig = ""
        self.chm13_pos = ""
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
        self.GRCh38_status = info[7]
        self.chm13_contig = info[8]
        self.chm13_pos = info[9]
        self.chm13_status = info[10]
        self.GRCh38 = info[11]
        self.chm13 = info[12]
        self.vaf = float(info[13])
        self.depth = int(info[14])
        self.normal_vaf = float(info[15])
        self.normal_depth = int(info[16])
        self.mutation_id = info[17]
        self.haplotype = info[18]
        self.kmer_ratio = float(info[19])
        self.liftoff_gene = info[20]
        self.liftoff_info = info[21]
        self.cgc = info[22]
        self.cmrg = info[23]
        self.rmsk = info[24]
        self.rmsk_type = info[25]
        self.contig_size = int(info[26])
        self.misassembly = info[27]
        self.centromere = info[28]
        self.segdup = info[29]
        self.segdup_similarity = info[30]
        self.other = info[31]
    
    
    def output(self):
        print(self.contig, self.start, self.end, self.ref, self.alt, self.GRCh38_contig, self.GRCh38_pos, self.GRCh38_status,
              self.chm13_contig, self.chm13_pos, self.chm13_status, self.GRCh38, self.chm13, self.vaf, self.depth, self.normal_vaf, self.normal_depth, self.mutation_id, self.haplotype, self.kmer_ratio,
              self.liftoff_gene, self.liftoff_info, self.cgc, self.cmrg, self.rmsk, self.rmsk_type, self.contig_size, self.misassembly, self.centromere, self.segdup, self.segdup_similarity, self.other, sep="\t")
    
    def complement(self):
        comp_dict = {"A": "T",
                     "C": "G",
                     "G": "C",
                     "T": "A"}
        
        return comp_dict[self.ref], comp_dict[self.alt]

def get_haplotype(contig, haplotype):
    if haplotype == "homozygous":
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

            comp_ref, comp_alt = info.complement()
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

            elif prev_info.ref == info.ref and prev_info.alt == info.alt:
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
                            prev_info.output()
                        else:
                            info.output()
                        prev_info = SNVinfo()
                continue
            elif prev_info.ref == comp_ref and prev_info.alt == comp_alt:
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
                            prev_info.output()
                        else:
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

            comp_ref, comp_alt = info.complement()
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

            elif prev_info.ref == info.ref and prev_info.alt == info.alt:
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
                            prev_info.output()
                        else:
                            info.output()
                        prev_info = SNVinfo()
                continue
            elif prev_info.ref == comp_ref and prev_info.alt == comp_alt:
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
                            prev_info.output()
                        else:
                            info.output()
                        prev_info = SNVinfo()
                continue

            prev_info = info 

        if prev_info.contig != "":
            prev_info.output()

if __name__ == "__main__":
    main()
