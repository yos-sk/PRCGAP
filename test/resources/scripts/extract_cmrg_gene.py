#!/usr/bin/env python3

gene_set = set()
with open("41587_2021_1158_MOESM4_ESM.tsv", "r") as f:
    for i, line in enumerate(f):
        if i == 0: continue
        items = line.rstrip("\n").split("\t")
        if items[12] != "TRUE": continue
        if items[8].split(":")[0] in ["chrX", "X", "chrY", "Y", "chrMT", "MT"]: continue
        gene_set.add(items[1])

for gene in gene_set:
    print(gene)