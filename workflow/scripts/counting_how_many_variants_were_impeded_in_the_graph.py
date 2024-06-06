#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:12:45 2024

@author: yomna
"""
import sys
samples = 3
#int(sys.argv[1])
input_file = "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/results/merged_vg_combined_table.txt"
#sys.argv[2]
#output_dir=os.path.dirname(sys.argv[3])
number_of_variants=0
with open(input_file) as file:
   for r in file:
        z = r.split('\t')
        if "Chromosom" not in z:
            genotype = z[5:5+samples]
            allele=[]
            for g in genotype:
                if g != "0" and g != ".":
                   allele.append(g)
            number_of_variants += len(set(allele))
print(number_of_variants)             
