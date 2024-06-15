#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:12:45 2024

@author: yomna
"""
import sys
samples = 3
#int(sys.argv[1])
input_file = "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/merge_concat_medaka_sniffles/merge_medaka_sniffles.vcf"
#sys.argv[2]
#output_dir=os.path.dirname(sys.argv[3])
number_of_variants=0
with open(input_file) as file:
   for r in file:
        if "#" not in r:
            z = r.split('\t')
            var = z[4]
            number_of_variants += len(var)
print(number_of_variants)             
