#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 12:21:01 2024

@author: yomna
"""

import sys

def read_agp(agp_infile):
    exclude_contigs = set()
    with open(agp_infile, "r") as agp:
        for line in agp:
            if line.startswith("TgMe49_00") or line.startswith("TgMe49_API"):
                parts = line.split("\t")
                if parts[4] == "W":
                    exclude_contigs.add(parts[5])
    return exclude_contigs
def check_variant_presence(files_to_check, id_to_check):
    ids2=set()
    #files_to_check = (medaka_vcf, sniffles_vcf)
    for file_path in files_to_check:
        with open(file_path, "r") as file:
            for line in file:
                if line.startswith('#'):
                    continue  # Skip header lines
                parts = line.strip().split("\t")
                ids2.add(parts[2])  # Assuming the ID is in the third column (0-based index)
    for id1 in id_to_check.split(","):
        for i in id1.split(":"):
            if i in ids2:
                return True  # Found the variant, exit the function immediately

    return False  


def main():
    combined_table = sys.argv[1]
    #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/results/merged_vg_combined_table.txt"
    #
    num_samples = int(sys.argv[2])
    #1
    #
    samples =sys.argv[3:3+num_samples]
    #["2015T"]
    #
    medaka_files = sys.argv[3+num_samples:3+2*num_samples]
    #["/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/medaka/medaka_2015T/medaka.annotated_with_VarType_new_head.assignedID.vcf"]
    #
    sniffles_files = sys.argv[3+2*num_samples:3+3*num_samples]
    #["/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Sniffles/2015T/sniffles_2015T_with_reference_corrected.newHead_assignedID.vcf"]
    #
    agp_file = sys.argv[-2]
    #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/Data/t-gondii-me49-new-genome-annotation_results/pseudo.pseudochr.agp"
    #
    outfile = sys.argv[-1]
    #"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/trail.txt"
    #

    exclude_contigs = read_agp(agp_file)
    sample_dict = dict(zip(samples, zip(medaka_files, sniffles_files)))
    l=0
    with open(combined_table, 'r') as f, open(outfile, 'w') as out:
        headers = next(f).strip()  # Presuming the first line is headers
        out.write(headers + "\n")
        for line in f:
            l=l+1
            print(l)
            parts = line.strip().split("\t")
            id_to_check = parts[2].split("=")[1]
            chrom=parts[0]
            if chrom in exclude_contigs:
                continue  # Skip processing this row altogether
            i=6
            for sample, (medaka_vcf, sniffles_vcf)in sample_dict.items(): 
                if parts[i] == "0":
                   variant_present = check_variant_presence( (medaka_vcf, sniffles_vcf), id_to_check)
                   if variant_present:
                       parts[i]= "."
                   else:
                       parts[i] = "0"
                   i=i+1    
            out.write("\t".join(parts) + "\n")

if __name__ == "__main__":
    main()
