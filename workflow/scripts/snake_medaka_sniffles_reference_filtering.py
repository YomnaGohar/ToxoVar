#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 12:21:01 2024

@author: yomna
"""

import sys
def check_variant_presence(sample_id, id_to_check):
    for id1 in id_to_check.split(","):
        for i in id1.split(":"):
            if i in sample_id:
                return True  # Found the variant, exit the function immediately

    return False  


def main():
    combined_table = sys.argv[1]
    num_samples = int(sys.argv[2])
    samples =sys.argv[3:3+num_samples]
    medaka_files = sys.argv[3+num_samples:3+2*num_samples]
    sniffles_files = sys.argv[3+2*num_samples:3+3*num_samples]
    outfile = sys.argv[-1]
    sample_dict = dict(zip(samples, zip(medaka_files, sniffles_files)))
    ids_dict = {}
    for sample, (medaka_vcf, sniffles_vcf) in sample_dict.items():
        ids_dict[sample] = set()
        for file_path in (medaka_vcf, sniffles_vcf):
            with open(file_path, "r") as file:
                for line in file:
                    if line.startswith('#'):
                        continue  # Skip header lines
                    
                    parts = line.strip().split("\t")
                    if len(parts) > 2: 
                       ids_dict[sample].add(parts[2])  # Add the 
    l=0
    with open(combined_table, 'r') as f, open(outfile, 'w') as out:
        for line in f:
            if not (line.startswith("#") or line.startswith("##")):
                l=l+1
                #print(l)
                parts = line.strip().split("\t")
                id_to_check = parts[2].split("=")[1]
                chrom=parts[0]
                i=9
                for sample in samples: 
                    if parts[i] == "0":
                       variant_present = check_variant_presence(ids_dict[sample] , id_to_check)
                       if variant_present:
                           parts[i]= "."
                       else:
                           parts[i] = "0"
                       i=i+1    
                out.write("\t".join(parts) + "\n")
            else:
                out.write(line)
if __name__ == "__main__":
    main()
