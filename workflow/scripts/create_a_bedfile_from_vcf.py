#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 18:11:39 2025

@author: yomna
"""

import sys
import re
def vcf_to_bed(vcf_file, bed_file, window_size):
    """
    Convert a VCF file to a BED file, with a specified window size around each variant.

    Args:
        vcf_file (str): Path to the input VCF file.
        bed_file (str): Path to the output BED file.
        window_size (int): Number of bases to include around each variant.

    Returns:
        None
    """
    with open(vcf_file, 'r') as vcf, open(bed_file, 'w') as bed:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines

            cols = line.strip().split("\t")
            chrom = cols[0]  # Chromosome
            pos = int(cols[1])  # Position (1-based)
            label = cols[2]
            start = max(0, pos - window_size - 1)  # 0-based start
            info_field = cols[7]
            match = re.search(r"END=(\d+)", info_field)
            if "DEL" in label:
                alt= len(cols[3])
                end= pos+alt + window_size
            elif "INV" in cols[4]:
                  end =  int(match.group(1)) + window_size
            else:      
                 end = pos + window_size  # 1-based end

            bed.write(f"{chrom}\t{start}\t{end}\t{label}_{chrom}_{pos}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python vcf_to_bed_with_window.py <input.vcf> <output.bed> <window_size>")
        sys.exit(1)

    vcf_file = sys.argv[1]  # Input VCF file
    bed_file = sys.argv[2]  # Output BED file
    window_size = int(sys.argv[3])  # Window size around the variant

    vcf_to_bed(vcf_file, bed_file, window_size)
    print(f"Converted {vcf_file} to {bed_file} with a window size of {window_size} bases.")
