#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to count SV types from a VCF file (e.g., from Sniffles).
Usage: python count_sv_types.py input.vcf output.csv
"""

import sys
import pandas as pd

def count_svtypes(vcf_file):
    svtype_counts = {}

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info_field = fields[7]
            info_items = info_field.split(';')
            for item in info_items:
                if item.startswith("SVTYPE="):
                    svtype = item.split("=")[1]
                    svtype_counts[svtype] = svtype_counts.get(svtype, 0) + 1
                    break  # SVTYPE found, move to next line

    return svtype_counts

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python count_sv_types.py <input.vcf> <output.csv>")
        sys.exit(1)

    vcf_input = sys.argv[1]
    output_csv = sys.argv[2]

    counts = count_svtypes(vcf_input)
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['Count'])
    df.index.name = 'SVTYPE'
    df.to_csv(output_csv)


