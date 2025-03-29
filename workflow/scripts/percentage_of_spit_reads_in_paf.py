#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 22:29:54 2025

@author: yomna
"""

import pandas as pd
import sys
# Function to check full containment
def is_fully_contained(start1, end1, start2, end2):
    return start1 <= start2 and end1 >= end2

# Read BED file and process
bed_file = sys.argv[1]#"/home/yomna/hpc_project/ToxoVar/analysis/Sniffles/2020T/sniffles_2020T_with_reference_corrected.bed"  # Replace with your BED file path
paf_file = sys.argv[2] #"/home/yomna/hpc_project/ToxoVar/analysis/Sniffles/2020T_trial/alignments.paf"  # Replace with your PAF file path

# Load BED positions
bed_positions = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "label"])

# Load only required columns from PAF file
use_columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]  # Column indices for required fields
paf_columns = [
    "query_name", "query_length", "query_start", "query_end",
    "strand", "target_name", "target_length", "target_start", "target_end",
    "residue_matches", "alignment_block_length", "mapping_quality", "NM", "ms", "AS", "nn", "tp"
]
paf_data = pd.read_csv(paf_file, sep="\t", header=None, usecols=use_columns, names=paf_columns)

# Exclude secondary alignments (tag tp:A:S)
paf_data = paf_data[paf_data["tp"] != "tp:A:S"]

# Group alignments by query_name to include all split parts of the same read
grouped_alignments = paf_data.groupby("query_name")

results = []

for index, row in bed_positions.iterrows():
    chrom, start, end, label = row["chrom"], row["start"], row["end"], row["label"]
    chrom
    # Filter PAF entries for the current region
    region_reads = paf_data[(paf_data["target_name"] == chrom) &
                            (paf_data["target_start"] < end) &
                            (paf_data["target_end"] > start)]

    total_reads = len(region_reads["query_name"].unique())
    split_reads = 0
    hits = 0

    for query_name in region_reads["query_name"].unique():
        # Get all alignments for this query from the grouped data
        query_alignments = grouped_alignments.get_group(query_name)

        if len(query_alignments) > 1:  # Indicates the read is split
            split_reads += 1

            # Sort alignments: primary alignment (longest alignment_block_length) first
            query_alignments = query_alignments.sort_values(by="alignment_block_length", ascending=False)

            # Take the primary and supplementary alignments
            primary_alignment = query_alignments.iloc[0]
            supplementary_alignments = query_alignments.iloc[1:]

            for _, supp_alignment in supplementary_alignments.iterrows():
                # Check if the supplementary alignment is fully contained within the primary alignment
                if is_fully_contained(
                    primary_alignment["target_start"],
                    primary_alignment["target_end"],
                    supp_alignment["target_start"],
                    supp_alignment["target_end"]
                ) and primary_alignment["strand"] != supp_alignment["strand"]:
                    hits += 1
                    break

    # Calculate proportions
    split_proportion = split_reads / total_reads*100 if total_reads > 0 else 0
    hits_proportion = hits / total_reads *100 if total_reads > 0 else 0
    # Save results for this region
    results.append({
        "label": label,
        "chrom": chrom,
        "start": start,
        "end": end,
        "total_reads": total_reads,
        "split_reads": split_reads,
        "split_proportion": split_proportion,
        "hits": hits,
        "hits_proportion": hits_proportion  
    })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to a CSV file
results_df.to_csv(sys.argv[3], index=False)

print("Analysis completed. Results saved to output_results.csv.")