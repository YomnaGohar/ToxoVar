#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 23:14:18 2025

@author: yomna
"""

import pysam

def group_reads_by_cigar(bam_file_before, top_n=10):
    # Open BAM files
    bam_before = pysam.AlignmentFile(bam_file_before, "rb")

    # Region of interest
    contig_b = "contig_46"
    start_b = 7298291
    end_b = 7327278

    def get_top_coverage_reads(grouped_reads, start, end, top_n):
        def calculate_coverage(read):
            # Adjust start and end to include clipping
            adjusted_start = read.reference_start - sum(length for op, length in read.cigartuples if op in {4, 5})  # Hard clipping
            adjusted_end = read.reference_end + sum(length for op, length in read.cigartuples if op in {4, 5})  # Soft clipping

            overlap_start = max(adjusted_start, start)
            overlap_end = min(adjusted_end, end)
            return max(0, overlap_end - overlap_start)

        reads_with_coverage = []

        for read in grouped_reads:
            coverage = calculate_coverage(read) / (end - start) * 100
            reads_with_coverage.append({
                "read_name": read.query_name,
                "read_length": read.query_length,
                "coverage": coverage,
                "start": read.reference_start,
                "end": read.reference_end
            })

        # Sort reads by coverage (descending), then by read length (descending)
        sorted_reads = sorted(reads_with_coverage, key=lambda x: (x["coverage"], x["read_length"]), reverse=True)
        return sorted_reads[:top_n]

    # Process BAM file and find the top reads by coverage and length
    grouped_reads = [read for read in bam_before.fetch(contig_b, start_b, end_b)]
    top_reads = get_top_coverage_reads(grouped_reads, start_b, end_b, top_n)

    # Close BAM files
    bam_before.close()

    return top_reads



# # 2015T
# bam_file_before = "/home/yomna/hpc_project/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/mapped_reads_to_scafs_sorted.bam"
# bam_file_after = "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/additional_validation/mapping_2015_22_to_update_reference/mapped_sorted.bam"
# max_coverage_before, max_coverage_after=group_reads_by_cigar(bam_file_before, bam_file_after)
# #2020T
# bam_2020_before="/home/yomna/hpc_project/ToxoVar/analysis/Sniffles/2020T/calls_to_ref_sorted.bam"
# bam_2020_after="/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/additional_validation/mapping_2020T_to_update_reference/mapped_sorted.bam"
# before_groups_2020, after_groups_2020=group_reads_by_cigar(bam_2020_before, bam_2020_after)
#2000B
bam_2000_before="/home/yomna/hpc_project/ToxoVar/analysis/Sniffles/2000B/calls_to_ref_sorted.bam"
before_groups_2000=group_reads_by_cigar(bam_2000_before)
