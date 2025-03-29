#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 19:36:40 2025

@author: yomna
"""

import pysam

def extract_reads_with_genes(bam_file, gene1, gene2, output_fasta):
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Dictionary to store read names mapping to each gene
    gene_reads = {gene1: set(), gene2: set()}
    
    # Iterate through all reads in the BAM file
    for read in bam.fetch():
        if gene1 in read.query_name:
            gene_reads[gene1].add(read.reference_name)
        if gene2 in read.query_name:
            gene_reads[gene2].add(read.reference_name)
    
    # Find reads that map to both genes
    common_reads = gene_reads[gene1].intersection(gene_reads[gene2])
    
    # Write the reads to a FASTA file
    with open(output_fasta, "w") as fasta:
        for read in bam.fetch():
            if read.query_name in common_reads:
                fasta.write(f">{read.query_name}\n{read.query_sequence}\n")
    
    bam.close()

# Example usage
bam_file = "/home/yomna/hpc_project/ToxoVar/analysis/Sniffles/2000B/ROP8_ROP2A/alignment_of_genes_to_contig_46_7298291_7327278_sorted.bam"  # Path to your BAM file
gene1 = "TGME49_215750.R866"
gene2 = "TGME49_275310.R2947"
output_fasta = "common_reads.fasta"

extract_reads_with_genes(bam_file, gene1, gene2, output_fasta)
