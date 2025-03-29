#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 23:32:44 2025

@author: yomna
"""
import sys

def parse_vcf_sniffles(file_path):
    """Parse a VCF file and return a list of variants as tuples."""
    variants = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip header lines
            cols = line.strip().split("\t")
            chrom = cols[0]       # Chromosome/contig name
            pos = int(cols[1])    # Start position
            info = cols[7]        # INFO field
            var_type = None       # Variant type (e.g., DEL, INS, INV, etc.)            
            # Extract END position and SVTYPE from the INFO field
            end = pos  # Default for single-point variants
            for entry in info.split(";"):
                if entry.startswith("END="):
                    end = int(entry.split("=")[1])
                if entry.startswith("SVTYPE="):
                    var_type = entry.split("=")[1]
            
            # Append the parsed variant
            variants.append((chrom, pos, end, var_type, line.strip()))
    
    return variants

def parse_vcf_medaka(file_path):
    """Parse a VCF file and return the header and a list of variants as tuples."""
    variants = []
    header = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.strip())  # Store header lines
                continue           
            cols = line.strip().split("\t")
            chrom = cols[0]       # Chromosome/contig name
            pos = int(cols[1])    # Start position
            info = cols[7]        # INFO field
            var_type = info.split("-")[-3]  # Extract variant type            
            if var_type == "DEL":
                end = pos + int(info.split("-")[-1])  # Calculate end position
            else:
                end = pos + 1           
            variants.append((chrom, pos, end, var_type, line.strip()))    
    return header, variants

def parse_bed_file(bed_file):
    """Parse a BED file and return a list of regions to exclude."""
    exclude_regions = []
    with open(bed_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            chrom = cols[0]       # Chromosome/contig name
            start = int(cols[1])  # Start position
            end = int(cols[2])    # End position
            exclude_regions.append((chrom, start, end))
    return exclude_regions

def filter_vcf(file1, file2, bed_file, output_file, log_file):
    """Filter variants based on overlap and a specified region to exclude."""
    # Parse the VCF files
    header, vcf1 = parse_vcf_medaka(file1)
    vcf2 = parse_vcf_sniffles(file2)
    exclude_regions = parse_bed_file(bed_file)
    
    filtered_variants = []
    for chrom1, pos1, end1, var_type1, line1 in vcf1:
        # Check for overlap with file2 variants
        overlap = False
        for chrom2, pos2, end2, var_type2, _ in vcf2:
            if chrom1 == chrom2 and not (end1 < pos2 or end2 < pos1):
                overlap = True
                with open(log_file, "a") as log:
                    log.write(f'Overlap in {var_type1} with Sniffles detected and it will be removed: {chrom1}, {pos1}-{end1} overlaps with {chrom2}, {pos2}-{end2}\n')
                break
        
        # Check if the variant falls within the excluded regions
        for exclude_chrom, exclude_start, exclude_end in exclude_regions:
            if chrom1 == exclude_chrom and not (end1 < exclude_start or pos1 > exclude_end):
                overlap = True
                with open(log_file, "a") as log:
                    log.write(f'Variant on excluded region detected and removed: {chrom1}, {pos1}-{end1}, Type: {var_type1}\n')
                break
        
        # Keep the variant if there's no overlap or exclusion
        if not overlap:
            filtered_variants.append(line1)
    
    # Write the output VCF with the header
    with open(output_file, "w") as out:
        for line in header:  # Write header lines
            out.write(line + "\n")
        for variant in filtered_variants:  # Write filtered variants
            out.write(variant + "\n")

# Example usage
file1 = sys.argv[1]  # Medaka VCF file
file2 = sys.argv[2]  # Sniffles VCF file
bed_file = sys.argv[3]  # BED file for excluded regions
log_file = sys.argv[4]  # Log file for excluded variants
output = sys.argv[5]  # Output filtered VCF file

filter_vcf(file1, file2, bed_file, output, log_file)
