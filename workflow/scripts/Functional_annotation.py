#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:51:29 2025

@author: yomna
"""

import pandas as pd
import sys
import gzip
import matplotlib.pyplot as plt
# Define file paths
vep_file = "/home/yomna/hpc_project/ToxoVar/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vep"
vcf_file = "/home/yomna/hpc_project/ToxoVar/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized.vcf" #sys.argv[2]
# File path for the uploaded GFF3 file
gff_file_path = "/home/yomna/hpc_project/ToxoME49_2015_2022/combined_with_pipline_after_rebasecalling_and_herro_using_all_read_that_donot_map_to_human_or_mouse/after_mit_removal/assembly/nano-corr/toxoplasma-gondii-annotation_results/scaffold.out.modified.with.gffread.for.vep.gff3.gz"

# Parse GFF file to extract gene descriptions
gff_data = []
with gzip.open(gff_file_path, 'rt') as gff:
    for line in gff:
        if line.startswith("#"):
            continue  # Skip comments
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue  # Skip malformed lines

        chrom, source, feature_type, start, end, score, strand, phase, attributes = cols
        if feature_type != "gene":
            # Extract gene ID and description from attributes field
            attributes_dict = {attr.split("=")[0]: attr.split("=")[1] for attr in attributes.split(";") if "=" in attr}
            gene_id = attributes_dict.get("Derives_from", None)
            description = attributes_dict.get("product", None)
    
            if gene_id and description:
                gff_data.append([gene_id, description])

# Convert GFF data into a DataFrame
gff_df = pd.DataFrame(gff_data, columns=["Feature", "Description"])
# Read VCF file and extract relevant information
vcf_data = []
with open(vcf_file, 'r') as vcf:
    for line in vcf:
        if line.startswith("#"):  # Skip header lines
            continue
        cols = line.strip().split('\t')
        chrom, location = cols[:2]
        ID_before=cols[2]
        ID= cols[7]
        genotype_info = cols[11:]  # Assuming 2000B and 2020T are in these columns
        vcf_data.append([chrom, location,ID_before, ID] + genotype_info)

vcf_df = pd.DataFrame(vcf_data, columns=["Chrom", "Location","ID_before_Normalization", "ID", "2020T","2000B"])

# Read VEP file and extract relevant information
vep_data = []
with open(vep_file, 'r') as vep:
    for line in vep:
        if line.startswith("#"):  # Skip header lines
            continue
        cols = line.strip().split('\t')
        ID = cols[0]
        start_end =  cols[1]
        allele, gene, feature, feature_type, consequence = cols[2:7]
        cDNA_position, CDS_position, protein_position = cols[7:10]
        amino_acids, codons = cols[10:12]
        vep_data.append([ ID, start_end,allele, gene, feature, feature_type, consequence, 
                         cDNA_position, CDS_position, protein_position, amino_acids, codons])

vep_df = pd.DataFrame(vep_data, columns=[ "ID_before_Normalization","start_end", "Allele", "Gene", "Feature", "Feature_type", 
                                          "Consequence", "cDNA_position", "CDS_position", "Protein_position", 
                                          "Amino_acids", "Codons"])

# Merge VEP and VCF data on Chrom, Location, and ID
merged_df = pd.merge(vep_df, vcf_df, on=[ "ID_before_Normalization"], how="left")
merged_df=pd.merge(merged_df, gff_df, on=[ "Feature"], how="left")
# Display the final table
table_path =  "/home/yomna/hpc_project/ToxoVar/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized_vep_table.csv" #sys.argv[3]
merged_df.to_csv(table_path, index=False)
#remove large variants
#contig_2-2235314-DEL-0-520
#contig_35-699698-DEL-0-203
#contig_35-699715-DEL-0-202
#contig_35-1340293-DEL-0-1180
#contig_42-608619-DEL-0-52
#contig_46-7310767-INS-0-12896
filtered_df = merged_df[~merged_df["ID_before_Normalization"].isin([
    "contig_2-2235314-DEL-0-520", "contig_35-699698-DEL-0-203", 
    "contig_35-699715-DEL-0-202", "contig_35-1340293-DEL-0-1180", 
    "contig_42-608619-DEL-0-52", "contig_46-7310767-INS-0-12896"
])]
filtered_df["2020T"] = pd.to_numeric(filtered_df["2020T"], errors="coerce")
filtered_df["2000B"] = pd.to_numeric(filtered_df["2000B"], errors="coerce")

effect_counts_2020T = filtered_df[filtered_df["2020T"] != 0]["Consequence"].value_counts()
effect_counts_2000B = filtered_df[filtered_df["2000B"] != 0]["Consequence"].value_counts()

# Combine both datasets into a single DataFrame for plotting
effect_counts_df = pd.DataFrame({"2020T": effect_counts_2020T, "2000B": effect_counts_2000B}).fillna(0)

# Plot the bar chart
plt.figure(figsize=(12, 6))
effect_counts_df.plot(kind="bar", stacked=False, width=0.8)
plt.xlabel("Effect Type (Consequence)")
plt.ylabel("Count")
plt.title("Variant Effects in 2020T and 2000B")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Sample")
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.tight_layout()
plt.savefig("/home/yomna/hpc_project/ToxoVar/analysis/Graph/graph_construction/results/with_2015T_filtering/merged_vg_combined_table_placed_ref_for_igv_no_2015_normalized_vep_table_without_SV.pdf" )
plt.show()

