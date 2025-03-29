#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:47:52 2025

@author: yomna
"""

import pandas as pd
import re
from pathlib import Path
#import matplotlib.pyplot as plt
#import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np
# Path to files
vcf_file =  sys.argv[1] #"/home/yomna/hpc_project/ToxoVar/analysis/medaka/medaka_2020T/medaka.annotated_with_VarType.vcf"# 
pileup_file = sys.argv[2]# "/home/yomna/hpc_project/ToxoVar/analysis/pileup/2020T_mpileup.txt"
file_path = Path(vcf_file)
parent_directory = file_path.parent
#valid_vcf_out =  sys.argv[3]#"/home/yomna/hpc_project/ToxoVar/analysis/medaka/medaka_2020T/medaka.annotated_with_VarType_valid.vcf"
# def add_mpileup_to_comb_medaka(pileup_file,medaka_dict):
#     with open(pileup_file) as mpileup_2015:           
#         for line_2 in mpileup_2015: # to read the file line by line
#             arr_mp= line_2.strip().split("\t") # line.strip() is to remove white spaces from the line. split fuction will create the list  
#             chom_mp = arr_mp[0]
#             pos_mp = arr_mp[1]
#             if (chom_mp,pos_mp) in medaka_dict: 
#                if medaka_dict[(chom_mp,pos_mp)] == "DEL":
#                   l_deletions= len(medaka_dict[(chom_mp,pos_mp)][2]) -1  
#                if len(arr_mp) == 4:
#                    depth_mp=arr_mp[3]
#                    Toxo=depth_mp
#                    Toxo_dict_i = {c: int(x) for c, x in (subString.split("=") for subString in arr_mp[3].split(";"))}
#                    medaka_dict[(chom_mp,pos_mp)][-3]=[Toxo,Toxo_dict_i]
#                    sum_i=sum(Toxo_dict_i.values())
#                    medaka_dict[(chom_mp,pos_mp)][-1]=sum_i
#                else: # if all reads were filtered outs
#                    medaka_dict[(chom_mp,pos_mp)][-3]=["",{}]
#                    medaka_dict[(chom_mp,pos_mp)][-1]=0

def add_mpileup_to_comb_medaka(pileup_file, medaka_dict):
    # Create a dictionary to store mpileup information
    mpileup_dict = {}

    # Read the mpileup file and populate the dictionary
    with open(pileup_file) as mpileup_2015:
        for line in mpileup_2015:
            arr_mp = line.strip().split("\t")
            chrom_mp, pos_mp = arr_mp[0], arr_mp[1]
            #print(type(pos_mp))
            if len(arr_mp) >= 4:  # Ensure depth information exists
                depth_mp = arr_mp[3]
                Toxo_dict_i = {c: int(x) for c, x in (subString.split("=") for subString in depth_mp.split(";"))}
                mpileup_dict[(chrom_mp, pos_mp)] = Toxo_dict_i
            else:
                mpileup_dict[(chrom_mp, pos_mp)] = {}
    #print(mpileup_dict)
    # Update medaka_dict based on mpileup_dict
    for (chrom, pos), variant_info in medaka_dict.items():
        if (chrom, pos) in mpileup_dict:
            if variant_info[5] == "DEL":  # If it's a deletion
                deletion_length = len(variant_info[1]) - 1
                deletion_info = []

                # Add all relevant positions for the deletion range
                for offset in range(1,deletion_length + 1):
                    current_pos = str(int(pos) + offset)  # Convert to string to match mpileup_dict keys
                    #print((chrom, current_pos))
                    if (chrom, current_pos) in mpileup_dict:
                        deletion_info.append(mpileup_dict[(chrom, current_pos)])
                    else:
                        deletion_info.append({})  # Add empty dict if no data exists
                
                medaka_dict[(chrom, pos)][-3] = deletion_info
                medaka_dict[(chrom, pos)][-1] = [sum(d.values()) for d in deletion_info]
            else:  # For non-deletion variants
                medaka_dict[(chrom, pos)][-3] = [mpileup_dict[(chrom, pos)]]
                medaka_dict[(chrom, pos)][-1] = sum(mpileup_dict[(chrom, pos)].values())


def check_deletions_in_snps(Toxo_dict_i,total_dp,var_type):
    remove = False
    if var_type == "SNP":
        s_pos=0
        s_neg= 0
        if "(" in Toxo_dict_i[0]:
            s_pos= Toxo_dict_i[0]["("]
        if ")" in Toxo_dict_i[0]:
            s_neg= Toxo_dict_i[0][")"]
        total = s_pos+s_neg
        perc = total /total_dp
        if perc >= 0.5:
            remove = True    
    return remove
def add_eval_1 ( Toxo_dict_i,var_type,depth=10,AF=50,GQ=10):
         qual=int(round(float(medaka_dict[i][3])))
         total_dp=medaka_dict[i][-2] #depth from medaka
         hq_dp=medaka_dict[i][-1] #depth from pileup
         #if Alt_base.lower() in Toxo_i:
         remove = check_deletions_in_snps(Toxo_dict_i,hq_dp,var_type)      
        #print("sum_i: ", sum_i)
         base_sum=[]
         if var_type =="SNP": 
    #      print("SNP: keys, altbase", keys[0].lower(), ", ", Alt_base[b].lower() )
            for keys in Toxo_dict_i[0]:
                if keys[0].lower() == Alt_base.lower():
                   base_sum.append(Toxo_dict_i[0][keys])  
            allele_freq = sum(base_sum) / hq_dp * 100 if hq_dp > 0 else 0       
         if var_type =="INS":
     #     print("INS: keys, altbase", keys.lower(), ", ", Alt_base[b].lower() )
            for keys in Toxo_dict_i[0]:
                if keys.lower() == Alt_base.lower():
                   base_sum.append(Toxo_dict_i[0][keys])
            allele_freq = sum(base_sum) / hq_dp * 100 if hq_dp > 0 else 0 
         if var_type == "DEL":
            af = []
            for positions, depth in zip(Toxo_dict_i, hq_dp):
                base_sum.extend(
                    positions[keys] for keys in positions if keys in ["(", ")"]
                )
                # Calculate allele frequency
                allele_freq = sum(base_sum) / depth * 100 if depth > 0 else 0
                af.append(allele_freq)
            allele_freq = min(af) if af else 0
         if var_type == "MNP" or  var_type == "MIXED":
            allele_freq = np.nan
         if total_dp >= depth and qual >= GQ and (not remove):                
             medaka_dict[i].append(1)
             medaka_dict[i][8].extend([allele_freq])
         else: 
             medaka_dict[i].append(0)
             medaka_dict[i][8].extend([allele_freq])
         if remove:
            _, ref, alt, q, _, vartype, _,_, pileup,dp,depth_hq,valid = medaka_dict[i]
            allele_fre= pileup[-1]
            p = pileup[0]
            snps_overlaping_with_deleltions.append({"CHROM": i[0], "POS": int(i[1]), "DP_total": dp, "DP_hq": depth_hq,"GQ": float(q), "Variant_Type": vartype,"Valid": valid,"AF":allele_fre,"pileup":p,"alt":alt})               


# #old way of making valid and invalid    
# def add_eval (Toxo_i, Toxo_dict_i,var_type):
#          qual=float(medaka_dict[i][3])
#          if Alt_base.lower() in Toxo_i:
#             sum_i=medaka_dict[i][-1] #This will decide which 
#             #print("sum_i: ", sum_i)
#             base_sum=[]
#             for keys in Toxo_dict_i:
#                 if var_type =="SNP": 
#               #      print("SNP: keys, altbase", keys[0].lower(), ", ", Alt_base[b].lower() )
#                     if keys[0].lower() == Alt_base.lower():
#                         base_sum.append(Toxo_dict_i[keys])
                        
#                 if var_type =="INS":
#                #     print("INS: keys, altbase", keys.lower(), ", ", Alt_base[b].lower() )
#                     if keys.lower() == Alt_base.lower():
#                        base_sum.append(Toxo_dict_i[keys])
                       
#                 if var_type =="DEL":
#                 #    print("DEL:")
#                     alt_base_pattern = re.escape(Alt_base.lower()) + '-+'
#                     alt_base_regex = re.compile(alt_base_pattern, re.IGNORECASE)
#                     if alt_base_regex.match(keys.lower()):
#                  #       print("DEL: keys, altbase","match:" , ", ", keys.lower(), ", ", Alt_base[b].lower())
#                         base_sum.append(Toxo_dict_i[keys])
#                     #else: 
#                         #print("no match")
#                 #print("base_sum", base_sum)
#             base_sum_c=sum(base_sum)/sum_i*100      
#             if sum_i >= depth and base_sum_c >= allele_freq and qual >= 1:                
#                 medaka_dict[i].append(1)
#                 medaka_dict[i][8].extend([base_sum_c])
#             else: 
#                 medaka_dict[i].append(0)
#                 medaka_dict[i][8].extend([base_sum_c])
#          else:
#              medaka_dict[i].append(0)
#              medaka_dict[i][8].extend([0])

def extract_info_fields(info_str):
    # Extract DP (depth) from INFO
    dp_match = re.search(r"DP=(\d+)", info_str)
    dp = int(dp_match.group(1)) if dp_match else None
    # Extract Variant Type from INFO
    vartype_match = re.search(r"VARTYPE=([A-Z]+)", info_str)
    vartype = vartype_match.group(1) if vartype_match else "Unknown"
    return dp, vartype

def extract_valid(medaka_infile, valid_vcf_out, valid_entries):
    with open(medaka_infile, 'r') as infile, open(valid_vcf_out, 'w') as outfile:
        used_valid=[]
        for line_3 in infile:
            if line_3.startswith('#'):
                outfile.write(line_3)  # Write header lines to the output file
            else:
                columns = line_3.strip().split('\t')
                chromosom = columns[0]
                position =  columns[1]
                alt=columns[4]
                qual=columns[5]
                for v in valid_entries:
                    if (v[0]==chromosom and v[1]==position and v[2]== alt and v[3]==qual):
                        if (v[0],v[1]) not in used_valid:
                            outfile.write(line_3)  # Write matching lines to the output file
                            used_valid.append((v[0],v[1]))
            
                       
# Step 1: Parse the VCF file into a DataFrame and create a set of relevant positions
# Parse the VCF file
with open(vcf_file) as medaka:
    medaka_dict={}
    for line in medaka: # to read the file line by line
        if not (line.startswith("#")):
            arr = line.strip().split("\t") # line.strip() is to remove white spaces from the line. split fuction will create the list  
            chom = arr[0]
            pos = arr[1]
            #vartype=arr[7].split("=")[-1]
            depth_total,vartype =extract_info_fields(arr[7])
            arr[7] = vartype
            if not (chom,pos) in medaka_dict: # it's important to include this if statement because medaka_comp has repeated lines
                medaka_dict[(chom,pos)]=arr[2:]    
                medaka_dict[(chom,pos)].append([{}])
                medaka_dict[(chom,pos)].append(depth_total) #add total depth
                medaka_dict[(chom,pos)].append(0)# this will be a placeholder for high quality depth
            elif float(medaka_dict[(chom,pos)][3]) < float(arr[5]): 
                 #print(arr)
                 medaka_dict[(chom,pos)]=arr[2:]
                 medaka_dict[(chom,pos)].append([{}]) # this will contain pileup information
                 medaka_dict[(chom,pos)].append(depth_total) #add total depth
                 medaka_dict[(chom,pos)].append(0)  # this will contain pileup depth_hq from pileup

add_mpileup_to_comb_medaka (pileup_file,medaka_dict)                 
depth= int(sys.argv[3])#10
AF= int(sys.argv[4])#50
GQ= int(sys.argv[5])#10
snps_overlaping_with_deleltions=[]

for i in medaka_dict: 
   # print("entry:", i)
    Alt_base=medaka_dict[i][2]        
    Toxo_dict_i=medaka_dict[i][8]
    var_type=medaka_dict[i][5]
    #add_eval(Toxo_i, Toxo_dict_i,var_type)
    add_eval_1(Toxo_dict_i,var_type,depth=depth,AF=AF,GQ=GQ)
    
vcf_data = []

# Parse the VCF file
for i in medaka_dict:
    _, ref, alt, q, _, vartype, _,_, pileup,dp,depth_hq,valid = medaka_dict[i]
    allele_fre= pileup[-1]
    p = pileup[:-1]
    vcf_data.append({"CHROM": i[0], "POS": int(i[1]), "DP_total": dp, "DP_hq": depth_hq,"GQ": float(q), "Variant_Type": vartype,"Valid": valid,"AF":allele_fre,"pileup":p,"ref":ref,"alt":alt})

df = pd.DataFrame(vcf_data)
df.to_csv(parent_directory/"table_with_variants_info.csv",index=False)
#df=pd.read_csv("/home/yomna/hpc_project/ToxoVar/analysis/medaka/medaka_2000B/table_to_make_medaka_plot.csv")
#df=pd.read_csv("/home/yomna/hpc_project/ToxoVar/analysis/medaka/medaka_2020T/table_to_make_medaka_plot.csv")
df_snp = pd.DataFrame(snps_overlaping_with_deleltions)
df_snp.to_csv(parent_directory/"table_snp_in_deletions_info.csv")


    
    
variant_types = df["Variant_Type"].unique()

# Get global ranges for consistent axes
# Get global ranges for consistent axes
# dp_min, dp_max = df_filtered["DP_total"].min(), df_filtered["DP_total"].max()
# gq_min, gq_max = df_filtered["GQ"].min(), df_filtered["GQ"].max()
# af_min, af_max = df_filtered["AF"].min(), df_filtered["AF"].max()

# for vartype in variant_types:
#     subset = df_filtered[df_filtered["Variant_Type"] == vartype]
    
#     # Create a 3D plot
#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # Plot data points
#     for _, row in subset.iterrows():
#         color = "red" if row["Valid"] == 1 else "blue"
#         ax.scatter(row["DP_total"], row["GQ"], row["AF"], color=color, s=50, edgecolor="k", alpha=0.8)
    
#     # Set axis labels
#     ax.set_title(f"3D Plot of Allele Frequency vs Genotype Quality vs Depth ({vartype})")
#     ax.set_xlabel("Depth (DP)")
#     ax.set_ylabel("Genotype Quality (GQ)")
#     ax.set_zlabel("Allele Frequency (AF)")
    
#     # Set axis limits
#     ax.set_xlim(dp_min, dp_max)
#     ax.set_ylim(gq_min, gq_max)
#     ax.set_zlim(af_min, af_max)
    
#     # Show plot
#     plt.tight_layout()
#     #plt.savefig("3d.pdf")
#     plt.show()
    
# import plotly.express as px
# import plotly.io as pio

# # Set renderer
# pio.renderers.default = "browser"  # Change to "notebook" if using Jupyter

# # Loop through each Variant_Type
# for vartype in df_filtered["Variant_Type"].unique():
#     subset = df_filtered[df_filtered["Variant_Type"] == vartype]
    
#     # Check if the subset is empty
#     if subset.empty:
#         print(f"No data for Variant_Type: {vartype}")
#         continue

#     print(f"Creating plot for Variant_Type: {vartype}")

#     fig = px.scatter_3d(
#         subset,
#         x="DP_total",        # Depth (DP) on x-axis
#         y="GQ",              # Genotype Quality (GQ) on y-axis
#         z="AF",              # Allele Frequency (AF) on z-axis
#         color=subset["Valid"].replace({1: "Valid", 0: "Invalid"}),  # Map Validity to labels
#         color_discrete_map={"Valid": "red", "Invalid": "blue"},     # Define custom colors
#         title=f"3D Plot of Allele Frequency vs Genotype Quality vs Depth ({vartype})",
#         labels={"DP_total": "Depth (DP)", "GQ": "Genotype Quality (GQ)", "AF": "Allele Frequency (AF)"},
#         size=[5 for _ in range(len(subset))],  # Set dot size smaller
#     )

#     # Display the plot
#     fig.show()

# Round Depth (DP_total) to the nearest integer

df_filtered = df #[(df["DP_total"] <= 50) & (df["GQ"] <= 500)]
df_filtered["GQ"] = df_filtered["GQ"].round().astype(int)
df_filtered.loc[df_filtered["GQ"] >= 30, "GQ"] = 30

for vartype in variant_types:
    #vartype="INS"
    subset = df_filtered[df_filtered["Variant_Type"] == vartype]
    # Create a pivot table
    table=pd.pivot_table(
        subset,
        values="CHROM",  # Any column can be used; it's just counting rows
        index="GQ",       # Rows: GQ
        columns="DP_total",  # Columns: Depth (rounded)
        aggfunc="count",  # Count the number of occurrences
        fill_value=0      # Fill missing combinations with 0
    )
    
    table.to_csv( parent_directory / f"GQ_30_depth_total_{vartype}.csv")    
df_filtered_1 = df_filtered[(df_filtered["DP_hq"] !=0)]
for vartype in variant_types:
    subset = df_filtered_1[df_filtered_1["Variant_Type"] == vartype]
    subset["AF_50"] = subset["AF"] < 50  # True if AF > 50%    
    # Pivot table for proportions
    table = pd.pivot_table(
        subset,
        values="AF_50",      # Use the AF_50 column
        index="GQ",          # Rows: GQ
        columns="DP_total",  # Columns: Depth (rounded)
        aggfunc="mean",      # Calculate the proportion (mean of True/False)
        fill_value='NA'       # Fill missing combinations with 0
    )
    table.to_csv(parent_directory /f"GQ_30_depth_total_{vartype}_AF.csv")    
    
for vartype in variant_types:
    subset = df_filtered_1[df_filtered_1["Variant_Type"] == vartype]
    subset["AF_0"] = subset["AF"] == 0  # True if AF == 0  
    # Pivot table for proportions
    table = pd.pivot_table(
        subset,
        values="AF_0",      # Use the AF_50 column
        index="GQ",          # Rows: GQ
        columns="DP_total",  # Columns: Depth (rounded)
        aggfunc="mean",      # Calculate the proportion (mean of True/False)
        fill_value='NA'       # Fill missing combinations with 0
    )
    table.to_csv(parent_directory / f"GQ_30_depth_total_{vartype}_AF_0.csv")
    
valid_variations = []    
for e in medaka_dict:  
    if medaka_dict[e][-1] == 1:
        c=e[0]
        p=e[1]
        alt=medaka_dict[e][2]
        qual=medaka_dict[e][3]
        valid_variations.append((c,p,alt,qual))
valid_vcf_out=  parent_directory / "medaka.annotated_with_VarType_valid.vcf"
extract_valid(vcf_file, valid_vcf_out, valid_variations)     

alleles_that_doesnt_exist = df[df["AF"] == 0]
alleles_that_doesnt_exist.to_csv(parent_directory / "variants_that_doesnt_exist_in_pileup.csv",index=False)
