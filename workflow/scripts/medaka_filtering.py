#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:47:52 2025

@author: yomna
"""

import pandas as pd
import re
from pathlib import Path
import sys
import numpy as np
# Path to files
vcf_file =  sys.argv[1]
pileup_file = sys.argv[2]
file_path = Path(vcf_file)
parent_directory = file_path.parent

def add_mpileup_to_comb_medaka(pileup_file, medaka_dict):
    mpileup_dict = {}
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


def check_deletions_in_snps(Toxo_dict_i,total_dp,var_type,snp_del_overlap):
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
        if perc >= snp_del_overlap:
            remove = True    
    return remove
def add_eval_1 ( Toxo_dict_i,var_type,depth=10,GQ=10):
         qual=int(round(float(medaka_dict[i][3])))
         total_dp=medaka_dict[i][-2] #depth from medaka
         hq_dp=medaka_dict[i][-1] #depth from pileup
         remove = check_deletions_in_snps(Toxo_dict_i,hq_dp,var_type,snp_del_overlap)      
         base_sum=[]
         if var_type =="SNP": 
            for keys in Toxo_dict_i[0]:
                if keys[0].lower() == Alt_base.lower():
                   base_sum.append(Toxo_dict_i[0][keys])  
            allele_freq = sum(base_sum) / hq_dp * 100 if hq_dp > 0 else 0       
         if var_type =="INS":
            for keys in Toxo_dict_i[0]:
                if keys.lower() == Alt_base.lower():
                   base_sum.append(Toxo_dict_i[0][keys])
            allele_freq = sum(base_sum) / hq_dp * 100 if hq_dp > 0 else 0 
         if var_type == "DEL":
            af = []
            for positions, d in zip(Toxo_dict_i, hq_dp):
                base_sum.extend(positions[keys] for keys in positions if keys in ["(", ")"])
                allele_freq = sum(base_sum) / d * 100 if d > 0 else 0
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
            
                       
# Parse the VCF file
with open(vcf_file) as medaka:
    medaka_dict={}
    for line in medaka: # to read the file line by line
        if not (line.startswith("#")):
            arr = line.strip().split("\t")  
            chom = arr[0]
            pos = arr[1]
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
depth= int(sys.argv[3])
GQ= int(sys.argv[4])
snp_del_overlap = float(sys.argv[5])

for i in medaka_dict: 
    Alt_base=medaka_dict[i][2]        
    Toxo_dict_i=medaka_dict[i][8]
    var_type=medaka_dict[i][5]
    add_eval_1(Toxo_dict_i,var_type,depth=depth,GQ=GQ)
    
       
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

