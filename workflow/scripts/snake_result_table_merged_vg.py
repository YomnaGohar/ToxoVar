#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:19:23 2024

@author: marie
"""

import csv
import sys

merged_med_snif_vcf=sys.argv[1]
vg_calls_2015T=sys.argv[2]
vg_calls_2020T=sys.argv[3]
vg_calls_2000B=sys.argv[4]

vg_list=[vg_calls_2015T,vg_calls_2020T,vg_calls_2000B]
resulttable=sys.argv[5]

def read_merged_vcf(merged_vcf):
    merged_vcf_dict={}
    with open (merged_vcf, "r") as infile:
        for line in infile: 
            if not line .startswith("#"):
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                vartype=splitline[7].split("-")[2]
                v_list=splitline[9:12]
                
                #A "." in the merged vcf indicate that there was a problem in getting the genotype for this sample 
                #'due to duplicated lines in concat vcf, to work with numbers instead of strings its converted into a -2 
                if (".") in v_list:
                    v_list = ['-' if element == '.' else element for element in v_list]
                
                v_list = [(element) for element in v_list]
                
                if not (chrom,pos) in merged_vcf_dict:
                    #row structure: chrom, pos, ref, alt, vartype, valid/notvalid infor from concat vcf, valid/not valid info from merged vcf, valid/not_valid info from vg 
                    merged_vcf_dict[(chrom,pos)]=[[splitline[3], splitline[4], vartype], v_list, [0,0,0]] #do I need this last list?
                    
    return merged_vcf_dict

def read_variantcalls_per_sample(vg_infile, sample,dict_pos, combi_dict):
# read in the information of the concatenated files used for merging and variant calls from vg: 
# it is checkt which of the alt alles from the merged vcf file is called for the variant 

    with open (vg_infile, "r") as infile:
        cnt=0
        vg_called_not_in_dic=0
        vg_called_diff_allele=0
        
        print("analyzing: ", infile)
        for line in infile: 
            if not line .startswith("#"):
                cnt=cnt+1
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                alt=splitline[4] #altternative allele from the variant 
                genotype=splitline[9].split(":")[0]
                if (chrom, pos) in combi_dict.keys(): #if  variant is in the merged vcf dictionary 
                    alt_allel_combi= combi_dict[(chrom, pos)][0][1] #get the different alt allele(s) called at this position after merging with Pangenie
                    if "," in alt_allel_combi: # if there is more than one alt allele, get a list of those.
                        alt_list=alt_allel_combi.split(",")
                        #check which alt allel is called for the variant: if 1 than first, if two = second, 
                        #if 3= third, if the alt allele in concat vcf is different from those after merging, -1 is printed 
                        position = [i + 1 for i, x in enumerate(alt_list) if x == alt] #check if ALT allele "alt" is in the list
                        if position:
                          #  print("Concat ALT allele:", alt, " is identical to position" , position, "in", alt_list)
                            position=int(position[0]) 
                        else:
                            position="-"
                            print("alt_allele", alt, "not found in", alt_list, "Variant: ",  chrom, pos)
                            vg_called_diff_allele += 1
                    elif alt == alt_allel_combi: #if one alt allele at this variant: identical: set position to 1, if different set position to -1
                        position=1
                    else: 
                        print("different ALT allele", alt," in concat compared to merged", alt_allel_combi, "Variant:", chrom, pos)
                        position="-"
                        vg_called_diff_allele += 1
                    #print(position)
                    if(genotype==".") or (genotype=="0"):
                        position=genotype 
                    combi_dict[(chrom, pos)][dict_pos][sample]=position
                else: 
                    print("variant not in dictionary: ", chrom, pos, sample, "\n", alt)
                    vg_called_not_in_dic += 1 
    
    print(" calls that are not in merged dictionary:", vg_called_not_in_dic)
    print(" calls with different alt allele:", vg_called_diff_allele)
    print("total entries in", infile, ": ", cnt)                
    return combi_dict            

def write_table(output_dic, outfile):
    #with open ("outtable.tsv", "w"):
    output_file = outfile

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
    
        # Write header if needed
        writer.writerow(['Chromosom', 'Position', 'Ref', 'ALT', 'Vartype', '2015T_merged', '2020T_merged','2000B_merged', '2015T_vg', '2020T_vg', '2000B_vg', "gene_id" ,'gene_info'])  
        # Adjust column names as needed
            
        # Write each entry as a comma-separated line
        for key, values in output_dic.items():
            # Flatten the nested lists
            #print(values)
            flattened_values = [item for sublist in values for item in sublist]
            #print(flattened_values)
            writer.writerow([key[0], key[1]] + flattened_values)

##MAIN:
result_dict={}
merged_dict=read_merged_vcf(merged_med_snif_vcf)
for i in range(0,len(vg_list)):
    result_dict=read_variantcalls_per_sample(vg_list[i], i, 2 ,merged_dict)

write_table(result_dict, resulttable)


