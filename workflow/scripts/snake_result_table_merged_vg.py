#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:19:23 2024

@author: marie
"""

import csv
import sys

merged_med_snif_vcf=sys.argv[1]
vg_list=sys.argv[2:-1]
resulttable=sys.argv[-1]
def read_merged_vcf(merged_vcf,number_samples):
    merged_vcf_dict={}
    with open (merged_vcf, "r") as infile:
        for line in infile:  #ID column also need to be added to make the check if the variant was called easier
            if not line .startswith("#"):
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                vartype=",".join([x.split("-") [2] for x in splitline[7].split(",")]) #splitline[7].split("-")[2]  doesn't take into accout different types of variants overlapping
                v_list=splitline[9:] #the samples start from column 9                
                #A "." in the merged vcf indicate that there was a problem in getting the genotype for this sample 
                #'due to duplicated lines in concat vcf, to work with numbers instead of strings its converted into a -2 
                if (".") in v_list:
                    v_list = ['-' if element == '.' else element for element in v_list]
                
                v_list = [(element) for element in v_list]
                
                if not (chrom,pos) in merged_vcf_dict:
                    v=[0]*number_samples
                    #row structure: chrom, pos, ref, alt, vartype, valid/notvalid infor from concat vcf, valid/not valid info from merged vcf, valid/not_valid info from vg 
                    merged_vcf_dict[(chrom,pos)]=[[splitline[7],splitline[3], splitline[4], vartype], v_list, v] #do I need this last list?
            else:
                splitline=line.strip().split("\t")
                sample_names=splitline[9:]
    return (merged_vcf_dict,sample_names)

def read_variantcalls_per_sample(vg_infile, sample,dict_pos, combi_dict):
# read in the information of the concatenated files used for merging and variant calls from vg: 
# it is checkt which of the alt alles from the merged vcf file is called for the variant 

    with open (vg_infile, "r") as infile:
        cnt=0
        vg_called_not_in_dic=0
        vg_called_diff_allele=0
        
        print("analyzing: ", infile)
        for line in infile: 
            if not line.startswith("#"):
                cnt=cnt+1
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                alt=splitline[4] #altternative allele from the variant 
                genotype=splitline[9].split(":")[0]
                if (chrom, pos) in combi_dict.keys(): #if  variant is in the merged vcf dictionary 
                    alt_allel_combi= combi_dict[(chrom, pos)][0][2] #get the different alt allele(s) called at this position after merging with Pangenie
                    #print(combi_dict[(chrom, pos)])
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

def write_table(output_dic, outfile,sample_names_m):
    #with open ("outtable.tsv", "w"):
    output_file = outfile

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file,delimiter='\t')    
        # Write header if needed
        writer.writerow(['Chromosom', 'Position','ID', 'Ref', 'ALT', 'Vartype']+sample_names_m+['vg']*len(vg_list))  
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
number_samples=len(vg_list)
merged_dict,sample_names_m=read_merged_vcf(merged_med_snif_vcf,number_samples)
for i in range(0,len(vg_list)):
    result_dict=read_variantcalls_per_sample(vg_list[i], i, 2 ,merged_dict)

write_table(result_dict, resulttable,sample_names_m)


