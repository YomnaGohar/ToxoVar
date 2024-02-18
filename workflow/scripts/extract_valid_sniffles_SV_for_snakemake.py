#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Jun  2 09:56:12 2023

@author: Marie
"""
import pandas as pd
import sys
#Parser for validating SV and counting the different types of SV within the sniffles results: 
##Extract information from sniffles vcf file##

def check_valid(inline):
    valid=False
    if inline[6]=="GT": #Exclude variant calls with GT in Genotype column 
        #print("In Line: ", inline, "Filter is =", inline[6])
        return (valid)
    if (inline[9].split(":")[0]=="0/0"): #Exclude variant calls with 0/0 calls 
        #print("In line: ", inline, "Genotype is =", inline[9])
        return valid
    if (inline[9].split(":")[0]=="./."): #Exclude variant calls with undefined Genotype call 
        #print("In line: ", inline, "Genotype is =", inline[9])
        return valid
    #print(inline, "is valid")
    valid=True
    return valid
    

def make_dir(infile, valid_sniffles_dic,allSV_dict): #Make directory from the sniffles data
    with open (infile) as sniffles:
        for line in sniffles: 
            if not line.startswith("#"):
               split_line=line.strip().split("\t")
               split_line[7]=split_line[7].strip().split(";")
               chrom=split_line[0]
               pos=split_line[1]
               key=(chrom, pos)
               if not key in allSV_dict:
                   allSV_dict[key]=[]
                   allSV_dict[key].append(split_line[2:])
               else:
                    del allSV_dict[key]               
               if check_valid(split_line) == True:
                   if not key in valid_sniffles_dic:
                       valid_sniffles_dic[key]=[]  
                       valid_sniffles_dic[key].append(split_line[2:])
                   else:
                        del valid_sniffles_dic[key]
        return valid_sniffles_dic, allSV_dict
    

def extract_sv_type(sublist, sv_list):   #Get SV type
    if len(sublist) >0:
        sv=sublist[5][1].split("=")[1]
        sv_list[0].append(sv)
        return sv_list
 
def getSV_list(input_dict, outlist): #
    for key in input_dict.keys():
        sublist_i=input_dict[key][0]
        outlist=extract_sv_type(sublist_i, outlist)
    return outlist

###extract valid SV###
def extract_valid(sniffles_infile, valid_vcf_out, valid_entries):    
    with open(sniffles_infile, 'r') as infile, open(valid_vcf_out, 'w') as outfile:
        for line_3 in infile:
            if line_3.startswith('#'):
                outfile.write(line_3)  # Write header lines to the output file
            else:
                columns = line_3.strip().split('\t')
                chromosom = columns[0]
                position =  columns[1]
                for v in valid_entries:
                    if (v[0]==chromosom and v[1]==position):
                        outfile.write(line_3) 


###Sniffles infiles ###

sniffles_2015T=sys.argv[1]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/snakemake_Yomna/analysis/Sniffles/2015T/sniffles_2015T_with_reference.vcf"

file_list = [sniffles_2015T]

###outfiles####

valid_vcf_2015=sys.argv[2]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/snakemake_Yomna/analysis/Sniffles/2015T/sniffles2-2_2015T_with_reference_valid.vcf"

outfile_list=[valid_vcf_2015]

statistics=sys.argv[3]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/snakemake_Yomna/analysis/Sniffles/2015T/statistics.csv"
####MAIN#####

sniffles_dic = {}
all_SV_dict={}  # Dictionary to store the combined lines
for file in file_list:
    sniffles_dic, all_SV_dict=make_dir(file, sniffles_dic, all_SV_dict)
     
valid_2015=[]
    
for e in sniffles_dic.keys():
    valid_2015.append(e)


valid_list=[valid_2015]

for i in range(0,len(file_list)):                     
    extract_valid(file_list[i], outfile_list[0], valid_list[0])

###COUNTING###
all_sv_list=[[]]
valid_SV_list=[[]]


all_SV_types=getSV_list(all_SV_dict, all_sv_list)
valid_SV_types=getSV_list(sniffles_dic, valid_SV_list)

sv_list=[all_SV_types, valid_SV_types]
sv_filt=["All","Valid"]
count_dic={}  
for i,f in zip(sv_list,sv_filt):
    subl=i[0]
    element_count={}
    for element in subl:
        if element in element_count:
            element_count[element] +=1
        else: 
            element_count[element] =1
    count_dic[f]=element_count
table=pd.DataFrame(count_dic).transpose()
table.to_csv(statistics)


