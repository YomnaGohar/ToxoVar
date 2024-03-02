#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 10:41:49 2024

@author: marie
"""
##Plots Venn-Diagramm for each of the files, merged vcf and vg calls 
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import csv
import sys

csv_file_path = sys.argv[1]

outfile={"merged":sys.argv[2], "vg": sys.argv[3]}


def get_vt(REF,ALT,VT):
    if (VT == "INS"):
        if len(ALT)-len(REF) < 50:
            VT="sINS"
    elif (VT == "DEL"):
        if len(REF)-len(ALT) < 50:
            VT = "sDEL"
    elif (VT == "SNV"):
        VT = "SNP"
    return VT

def read_csv(csv_file_path):
    per_file_vcf={"merged":{"SNP": {}, "sINS":{}, "INS":{}, "sDEL":{}, "DEL":{}, "ALL":{}}, 
                    "vg":{"SNP": {}, "sINS":{}, "INS":{}, "sDEL":{}, "DEL":{}, "ALL":{}}}
    
    # Open the CSV file and read its content
    with open(csv_file_path, 'r') as file:
        # Create a CSV reader object
        csv_reader = csv.reader(file)
        next(csv_reader)
    
        # Iterate through each row in the CSV file
        for row in csv_reader:
           # print(row)
            # Process each row as needed
            Vtype=get_vt(row[2], row[3],row[4])
            
            # per_file
            #code_concat="".join(row[5:8])
            code_merged="".join(row[5:8])
            code_vg="".join(row[8:11])
            
            for key in per_file_vcf.keys():
                if key=="merged":
                    code = code_merged
                elif key == "vg":
                    code = code_vg
                
                if code not in per_file_vcf[key]["ALL"].keys():
                    per_file_vcf[key]["ALL"][code]=1
                else:
                    per_file_vcf[key]["ALL"][code]+=1
    
                if code not in per_file_vcf[key][Vtype].keys():
                    per_file_vcf[key][Vtype][code]=1
                else:
                    per_file_vcf[key][Vtype][code]+=1
        
          
        return(per_file_vcf)
        
def get_catogorized_dict(vcf_dictionary):
    
    total_categorized_dict={}
    for key in vcf_dictionary:
        print("new file")
        total_categorized_dict[key]=categorize_calls_3(vcf_dictionary[key])
    return total_categorized_dict

def contains_only_0_and_1(input_str):
        set_code=set(input_str)
        if len(set_code) >= 2:
            if(len(set_code) == 2 and "0" in set_code):
                return False
            else: 
                return True
            return True
        else:
            return False
    
def modify_code(code):
    #code with "-" uncalled gt for the sample, set to 0 counted like not called for this sample
    char_code=""
    for char in code:
        if char == '-':
            char_code += '0'
        else: 
            char_code += char
    return char_code
  
def categorize_calls_3(file_dict):
    categorized_vcf = {"SNP":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0}, 
                       "sINS":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0}, 
                       "sDEL":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0}, 
                       "INS":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0},
                       "DEL":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0},
                       "ALL":{'a00': 0, '0a0': 0, '00a': 0, 'aa0': 0, 'a0a': 0, '0aa': 0, 'aaa': 0, "other":0}}

    for vartype in file_dict.keys():
        for code, count in file_dict[vartype].items():
            non_01_code=contains_only_0_and_1(code) #code with - in it (sniffles not able to determine gt for sample) or variants with different alt alleles called for differnet samples
            if non_01_code == True: 
                code=modify_code(code) #- in code is set to 0
            
            # Check if the code matches any of the desired patterns
            if code[0] != '0' and code[1]=='0' and code[2]=='0': #Code 100,200,300
                categorized_vcf[vartype]['a00'] += count
            elif code[0] == '0' and code[1] != '0' and code[2] == '0': #Code 010 020 030
                categorized_vcf[vartype]['0a0'] += count
            elif code[0] == '0' and code[1]=='0' and code[2] != '0': #code 001,002,003
                categorized_vcf[vartype]['00a'] += count
            elif code[0] != '0' and code[1]!="0" and code[2] == '0': 
                if(code[0] == code[1]): #code 110,220,330 
                    categorized_vcf[vartype]['aa0'] += count
                else: 
                    #code ab0 (codes were first and second called allele are different for two samples), divide into two seperate codes 
                    categorized_vcf[vartype]['a00'] +=count
                    categorized_vcf[vartype]['0a0'] +=count
                   # print("code:",code, "with" ,count, "counts added to pattern: a00 and 0a0")
            elif code[0] != '0' and code[1]=="0" and code[2] != '0': #code: 101,202,303
                if(code[0] == code[2]):
                    categorized_vcf[vartype]['a0a'] += count
                else: 
                    #a0b first and third alt different between two samples 
                    categorized_vcf[vartype]['a00'] +=count
                    categorized_vcf[vartype]['00a'] +=count
                    #print("code:",code, "with" ,count, "counts added to pattern: a00 and 00a")
            elif code[0] == '0' and code[1]!="0" and code[2] != '0': #011 022 033
                if(code[1] == code[2]):
                    categorized_vcf[vartype]['0aa'] += count
                else: 
                    #0ab codes were second and third alt allele are differnet between two samples 
                    categorized_vcf[vartype]['00a'] +=count
                    categorized_vcf[vartype]['0a0'] +=count
                    #print("code:",code, "with" ,count, "counts added to pattern: 00a and 0a0")
            elif code[0] != '0' and code[1]!="0" and code[2] != '0': 
                if(code[0] == code[1] == code[2] ):
                    categorized_vcf[vartype]['aaa'] += count #111, 222, 333
                elif(code[0] != code[1] and code[0] == code[2]): #"aba"
                    categorized_vcf[vartype]['a0a'] +=count
                    categorized_vcf[vartype]['0a0'] +=count
                elif(code[0] != code[1] and code[0] != code[2] and code[1]==code[2]): #"abb"
                    categorized_vcf[vartype]['a0a'] +=count
                    categorized_vcf[vartype]['0a0'] +=count
                elif(code[0] != code[1] and code[0] != code[2] and code[1] != code[2]): #"abc"
                    categorized_vcf[vartype]['a00'] +=count
                    categorized_vcf[vartype]['00a'] +=count
                    categorized_vcf[vartype]['0a0'] +=count
                    
                else:
                    print("code", code, "could not be added")
                    categorized_vcf[vartype]["other"] +=1
            else:
                #variants not recalled in vg 
                print(f"{count} of {vartype} not recalled in vg (000)")

        
    return categorized_vcf

def plot_venn_diagramm(single_categorized_vcf,key, outfile):
    
    #Plotting a 2x3 plot with all 6 diagramms (SNP, sINS, INS, sDEL, DEL and ALL)
    fig, axes = plt.subplots(nrows=2, ncols=(len(single_categorized_vcf) + 1) // 2, figsize=(15, 10))

    # Flatten the 2D array of subplots for easier indexing
    axes = axes.flatten()

    # Iterate through each vartype and corresponding subplot
    for i, (vartype, ax) in enumerate(zip(single_categorized_vcf.keys(), axes)):
        counts = single_categorized_vcf[vartype]
        ax.set_title(f"{key} overlap: {vartype}")
        venn3(subsets=(
            counts['a00'], counts['0a0'], counts['aa0'],
            counts['00a'], counts['a0a'], counts['0aa'],
            counts['aaa']),
            set_labels=('valid_2015T', 'valid_2020T', 'valid_2000B'),
            alpha=0.5,
            ax=ax,
            normalize_to=0.5)

    # Remove any remaining empty subplots
    for i in range(len(single_categorized_vcf), len(axes)):
        fig.delaxes(axes[i])

    # Adjust layout
    plt.tight_layout()

    #save_plots
    plt.savefig(outfile, bbox_inches='tight')
    # Show the plot
   # plt.show() 

       
per_file_dict=read_csv(csv_file_path)
per_file_categorized=get_catogorized_dict(per_file_dict)
for key in per_file_categorized.keys():
    plot_venn_diagramm(per_file_categorized[key], key, outfile[key])

            