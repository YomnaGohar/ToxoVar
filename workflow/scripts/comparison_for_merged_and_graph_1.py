#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:54:49 2024

@author: yomna
"""
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import os
samples = 2
#int(sys.argv[1])
input_file = "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref.vcf"
#sys.argv[2]
#
output_dir=os.path.dirname(sys.argv[3])

# Initialize counts
counts = [[0 for _ in range(5)] for _ in range(samples)]
vartypes = ["SNV", "sIns", "sdel", "INS", "DEL"]
counts_vartype_disagree = {vt: [[0 for _ in range(7)] for _ in range(samples)] for vt in vartypes} #[count of (0,.),(.,0),(x,.),(.,x),(x,0),(0,x),(x1,x2)]
counts_vartype_agree = {vt: [[0 for _ in range(3)] for _ in range(samples)] for vt in vartypes} #[count of (.,.),(0,0),(x,x))
vartype_count_merged={vt: [0 for _ in range(samples)] for vt in vartypes}
vartype_count_vg={vt: [0 for _ in range(samples)] for vt in vartypes}
vartype_overlap_vg={vt: 0 for vt in vartypes}

def update_vartype_disagree(vtype, i, category):
    counts_vartype_disagree[vtype][i][category] += 1

def update_vartype_agree(vt, i, category):
        counts_vartype_agree[vtype][i][category] += 1
def update_counts(merged, vg, i, vtype):
    if merged == "." or merged == "-":
        if vg == ".":
            counts[i][4] += 1  # category 5 .,.
            update_vartype_agree(vtype, i, 0)
        else:
            counts[i][2] += 1  # category 3 .,1
            if vg != "0": #if the allele is not a reference
                update_vartype_disagree(vtype, i, 3) #(.,x)
                vartype_count_vg[vtype][i] +=1
                #print("merged == . and vg  != 0")
                #print(i,vg)
            else:
                update_vartype_disagree(vtype, i, 1) #((.,0))
    else:
        if merged == vg:
            counts[i][0] += 1  # category 1 1,1
            if merged == "0":
               update_vartype_agree(vtype, i, 1)  
            else:
                update_vartype_agree(vtype, i, 2)
                vartype_count_merged[vtype][i] +=1
                vartype_count_vg[vtype][i] +=1
        elif vg == ".":
            counts[i][3] += 1  # category 4 1,.
            if merged == "0":
                update_vartype_disagree(vtype, i, 0)  
            else:
                update_vartype_disagree(vtype, i, 2)
                vartype_count_merged[vtype][i] +=1
        else:
            counts[i][1] += 1  # category 2 1,0
            if merged == "0" and int(vg) > 0:
               update_vartype_disagree(vtype, i,5)
               vartype_count_vg[vtype][i] +=1
               #print("merged == 0 and int(vg) > 0")
               #print(i,vg)
            elif int(merged) > 0 and vg == "0":
                 update_vartype_disagree(vtype, i,4)
                 vartype_count_merged[vtype][i] +=1
            elif int(merged) > 0 and int(vg) > 0:
                update_vartype_disagree(vtype, i,6)
                vartype_count_merged[vtype][i] +=1
                vartype_count_vg[vtype][i] +=1
                #print("int(merged) > 0 and int(vg) > 0")
                #print(i,vg)

with open(input_file) as file:
   for r in file:
       if r.startswith("#"):
          continue
       genotype=[]
       for i in range(samples):
            z = r.strip().split('\t')
            merged = z[i + 9]
            vg = z[i + 9 + samples]
            genotype.append(vg)
            variant_type=[i.split("-")[-3] for i in z[2].split(",")]
            if merged != "." and  merged != "-" and merged != "0":
                variant_type_1 = variant_type[int(merged) - 1] 
            elif vg != "." and vg != "0":
                  variant_type_1 = variant_type[int(vg) - 1] 
            else:
                  variant_type_1 = variant_type[0]
            if "INS" == variant_type_1:
                if "," in z[4]: 
                    if merged != "." and  merged != "-" and merged != "0":
                        vt_mod = z[4].split(",")[int(merged) - 1] 
                    elif vg != "." and vg != "0":
                          vt_mod = z[4].split(",")[int(vg) - 1] 
                    else:
                          vt_mod = z[4].split(",")[0]
                else:
                    vt_mod=z[4]
                if len(vt_mod) >= 50:
                    vtype="INS"
                else: 
                    vtype="sIns"
            elif  "DEL" in variant_type:
                vt_mod = z[3] #if it's deletion then we look at the length of the reference
                if len(vt_mod) >= 50:
                    vtype="DEL"
                else:
                     vtype="sdel"
            elif z[4] == "SNV":
                vtype="SNV"
            update_counts( merged, vg, i, vtype)
       if len(set(genotype)) == 1 and genotype[0] != "0" and  genotype[0] != "."  and genotype[0] != "-":
          vartype_overlap_vg[vtype] += 1
# Convert counts to percentages
percentages = []
for sample_counts in counts:
    total = sum(sample_counts)
    if total == 0:
        percentages.append([0] * 5)  # Avoid division by zero
    else:
        percentages.append([(count / total) * 100 for count in sample_counts])
# Transpose percentages for plotting
percentages_transposed = list(zip(*percentages))
#Categories
categories = [
    "Graph and merged analyses both called an allele, and their calls agree (e.g., 1,1)",
    "Graph and merged analyses both called an allele, but their calls disagree (e.g., 0,1)",
    "Merged analysis did not call any, but the Graph called an allele (e.g., .,1)",
    "Merged analysis did, but Graph did not call an allele (e.g., 1,.)",
    "Neither graph nor merged analysis called an allele (.,.)"]
# Plotting
x = np.arange(samples)  # Sample indices
width = 0.15  # Width of the bars
fig, ax = plt.subplots(figsize=(12, 6))
for i, (category, percentage) in enumerate(zip(categories, percentages_transposed)):
    ax.bar(x + i * width, percentage, width, label=category)
# Add labels, title, and legend
ax.set_xlabel('Samples')
ax.set_ylabel('Percentage')
ax.set_title('Combined Bar Plot of Categories for Each Sample (Percentage)')
ax.set_xticks(x + width * 2)
ax.set_xticklabels([f'Sample_{i+1}' for i in range(samples)])
ax.legend()
output_plot = os.path.join(output_dir, 'combined_bar_plot.pdf')
plt.savefig(output_plot, bbox_inches='tight')
#plt.show()

##################################################################################################################################
#disagree
# Calculate percentages
categories = ["(0,.)", "(.,0)", "(x,.)", "(.,x)", "(x,0)", "(0,x)", "(x1,x2)"]
percentages = {vt: [] for vt in vartypes}
for vt in vartypes:
    for sample in range(samples):
        percentages[vt].append([(count / total) * 100 for count in counts_vartype_disagree[vt][sample]])

# Plotting
for i in range(samples):
    fig, ax = plt.subplots(figsize=(14, 8))
    indices = np.arange(len(vartypes))  # Variant types indices
    width = 0.1  # Width of the bars

    for j, category in enumerate(categories):
        category_percentages = [percentages[vt][i][j] for vt in vartypes]
        ax.bar(indices + j * width, category_percentages, width, label=category)

    # Add labels, title, and legend
    ax.set_xlabel('Variant Types')
    ax.set_ylabel('Percentage')
    ax.set_title(f'Variant Types for Sample {i + 1}')
    ax.set_xticks(indices + width * (len(categories) - 1) / 2)
    ax.set_xticklabels(vartypes)
    ax.legend()
    output_plot = os.path.join(output_dir, f'sample_{i + 1}_bar_plot_disagree.pdf')
    plt.savefig(output_plot, bbox_inches='tight')
#    plt.show()
##################################################################################################################################
#agree
# Calculate percentages
categories = ["(.,.)","(0,0)","(x,x)"]
percentages = {vt: [] for vt in vartypes}
for vt in vartypes:
    for sample in range(samples):
        percentages[vt].append([(count / total) * 100 for count in counts_vartype_agree[vt][sample]])

# Plotting
for i in range(samples):
    fig, ax = plt.subplots(figsize=(14, 8))
    indices = np.arange(len(vartypes))  # Variant types indices
    width = 0.1  # Width of the bars

    for j, category in enumerate(categories):
        category_percentages = [percentages[vt][i][j] for vt in vartypes]
        ax.bar(indices + j * width, category_percentages, width, label=category)

    # Add labels, title, and legend
    ax.set_xlabel('Variant Types')
    ax.set_ylabel('Percentage')
    ax.set_title(f'Variant Types for Sample {i + 1}')
    ax.set_xticks(indices + width * (len(categories) - 1) / 2)
    ax.set_xticklabels(vartypes)
    ax.legend()
    output_plot = os.path.join(output_dir, f'sample_{i + 1}_bar_plot_agree.pdf')
    plt.savefig(output_plot, bbox_inches='tight')
#    plt.show()    
##########################################################################################################################################
#calculate overlap between samples in their variant types
