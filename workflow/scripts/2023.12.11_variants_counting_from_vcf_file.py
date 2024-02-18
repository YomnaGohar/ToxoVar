#!/usr/bin/env python

"""
Created on Wed Dec 13 20:40:39 2023
@author: yomna
"""

import pandas as pd
import sys

# Files:     
medaka_vcf = sys.argv[1] 
medaka_files = [medaka_vcf]
# Counting statistics from Vcf files

for file in range(0,len(medaka_files)):
    with open(medaka_files[file]) as medaka:
        medaka_dict={}
        variant_count_total = {"SNP":0,"INS":0,"DEL":0,"MNP":0,"MIXED":0}
        dupl_chrm_pos_nonIdentical = {"SNP":0,"INS":0,"DEL":0,"MNP":0,"MIXED":0}
        dupl_chrm_pos_Identical = {"SNP":0,"INS":0,"DEL":0,"MNP":0,"MIXED":0}
        unique_variants = {"SNP":0,"INS":0,"DEL":0,"MNP":0,"MIXED":0}
        for line in medaka: # to read the file line by line
             if not (line.startswith("#")):
                arr = line.strip().split("\t") # line.strip() is to remove white spaces from the line. split fuction will create the list  
                chom = arr[0]
                pos = arr[1]
                vartype=arr[7].split(";")[4].split("=")[1]
                if not vartype in variant_count_total: variant_count_total[vartype]=0 #important because some weird variant names arise and they are not in the dictionary
                if not(vartype in dupl_chrm_pos_nonIdentical): dupl_chrm_pos_nonIdentical[vartype]=0  #important because some weird variant names arise and they are not in the dictionary
                if not(vartype in dupl_chrm_pos_Identical): dupl_chrm_pos_Identical[vartype]=0  #important because some weird variant names arise and they are not in the dictionary
                if not(vartype in unique_variants): unique_variants[vartype]=0  #important because some weird variant names arise and they are not in the dictionary
                variant_count_total[vartype] +=1
                if not (chom,pos) in medaka_dict: # it's important to include this if statement because medaka_comp has repeated lines
                    medaka_dict[(chom,pos)]=[arr[2:]]  
                    unique_variants[vartype] +=1
                elif medaka_dict[(chom,pos)][0][2] == arr[4]:
                     dupl_chrm_pos_Identical[vartype] +=1
                     #print(medaka_dict[(chom,pos)])
                     #print(arr)
                     #print("\n")
                elif float(medaka_dict[(chom,pos)][0][3]) < float(arr[5]): 
                     vartype_1=medaka_dict[(chom,pos)][0][5].split(";")[4].split("=")[1]
                     dupl_chrm_pos_nonIdentical[vartype_1] +=1
                     medaka_dict[(chom,pos)]=[arr[2:]]   
                     unique_variants[vartype_1] -=1
                     unique_variants[vartype] +=1
                else:
                     dupl_chrm_pos_nonIdentical[vartype] +=1
 
 

    data = [variant_count_total, dupl_chrm_pos_nonIdentical, dupl_chrm_pos_Identical, unique_variants]
    df = pd.DataFrame(data, index=["Variant Count Total", "Dupl Chrom Pos NonIdentical", "Dupl Chrom Pos Identical", "Unique Variants"])
    print ("Results for sample:", medaka_vcf )
    print(df)
