#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 07:57:56 2024

@author: marie
"""
import pandas as pd
import sys
#script to add reference and unplaced filtering to the merged vcf column of the combined resulttable vcf using medaka and sniffles raw data
combined_table=sys.argv[1]
num_samples=int(sys.argv[2])
samples=sys.argv[3:3+num_samples]
print(samples)
medaka=sys.argv[3+num_samples:3+2*num_samples]
print(medaka)
#medaka_2020T=sys.argv[3]
#medaka_2000B=sys.argv[4]

#Note: Sniffles 2 raw variants are not corrected for Vartype and alt allele:
sniffles=sys.argv[3+2*num_samples:3+3*num_samples]
print(sniffles)
#sniffles_2020T=sys.argv[6]
#sniffles_2000B=sys.argv[7]
agp_file=sys.argv[-2]
outfile=sys.argv[-1]

#create dictionary of samples where name as keys and file as value
m={}
s={}
for i in  range(num_samples):
    m[samples[i]]=medaka[i]
    s[samples[i]]=sniffles[i]

df = pd.read_csv(combined_table, sep='\t', dtype="str")

def read_agp(agp_infile):
    exclude_contigs=[]
    with open (agp_infile, "r") as agp:
        for line in agp:
            splitline=line.split("\t")
            if splitline[0] == ("TgMe49_00") or splitline[0] == ("TgMe49_API"):
                if(splitline[4]=="W"):
                    exclude_contigs.append(splitline[5])
    return(exclude_contigs)


def read_vcf_to_dataframe(vcf_file):
    # Read the VCF file into a DataFrame
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

    # Set column names
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

    return vcf_df


def check_variant_presence_medaka(vcf_df, chrom, pos, alt, variant_present):
    #print(chrom, pos, alt)
    variant_present=((vcf_df["CHROM"]==chrom) & (vcf_df["POS"]==int(pos))).any() #&(vcf_df['ALT']==alt)).any()
    
    #if variant_present == True:
        #print("True", chrom, pos, alt)
    #print((vcf_df['CHROM'] == chrom) & (vcf_df['POS'] == pos)) # & (vcf_df['ALT'].str.contains(alt))).any())

    return variant_present

def check_variant_presence_sniffles(vcf_df, chrom, pos, alt, variant_present):
    variant_present=((vcf_df["CHROM"]==chrom) & (vcf_df["POS"]==int(pos)+1)).any()
    return variant_present




unplaced_seq=read_agp(agp_file)
indices_to_remove = []  
for index, row in df.iterrows():
    vartype = row['Vartype']
    alt = row['ALT']
    ref = row['Ref']
    chrom = row['Chromosom']
    pos = row['Position']
    
    if vartype == 'INS':
        length = len(alt)-len(ref)
        if length < 50:
            vartype= "sINS"
    elif vartype == 'DEL':
        length = len(ref)-len(alt)
        if length < 50:
            vartype="sDEL"
    else: 
        length=1
    for sample in samples:
        variant_present=False
        merged_col = sample + '_merged'
        if(chrom in unplaced_seq):
            #row[merged_col] = "."
            indices_to_remove.append(index)  # Append index to list

        #print(merged_col)
        elif vartype in ["sINS", "sDEL", "SNV"]:
            #medaka_df_name=f'medaka_{sample}_df' 
            medaka_df = read_vcf_to_dataframe(m[sample])
            
            if row[merged_col] == "0":
                variant_present = check_variant_presence_medaka(medaka_df, chrom, pos, alt, variant_present)
                
        else:
            #sniffles_df_name=f'sniffles_{sample}_df' 
            sniffles_df =  read_vcf_to_dataframe(s[sample])
            
            if row[merged_col] == "0":
                variant_present = check_variant_presence_medaka(sniffles_df, chrom, pos, alt, variant_present)
        
        if variant_present:
            df.at[index, merged_col] = '.'
df.drop(index=indices_to_remove, inplace=True)            
#print(df)
df.to_csv(outfile, index=False, sep='\t')        
        