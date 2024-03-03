#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 09:08:28 2023

@author: marie
"""
import pysam
#Homopolymercheck Medaka calls: 
# Are the SNPs and Indel calls from Medaka part of homopolymer structures: 
# The script tests if the base before a SNP or a insertion is identical to the base of the event reported in the VCF 

#IMPORTANT: Medaka is one base before the reference sequence when read in using pysam (pos-1 is the base medaka calls as reference base! Thats why base at pos-2 is the base of interest in the reference sequence)
    
medaka_infile="/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2015T_graph_Alignment/2015T_variants_MQ30_BQ20_vartype_filter_mod.vcf"
reference_seq="/home/yomna/Desktop/PhD_Yomna_Gohar/fasta/ME49_genome_from_Third_generation_sequencing_paper/ToxoME49_GenBank_JACEHA000000000.1.fsa"
n_bases=1 #number of bases before the event that should be identical (default 1 base)
  
# def homopolymer_check(vcf_data,reference, nbases, vartype):
#     homopolymer=0
#     no_homopolymer=0
    
#     for i in vcf_data.keys():
#         if (vcf_data[i][4]["VARTYPE"]==vartype):
#             pos=i[1]
#             chrom=i[0]
#             ref_seq_before=(reference_genome[chrom][int(pos)-(nbases+1):int(pos)]) 
#             #nbases + 1 is need because the reference sequence is 0 based and medaka vcf 1 based)
#             if all(char == ref_seq_before[1] for char in ref_seq_before):
#                 homopolymer +=1
#             else:
#                 no_homopolymer +=1
#         #else: 
#         #    print(vcf_data[i][4]["VARTYPE"])
            
#     total=homopolymer+no_homopolymer
#     print("total", vartype ,"calls=", total, "\n" ,"checking if " ,nbases, "bases before are identical to", 
#           vartype,"\n", f"homopolymeric {vartype}:", homopolymer, " (", round(homopolymer*100/total, 2), "% )" "\n" ,
#           f"no_homopolymeric {vartype}:", no_homopolymer, "(", round(no_homopolymer*100/total, 2),"% )", "\n")

def homopolymer_check_total(vcf_data,reference, vartype):
    homopolymer=0
    no_homopolymer=0
    print(vartype)
    for i in vcf_data.keys():
        if (vcf_data[i][4]["VARTYPE"]==vartype):
            pos=int(i[1])
            #print(pos)#medaka position 
            chrom=i[0]
            ins=(vcf_data[i][1]) #ins sequence 
            #print(ins)
            #len_ins=len(vcf_data[i][1]) #length of the insertion 
            ref_pos=pos-1 #insertion is added at this position
            ref_seq_before=(reference_genome[chrom][int(ref_pos)])
            #ref_seq_before=(reference_genome[chrom][ref_pos-len_ins:int(ref_pos)]) #get sequence from before the ins pos. check if it is identical to the inserted sequence, start:ref_pos-1 (exclusive)
            #print(ref_seq_before, "\n")
            
            #nbases + 1 is need because the reference sequence is 0 based and medaka vcf 1 based)
            if ins == ref_seq_before:
                homopolymer +=1
                print(ins)
                print(ref_seq_before)
                
            else:
                no_homopolymer +=1
                print(ins)
        #else: 
        #    print(vcf_data[i][4]["VARTYPE"])
            
    total=homopolymer+no_homopolymer
    print("total", vartype ,"calls=", total, "\n" ,"checking if inserted sequence is present in the reference before the event:", 
          vartype,"\n", f"homopolymeric {vartype}:", homopolymer, " (", round(homopolymer*100/total, 2), "% )" "\n" ,
          f"no_homopolymeric {vartype}:", no_homopolymer, "(", round(no_homopolymer*100/total, 2),"% )", "\n")
    
    
reference_genome = {}

# Parse reference genome and store it in the dictionary
with pysam.FastaFile(reference_seq) as fasta:
    for contig in fasta.references:
        reference_genome[contig] = fasta.fetch(contig)

vcf_data = {}

# Read VCF 
with open(medaka_infile, 'r') as file:
    for line in file:
         # Skip header lines
         if line.startswith('#'):
              continue

         # Split line into fields, store entry in (list)
         fields = line.strip().split('\t')
         chrom = fields[0]  # Chromosome
         pos = int(fields[1])  # Position
         ref = fields[3]  # Reference allele
         alt = fields[4]  # Alternate alleles
         qscore = fields[5]  # QScore
         filtercol = fields[6]  # filtercolumn

         # INFO field stored in dictionary
         info = {}
         info_fields = fields[7].split(';')
         #print(info_fields)
         info["DP"] = info_fields[0].split("=")[1]
         if info_fields[1].split("=")[0] == "DPS": #medaka vcf files
             info["DPS"] = info_fields[1].split("=")[1]
             info["VARTYPE"] = info_fields[4].split("=")[1]
         else:                                               # vg vcf files               
             #info["VARTYPE"] = info_fields[3].split("=")[1]    
             info["VARTYPE"] = info_fields[1]
         gtgq = fields[9]

         #each entry is combined into a list
         variant = [ref, alt, qscore, filtercol, info, gtgq]
         vcf_data[chrom, pos] = variant

#checking for homopolymeric sequence before SNPs and Ins:
#homopolymer_check(vcf_data, reference_genome, n_bases, "SNP")
#homopolymer_check(vcf_data, reference_genome, n_bases, "INS")

homopolymer_check_total(vcf_data, reference_genome, "INS")

