#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:34:25 2023

@author: marie
"""

import pysam
import sys

# Define paths to your VCF file and reference genome FASTA file
vcf_file =sys.argv[1]
# "/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/snakemake_Yomna/analysis/Sniffles/2015T/sniffles_2015T_with_reference_valid.vcf"
reference_genome_file = sys.argv[2]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/fasta/ME49_genome_from_Third_generation_sequencing_paper/ToxoME49_GenBank_JACEHA000000000.1.fsa"
outfile=sys.argv[3]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/snakemake_Yomna/analysis/Sniffles/2015T/sniffles_2015T_with_reference_valid_corrected.vcf"
# Dictionary to store reference genome sequences
reference_genome = {}

def get_mpileup(chromosome, position, mpileup_infile):
    #method for parsing the mpileup and returning the mpileup information of a specific Chromosom Position combination 
    #just for checking if changing is working correctly
    
    splitline=[]
    position=str(position)
    #print(position)
    with open (mpileup_infile, "r") as mpileup:
        for line2 in mpileup:
             splitline=line2.split("\t")
             if(splitline[0]==chromosome and splitline[1]==position):
                #print(splitline)
                return splitline
    
# Parse reference genome and store it in the dictionary
with pysam.FastaFile(reference_genome_file) as fasta:
    for contig in fasta.references:
        reference_genome[contig] = fasta.fetch(contig)

# Open Sniffles VCF file for reading and create an output VCF file for writing the corrected lines
with open(vcf_file, 'r') as vcf_in, open(outfile, 'w') as vcf_out:
    for line in vcf_in:
        if line.startswith('#'):
            #print(line)
            vcf_out.write(line) #Write Header Lines
        else:
            ##Lines are splitted by Tab, extracted are chromosom, position, reference and alternative Base and the info column
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, info, _ = fields[:9]
            var_type=(info.split(";")[1].split("=")[1])
            ## mpileup can be used to see where the SV is reported in mpileup
            #mpileup_info=[] #mpileup_info=get_mpileup(chrom,pos, mpileup) #print("mpileup_reference base is:", mpileup_info [2])
            
            ##Correcting the different Types of SV in the Sniffles VCF:
            #Deletions have an N in ALT BASE: The new Altbase should be the Base at position one base earlier in the reference sequence;
            #Position-1 is the new position of the SV, The base at Pos-1 is added to the reference and the Base at the Position is the new ALT allele
            
            if var_type=="DEL":
                
                pos=int(pos)-1
                pos=str(pos)
    
                ref=reference_genome[chrom][int(pos)-1]+ref
                alt=reference_genome[chrom][int(pos)-1]
                
            #Insertion have N in the Reference Base, in VCF the Base before the INsertion should be the Reference Base
            # The Position is corrected to Pos-1 base. The referencebase is the base at this position in the reference 
            # To the Altbase the base is also added
            
            if var_type=="INS":
                
                pos=int(pos)-1
                pos=str(pos)
    
                ref=reference_genome[chrom][int(pos)-1]
                alt=reference_genome[chrom][int(pos)-1]+alt
                
            #INVERSIONS:
            if var_type=="INV":
                ref=reference_genome[chrom][int(pos)]
                     
            if var_type=="BND":
                print("There is a BND in the valid file")
                
            if var_type=="DUP":
                print("There is a duplication in the valdi file")
                
            # Reconstruct the modified line and write it to the output VCF file
            fields[1] = pos
            fields[3] = ref
            fields[4] = alt
            
            vcf_out.write('\t'.join(fields) + '\n')
