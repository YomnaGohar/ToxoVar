#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 15:48:03 2023

@author: marie
script for filtering the valid snps and indels from the medaka calls, using mpileup information (reads are filtered for mapping and Base quality)
insert Path to medaka_Qual1 files, mpileup and outfiles 
"""

##automatic Filtering --> valid SNPs and Indels
# define functions
import re
import random
import pandas as pd
import sys
def add_mpileup_to_comb_medaka(mpileup,a,medaka_dict):
        #This function takes as an argument the mpileup file, sample name for this mpileup file and combined varients vcf 
        #as a form of dictionary (medaka_dict). for each line in the mpileup file, the function will check if it's present in the
        #combined vcf dictionary (medaka_dict), if yes it will update the dictionary with the information from 4th column in mpileup file. 
        #Input example, {('KE138874','964'): [['.','G','T','17.008', 'PASS', 'DPS=3,2;SNP;HOM;DP=5;VARTYPE=SNP',
        #  'GT:GQ', '1:17', '.:.']]}
        #Ouput example,{('KE138874','964'): [['.','G','T','17.008', 'PASS', 'DPS=3,2;SNP;HOM;DP=5;VARTYPE=SNP',
        #  'GT:GQ', '1:17', '.:.'],{'2015': ')=5;(=3\n'}]}
                                             
    with open(mpileup) as mpileup_2015:           
        for line_2 in mpileup_2015: # to read the file line by line
            arr_mp= line_2.strip().split("\t") # line.strip() is to remove white spaces from the line. split fuction will create the list  
            chom_mp = arr_mp[0]
            pos_mp = arr_mp[1]
            if (chom_mp,pos_mp) in medaka_dict: 
               if len(arr_mp) == 4:
                   depth_mp=arr_mp[3]
                   Toxo=depth_mp
                   Toxo_list=Toxo.split(";")  #creating a list from pileup data. ex, [A=5, C=6]
                   Toxo_dict=dict(subString.split("=") for subString in Toxo_list)
                   Toxo_dict_i=dict([c, int(x)] for c, x in Toxo_dict.items())
                   medaka_dict[(chom_mp,pos_mp)][a][0] =[Toxo,Toxo_dict_i]
                   sum_i=sum(Toxo_dict_i.values())
                   medaka_dict[(chom_mp,pos_mp)][a][1]=sum_i
   
def get_alt_idx(s_signal):
    #''' this method gets the index for the alternative base, to see from which sample, 
    #has which alt base'''
    ###only one index possible, only one sample --> might not be needed anymore ###
    idx="."
    signal_split=s_signal.split(":")[0]
    
    if(signal_split=="1"):
        idx=0
     
    return idx   
 
def check_filtering(sum_i, base_sum_c, filter_sample):
    if sum_i < 7 and base_sum_c >= 80.0:
        filter_sample[0] += 1
    elif sum_i >= 7 and base_sum_c < 80.0:
        filter_sample[1] += 1
    elif sum_i < 7 and base_sum_c < 80.0:
        filter_sample[2] += 1
    elif sum_i >= 7 and base_sum_c >= 80.0:
        filter_sample[4] += 1   
    return (filter_sample)

def add_eval (signal,a,b,Toxo, vartype, filtering_sample):
    #This function adds an evaluation score into the medaka dictionary, when the alternative 
    #base called is valid (1) or not (0), valid means that 80% of the reads are agreening and more than 7 are counted
    ### adjusted for insertions and Deletions, because those are not only 1 bp in length, compared to SNPs. 
    ###SNPs --> keys[0].lower() only used the first base to compare to the alt.base
    ###INS ---> using the complete key to compare the alt base to the mpileup informations (option to also look into smaller or larger insertions reported in mpileup)
    ###DEL --> Medaka reports the base before the deletion, so Altbase needs to be g-, g--- ... not g
    
    if signal != '.:.':
      #   print("ALTBASE: ", Alt_base)
         if Alt_base[b].lower() in Toxo:
            sum_i=sum(medaka_dict[i][a][0][1].values())
            #print("sum_i: ", sum_i)
            base_sum=[]
            Toxo_dict_i=medaka_dict[i][a][0][1]
            #print("Toxo_dict_i:", Toxo_dict_i)
            for keys in Toxo_dict_i:
             #   print("keys", keys)
                if vartype =="SNP":
              #      print("SNP: keys, altbase", keys[0].lower(), ", ", Alt_base[b].lower() )
                    if keys[0].lower() == Alt_base[b].lower():
                        base_sum.append(Toxo_dict_i[keys])
                        
                if vartype =="INS":
               #     print("INS: keys, altbase", keys.lower(), ", ", Alt_base[b].lower() )
                    if keys.lower() == Alt_base[b].lower():
                       base_sum.append(Toxo_dict_i[keys])
                       
                if vartype =="DEL":
                #    print("DEL:")
                    alt_base_pattern = re.escape(Alt_base[b].lower()) + '-+'
                    alt_base_regex = re.compile(alt_base_pattern, re.IGNORECASE)
                    if alt_base_regex.match(keys.lower()):
                 #       print("DEL: keys, altbase","match:" , ", ", keys.lower(), ", ", Alt_base[b].lower())
                        base_sum.append(Toxo_dict_i[keys])
                    #else: 
                        #print("no match")
                #print("base_sum", base_sum)
            base_sum_c=sum(base_sum)/sum_i*100
                
           # print("base_sum_c:", base_sum_c)
           # print("sum_i", sum_i, "\n")
                
            filtering_sample[vartype]=check_filtering(sum_i, base_sum_c, filtering_sample[vartype])
                
            if sum_i >= 7 and base_sum_c >= 80:
                medaka_dict[i][a][2]=1
            else: 
                medaka_dict[i][a][2]=0
         else:
             filtering_sample[vartype][3] +=1
             medaka_dict[i][a][2]=0 
            
         '''
         ##Version of including insertions of different length from the mpileup file not comparing only CTT for example 
         but also CT or CTTT as insertions at this positions##
         
         if vartype =="INS":
            print(keys.lower(), ", ", Alt_base[b].lower() )
            alt_base_pattern = re.escape(Alt_base[b].lower())
            if(len(Alt_base[b].lower()) > 2):
                alt_base_pattern = '(?:' + alt_base_pattern + '|' + re.escape(Alt_base[b].lower()[:2]) + ')'
            else:
                alt_base_pattern = re.escape(Alt_base[b].lower()) + '.*'
                alt_base_regex = re.compile(alt_base_pattern, re.IGNORECASE)
                if alt_base_regex.match(keys.lower()):
                    print("true:", keys.lower(), "add: ", Toxo_dict_i[keys])
                    base_sum.append(Toxo_dict_i[keys])
                else: 
                    print("no match: ", keys.lower())
         '''
                   
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
                            

#Files: 
medaka_vcf=sys.argv[1]

medaka_files=[medaka_vcf]    #Can also be run for one sample at the time by inserting only one file to the list

mpileup = sys.argv[2]


mpileup_files=[mpileup]

### outfiles: valid_vcf 
valid_vcf= sys.argv[3]
#statistics=snakemake.output["stats"]

valid_outfiles=[valid_vcf] #Can also be run for one sample at the time by inserting only one file to the list

#The with statement itself ensures proper acquisition and release of resources. for more information check https://www.geeksforgeeks.org/with-statement-in-python/
for file in range(0,len(medaka_files)):
    filtering_sample={"SNP":[0,0,0,0,0],"INS":[0,0,0,0,0],"DEL":[0,0,0,0,0],"MNP":[0,0,0,0,0],"MIXED":[0,0,0,0,0]}
    unique_variants={"SNP":0,"INS":0,"DEL":0,"MNP":0,"MIXED":0}
    with open(medaka_files[file]) as medaka:
        medaka_dict={}
        for line in medaka: # to read the file line by line
            if not (line.startswith("#")):
                arr = line.strip().split("\t") # line.strip() is to remove white spaces from the line. split fuction will create the list  
                chom = arr[0]
                pos = arr[1]
                vartype=arr[7].split(";")[4].split("=")[1]
                if not(vartype in unique_variants): unique_variants[vartype]=0
                if not(vartype in filtering_sample): filtering_sample[vartype]=[0,0,0,0,0]
                if not (chom,pos) in medaka_dict: # it's important to include this if statement because medaka_comp has repeated lines
                    medaka_dict[(chom,pos)]=[arr[2:],[['',{}],0,'','']]    
                    unique_variants[vartype] +=1
                elif float(medaka_dict[(chom,pos)][0][3]) < float(arr[5]):  
                     vartype_1=medaka_dict[(chom,pos)][0][5].split(";")[4].split("=")[1]
                     unique_variants[vartype_1] -=1
                     unique_variants[vartype] +=1
                     medaka_dict[(chom,pos)]=[arr[2:],[['',{}],0,'','']]
                     

                     

    # add mpileup information
    add_mpileup_to_comb_medaka (mpileup_files[file],1,medaka_dict)

    #evaluation: 
    cnt=0
    for i in medaka_dict: 
       # print("entry:", i)
        Alt_base=medaka_dict[i][0][2].split(',')
            
        ###get the signal for each sample 
        signal_i=medaka_dict[i][0][7] 
       
        ###get the mpileup part of the medaka dictionary for each sample 
        Toxo_i=medaka_dict[i][1][0][0].lower()
    
        ###get the index for the alternative base for each sample
        idx_i= get_alt_idx(signal_i)
    
        ####perform evaluation
        var_type=(medaka_dict[i][0][5].split(";")[2])
        add_eval(signal_i,1,idx_i,Toxo_i, var_type, filtering_sample)
        cnt=cnt+1
        
    #Get Valid SNPs and indels:
    valid_variations = []    
    not_valid_variations=[]          
    for e in medaka_dict:  
        if medaka_dict[e][1][2] == 1:
            c=e[0]
            p=e[1]
            alt=medaka_dict[e][0][2]
            qual=medaka_dict[e][0][3]
            valid_variations.append((c,p,alt,qual))
        
        else: 
            not_valid_variations.append(e)
    MapQ_30_baseQ_20={"SNP":unique_variants["SNP"]-filtering_sample["SNP"][3],
                      "INS":unique_variants["INS"]-filtering_sample["INS"][3],
                      "DEL":unique_variants["DEL"]-filtering_sample["DEL"][3]}
    Depth_7= {"SNP":MapQ_30_baseQ_20["SNP"]-filtering_sample["SNP"][0]-filtering_sample["SNP"][2],
              "INS":MapQ_30_baseQ_20["INS"]-filtering_sample["INS"][0]-filtering_sample["INS"][2],
              "DEL":MapQ_30_baseQ_20["DEL"]-filtering_sample["DEL"][0]-filtering_sample["DEL"][2]}  
    valid={"SNP":Depth_7["SNP"]-filtering_sample["SNP"][1],
              "INS":Depth_7["INS"]-filtering_sample["INS"][1],
              "DEL":Depth_7["DEL"]-filtering_sample["DEL"][1]}     
    #print("filtered_by_depth", "filtered_by_allelefreq", "filtered_by_both", "filtered_by_signal", "valid")
    # filtered_by_depth: excluded because of depth smaller 7 although allele frequency >= 80%, 
    # filtered_by_allelefreq: allele frequency smaller 80% although depth >= 7, 
    # filtered_by_both: variant filtered because of depth and allelefreq, 
    # filtered_by_signal:  .:. signal instead of a genotype --> medaka did not assign a genotype
    # valid: valid calls 
    #print(filtering_sample)
    # Creating a list of dictionaries
    data = [MapQ_30_baseQ_20, Depth_7, valid]

    # Creating the DataFrame
    df = pd.DataFrame(data, index=["at least one base MapQ>30&baseQ>20 ",
                                   "Depth >= 7", "Allele Freq >= 80"])
    # Display the DataFrame
    print(df)
                         

        
#writing outfiles for all three samples including valid calls:
    extract_valid(medaka_files[file], valid_outfiles[file], valid_variations)    


