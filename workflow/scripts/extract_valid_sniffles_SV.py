#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Jun  2 09:56:12 2023

@author: Marie
"""
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
    

def make_dir(infile, valid_sniffles_dic,allSV_dict,i): #Make directory from the sniffles data
    with open (infile) as sniffles:
        for line in sniffles: 
            if not line.startswith("#"):
               split_line=line.strip().split("\t")
               split_line[7]=split_line[7].strip().split(";")
               chrom=split_line[0]
               pos=split_line[1]
               key=(chrom, pos)
               if not key in allSV_dict:
                   allSV_dict[key]=[[],[],[]]
               allSV_dict[key][i].append(split_line[2:])
               if check_valid(split_line) == True:
                   if not key in valid_sniffles_dic:
                       valid_sniffles_dic[key]=[[],[],[]]  
                   valid_sniffles_dic[key][i].append(split_line[2:])                   
        return valid_sniffles_dic, allSV_dict
    

def extract_sv_type(sublist, sv_list,i):   #Get SV type
    if len(sublist) >0:
        sv=sublist[5][1].split("=")[1]
        sv_list[i].append(sv)
        return sv_list
 
def getSV_list(input_dict, outlist): #
    #sv_list=[]
    for key in input_dict.keys():
        for pos in range(0,3):
            if((input_dict[key])[pos]):
                for sub in range(len(input_dict[key][pos])):
                    sublist_i=((input_dict[key])[pos][sub])
                #print(sublist_i[7])
                #valid_SV.append()
                    outlist=extract_sv_type(sublist_i, outlist, pos)
    return outlist

###extract valid SV###
def extract_valid(sniffles_infile, valid_vcf_out, valid_entries):
    #print(sniffles_infile, valid_vcf_out, len(valid_entries))
    
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

sniffles_2000B="sniffles2-2_2000B_with_reference_Genbank.vcf"
sniffles_2020T="sniffles2-2_2020T_with_reference_Genbank.vcf"
sniffles_2015T="sniffles2-2_2015T_with_reference_Genbank.vcf"

file_list = [sniffles_2015T, sniffles_2020T, sniffles_2000B]

###outfiles####

valid_vcf_2015="sniffles2-2_2015T_with_reference_Genbank_valid.vcf"
valid_vcf_2020="sniffles2-2_2020T_with_reference_Genbank_valid.vcf"
valid_vcf_2000="sniffles2-2_2000B_with_reference_Genbank_valid.vcf"

outfile_list=[valid_vcf_2015,valid_vcf_2020,valid_vcf_2000]


####MAIN#####

sniffles_dic = {}
all_SV_dict={}  # Dictionary to store the combined lines
fcnt=0
for file in file_list:
    sniffles_dic, all_SV_dict=make_dir(file, sniffles_dic, all_SV_dict, fcnt)
    
    fcnt +=1 
 

valid_2015=[]
valid_2020=[]
valid_2000=[]
    
for e in sniffles_dic.keys():
    if sniffles_dic[e][0]:
        valid_2015.append(e)
    if sniffles_dic[e][1]:
        valid_2020.append(e)
    if sniffles_dic[e][2]:
        valid_2000.append(e)

valid_list=[valid_2015, valid_2020, valid_2000]

for i in range(0,len(file_list)):                     
    extract_valid(file_list[i], outfile_list[i], valid_list[i])


###COUNTING###
all_sv_list=[[],[],[]]
valid_SV_list=[[],[],[]]


all_SV_types=getSV_list(all_SV_dict, all_sv_list)
valid_SV_types=getSV_list(sniffles_dic, valid_SV_list)

sv_list=[all_SV_types, valid_SV_types]

for i in sv_list:
    count_dic={}  
    sublist_keys=["2015T", "2020T", "2000B"]
    for ele in range (len(i)):
        subl=i[ele]
        k=sublist_keys[ele]
        element_count={}
        for element in subl:
            if element in element_count:
                element_count[element] +=1
            else: 
                element_count[element] =1
        count_dic[k]=element_count

    for key, c_dict in count_dic.items():
        sum_all=0
        print(key +": ")
        for element , count in c_dict.items():
            print(element + ": "+ str(count))
            sum_all=sum_all + count
        print("total= ", sum_all)
        print()
