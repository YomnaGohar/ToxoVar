#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 08:20:27 2023

@author: marie

Script for filtering and counting the variant calls made by vg call --> regenotyped variants from the merged vcf using the merged vcf in vg call
"""

import pandas as pd
import sys

#Read in the VCF file per sample, VCF includes all call, variants (1) and reference calls (0) 
def read_vcf(vcf):
    vcf_dict={}
    with open (vcf, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                if not (chrom,pos) in vcf_dict: 
                    vcf_dict[(chrom,pos)]=[splitline[2:], []]   
                else: 
                    print(chrom, pos, "already in dictionary")
    return vcf_dict
def read_agp(agp_infile,excl_contig):
    print(excl_contig == "None")
    if excl_contig == "None": 
        exclude_contigs=[]
        with open (agp_infile, "r") as agp:
            for line in agp:
                splitline=line.split("\t")
                if splitline[0] == ("TgMe49_00") or splitline[0] == ("TgMe49_API"):
                    if(splitline[4]=="W"):
                        exclude_contigs.append(splitline[5])
    else:
        exclude_contigs = excl_contig    
    #print(exclude_contigs)             
    return(exclude_contigs)

# Method for filtering and validating calls
# Filtering criteria are: "1 Genotype (variant)", "PASS", Readdepth > 7, Allelefreq > 80% 
# Filtering results per step are printed in a DF (get_filter_table.py)
def ref_criteria_check(variant, vcf_dict,key, unplaced_seq):
    map_info=variant[7].split(":")
    if ("," in variant[2]):
        #print(variant)
        #print(map_info[2])
        allele_depth=(map_info[2].split(","))
        allele_depth=[int(s) for s in allele_depth]
        max_alt=(max(allele_depth[1:]))
        idx=allele_depth.index(max_alt)
        if (idx == 0): #how can idx be zero if max_alt is search in the alternative alleles
            variant[2]=variant[2].split(",")[0]
        else:
            variant[2]=variant[2].split(",")[idx-1]
       # print(variant[2])
    if variant[4] != "PASS":
        vcf_dict[key][1].append(f"notPass, {variant[4]}")
        vcf_dict[key][1][0]="."
    else:
        if int(map_info[1]) < depth:
            vcf_dict[key][1].append(f"lowReadDepth, {map_info[1]}")
            vcf_dict[key][1][0]="."
        else:
            DP=int(map_info[1])
            REF_DP=int(map_info[2].split(",")[0])
            ALT_freq=(REF_DP*100)/DP
            if ALT_freq < alfreq:
                r_alt_freq=round(ALT_freq,2)
                vcf_dict[key][1].append(f"lowAllelefreq, ({r_alt_freq})")
                vcf_dict[key][1][0]="."
            elif(key[0] in unplaced_seq):
                vcf_dict[key][1].append(f"unplaced, ({key[0]})" )
                vcf_dict[key][1][0]="."
                #cnt=cnt+1
            else:
                vcf_dict[key][1].append(["PASS", DP, round(ALT_freq,2)])
             
    return vcf_dict           

def variant_criteria_check(variant, vcf_dict,key, which_altallele, unplaced_seq):
    if("," in variant[2]):
        variant[2]=variant[2].split(",")[int(which_altallele)-1]
    if variant[4] != "PASS":
        vcf_dict[key][1].append(f"notPass, {variant[4]}")
    else:
        map_info=variant[7].split(":")
        if int(map_info[1]) < depth:
            vcf_dict[key][1].append(f"lowReadDepth, {map_info[1]}")
        else:
            DP=int(map_info[1])
            ALT_DP=int(map_info[2].split(",")[int(which_altallele)])
            ALT_freq=(ALT_DP*100)/DP
            if ALT_freq < alfreq:
                r_alt_freq=round(ALT_freq,2)
                vcf_dict[key][1].append(f"lowAllelefreq, ({r_alt_freq})" )
            elif(key[0] in unplaced_seq):
                vcf_dict[key][1].append(f"unplaced, ({key[0]})" )
            else:
                vcf_dict[key][1].append(["PASS", DP, round(ALT_freq,2)])
                vcf_dict[key][1][0]=1
    return vcf_dict
def validate_calls(vcf_dict, stats_outfile, unplaced_seq):
    for key in vcf_dict.keys():
        variant=vcf_dict[key][0]
        vcf_dict[key][1].append(0)
        which_altallele=variant[7].split(":")[0]
        #exclude calls thst did not PASS all filtering criteria
        if variant[7].split(":")[0] == "0":
            vcf_dict[key][1].append("ref_call")
            vcf_dict=ref_criteria_check(variant,vcf_dict,key, unplaced_seq)
            #if reference vcf_dict[key]: [[varianat],[0,"ref_call",["Pass",DP,freq]]]
        else:
            vcf_dict=variant_criteria_check(variant,vcf_dict,key,which_altallele, unplaced_seq)   
            #if variant:  vcf_dict[key]: [[varianat],[0,["Pass",DP,freq]]]
        VarType=categorise_call(variant)
        vcf_dict[key][1].append(VarType)            
    vcf_dict=get_filter_table(vcf_dict, stats_outfile)
    return (vcf_dict)
#Method for extracting a filtering table, showing in each row the different filtering criteria and how many variants are remaining after each step per varinat type
def get_filter_table(vcf_dict, stats_outfile):
    filtering_dict={"SNP": [0,0,0,0,0,0,0], "sINS": [0,0,0,0,0,0,0], "INS": [0,0,0,0,0,0,0],"sDEL": [0,0,0,0,0,0,0],"DEL": [0,0,0,0,0,0,0], "MIXED":[0,0,0,0,0,0,0]}
    #per Vartype: [total, ref_call, not_pass (lowad, lowdepth), readdepth < depth, allelefreq,  valid]
    ref_calls={"SNP":[0,0,0,0,0],"sINS": [0,0,0,0,0], "INS": [0,0,0,0,0],"sDEL": [0,0,0,0,0],"DEL": [0,0,0,0,0], "MIXED":[0,0,0,0,0]} #[refcall,not_pass (lowad, lowdepth), readdepth < depth, allelefreq]
    for variant in vcf_dict.keys():
        #total (valid+not_valid calls)
        filtering_dict[vcf_dict[variant][1][-1]][0] +=1
        if(vcf_dict[variant][1][0]==1): #yes if it is not a reference call and valid
            #valid calls 
            filtering_dict[vcf_dict[variant][1][2]][6] +=1
        elif (vcf_dict[variant][1][0]=='.'): # a not valid reference call because we done assign "." if the variant was not valid
           # print("filtered reference call: ", vcf_dict[variant][1])
            filtering_dict[vcf_dict[variant][1][-1]][1] +=1
            ref_calls[vcf_dict[variant][1][-1]][0] +=1
            if vcf_dict[variant][1][2].split(",")[0] =="notPass":
               ref_calls[vcf_dict[variant][1][-1]][1] += 1
            elif vcf_dict[variant][1][2].split(",")[0] =="lowReadDepth":
                 ref_calls[vcf_dict[variant][1][-1]][2] += 1
            elif vcf_dict[variant][1][2].split(",")[0] =="lowAllelefreq":    
                 ref_calls[vcf_dict[variant][1][-1]][3] += 1
            elif vcf_dict[variant][1][2].split(",")[0] =="unplaced":    
                 ref_calls[vcf_dict[variant][1][-1]][3] += 1                             
        else: #zero is for reference and variants that don't pass filtering
            #calls not passing criteria (information which one in the dictionary, used to assign it to the step where it was filtered out)
            filter_criterium=vcf_dict[variant][1][1].split(",")[0]
            if filter_criterium == "lowAllelefreq":
                filtering_dict[vcf_dict[variant][1][2]][4] +=1
                vcf_dict[variant][1][0]=("/")
            elif filter_criterium == "lowReadDepth":
                filtering_dict[vcf_dict[variant][1][2]][3] +=1
                vcf_dict[variant][1][0]=("/")
            elif filter_criterium == "notPass":
                filtering_dict[vcf_dict[variant][1][2]][2] +=1
                vcf_dict[variant][1][0]=("/")
            elif filter_criterium == "ref_call":
                #print(vcf_dict[variant])
                filtering_dict[vcf_dict[variant][1][-1]][1] +=1
                ref_calls[vcf_dict[variant][1][-1]][0] +=1
            elif filter_criterium == "unplaced":
                filtering_dict[vcf_dict[variant][1][2]][5] +=1
                vcf_dict[variant][1][0]=("/")

    #total starting point of filtering (all, non-ref + valid calls)
    total=[filtering_dict["SNP"][0], 
           filtering_dict["sINS"][0], 
           filtering_dict["INS"][0],
           filtering_dict["sDEL"][0], 
           filtering_dict["DEL"][0], 
           filtering_dict["MIXED"][0]]
    
    #total variant calls (total-"0" calls):  
    variant_calls=[total[0]-filtering_dict["SNP"][1], 
                   total[1]-filtering_dict["sINS"][1], 
                   total[2]-filtering_dict["INS"][1],
                   total[3]-filtering_dict["sDEL"][1], 
                   total[4]-filtering_dict["DEL"][1], 
                   total[5]-filtering_dict["MIXED"][1]]
    
    #variants filtered out because they dont "PASS" vg internal filtering
    passed=[variant_calls[0]-filtering_dict["SNP"][2], variant_calls[1]-filtering_dict["sINS"][2], variant_calls[2]-filtering_dict["INS"][2],variant_calls[3]-filtering_dict["sDEL"][2], variant_calls[4]-filtering_dict["DEL"][2], variant_calls[5]-filtering_dict["MIXED"][2]]
    
    #variants filtered out, because of readdepth < 7
    readdepth=[passed[0]-filtering_dict["SNP"][3], passed[1]-filtering_dict["sINS"][3], passed[2]-filtering_dict["INS"][3],passed[3]-filtering_dict["sDEL"][3], passed[4]-filtering_dict["DEL"][3], passed[5]-filtering_dict["MIXED"][3]]
    
    #variants filtered out that have allelefreq <80%
    allelefreq=[readdepth[0]-filtering_dict["SNP"][4], readdepth[1]-filtering_dict["sINS"][4], readdepth[2]-filtering_dict["INS"][4],readdepth[3]-filtering_dict["sDEL"][4], readdepth[4]-filtering_dict["DEL"][4], readdepth[5]-filtering_dict["MIXED"][4]]
    
    #variants filtered out which are on unplaced sequences
    placed=[allelefreq[0]-filtering_dict["SNP"][5], allelefreq[1]-filtering_dict["sINS"][5], allelefreq[2]-filtering_dict["INS"][5],allelefreq[3]-filtering_dict["sDEL"][5], allelefreq[4]-filtering_dict["DEL"][5], allelefreq[5]-filtering_dict["MIXED"][5]]

    data=[total,variant_calls,passed,readdepth,allelefreq, placed]
   
    df=pd.DataFrame(data, index=["total calls", "variants", "pass", "readdepth >= 7", "allele freq >= 80%", "placed sequence"], columns=["SNP", "sINS", "INS", "sDEL", "DEL", "MIXED"])
    
    with open (stats_outfile, "a") as stats_out:
        print(df, file=stats_out)  
        
    #reference filtered out because they dont "PASS" vg internal filtering
    ref=[ref_calls["SNP"][0], ref_calls["sINS"][0], ref_calls["INS"][0],ref_calls["sDEL"][0], ref_calls["DEL"][0], ref_calls["MIXED"][0]]
    passed_ref=[ref_calls["SNP"][0]-ref_calls["SNP"][1], ref_calls["sINS"][0]-ref_calls["sINS"][1], ref_calls["INS"][0]-ref_calls["INS"][1],ref_calls["sDEL"][0]-ref_calls["sDEL"][1], ref_calls["DEL"][0]-ref_calls["DEL"][1], ref_calls["MIXED"][0]-ref_calls["MIXED"][1]]
    
    #reference filtered out, because of readdepth < 7
    readdepth_ref=[passed_ref[0]-ref_calls["SNP"][2], passed_ref[1]-ref_calls["sINS"][2], passed_ref[2]-ref_calls["INS"][2],passed_ref[3]-ref_calls["sDEL"][2], passed_ref[4]-ref_calls["DEL"][2], passed_ref[5]-ref_calls["MIXED"][2]]
    
    #reference filtered out that have allelefreq <80%
    allelefreq_ref=[readdepth_ref[0]-ref_calls["SNP"][3], readdepth_ref[1]-ref_calls["sINS"][3], readdepth_ref[2]-ref_calls["INS"][3],readdepth_ref[3]-ref_calls["sDEL"][3], readdepth_ref[4]-ref_calls["DEL"][3], readdepth_ref[5]-ref_calls["MIXED"][3]]
    
    #reference filtered out which are on unplaced sequences
    placed_ref=[allelefreq_ref[0]-ref_calls["SNP"][4], allelefreq_ref[1]-ref_calls["sINS"][4], allelefreq_ref[2]-ref_calls["INS"][4],allelefreq_ref[3]-ref_calls["sDEL"][4], allelefreq_ref[4]-ref_calls["DEL"][4], allelefreq_ref[5]-ref_calls["MIXED"][4]]

    data_ref=[ref,passed_ref,readdepth_ref,allelefreq_ref, placed_ref]
   
    df=pd.DataFrame(data_ref, index=["Reference", "pass", "readdepth >= 7", "allele freq >= 80%", "placed sequence"], columns=["SNP", "sINS", "INS", "sDEL", "DEL", "MIXED"])
    with open (stats_outfile, "a") as stats_out:
        print("Filtering of reference variants from the graph", file=stats_out)
        print(df, file=stats_out)      
    return vcf_dict
#Method for deviding small Indels and large Indels in INS/DEL and sINS/sDEL
def categorise_call(variant):
    vartype=variant[5].split(";")[1]
    REF_len=len(variant[1])
    ALT_len=len(variant[2])    
    if vartype=="INS":
        if ALT_len-1 <50:
            vartype="sINS"
    elif vartype=="DEL":
        if (REF_len-1 <50):
            vartype="sDEL"    
    return vartype
#additional counting method; counts how many Variants of each type are in the file in total and devided by valid and not valid 
#extract q scores and readDepth that could be used for plotting 
def count_calls(vcf_dict, stats_out):
    categories = {
        'total': [0, 0, 0],
        'sINS': [0, 0, 0],
        'INS': [0, 0, 0],
        'sDEL': [0, 0, 0],
        'DEL': [0, 0, 0],
        'SNP': [0, 0, 0],
        "MIXED":[0,0,0],
    }
    
    total_qscores=[]
    total_readDepth=[]
    
    for key in vcf_dict:
        validInfo=str(vcf_dict[key][1][0])
        #validInfo=int(validInfo)
        vartype=vcf_dict[key][1][-1]
        
        categories['total'][0]+=1
        if validInfo=="1":
            total_qscores.append(float(vcf_dict[key][0][3]))
            total_readDepth.append(int(vcf_dict[key][0][7].split(":")[1]))
            categories['total'][1]+=1
        else:
            categories['total'][2]+=1
        
        category = categories[vartype]
        category[0] += 1
        category[1 if validInfo == "1" else 2] += 1
    with open (stats_out, "a")as stats_file:
        for category, counts in categories.items():
            print(f"{category} calls: {counts[0]} ({counts[1]}, {counts[2]})", file=stats_file)

            
def write_total(total_out, infile, vg_valid):
    valid_keys=[]
    for k in vg_valid.keys():
            valid_keys.append(k)
            
    with open(infile, "r") as vg_input, open(total_out, "w") as out:
        for line in vg_input:
            if line.startswith("#"):
                out.write(line)
            else:
                splitline=line.strip().split("\t")
                chrom = splitline[0]
                pos = splitline[1]
                variant=(chrom,pos)
                if variant in valid_keys:
                    if ("," in splitline[4]):
                        splitline=get_alt_allele_total(splitline, vg_valid[variant])
                    if (vg_valid[variant][1][0]=="." or vg_valid[variant][1][0]==("/")):
                        #print(vg_valid[variant])
                        col9=splitline[9].split(":")
                        col9[0]="."
                        splitline[9]=":".join(col9)

                    mod_line="\t".join(splitline)
                    out.write(mod_line)
                    out.write("\n")

def get_alt_allele_total(splitline, dict_entry):
    if ("," in splitline[4]):
        #print(splitline)
        splitline[4]=dict_entry[0][2] #ALT allele
        #modify gt column
        col9=splitline[9].split(":")
        gt=col9[0] #which gt is called, information to extract the suitable information for this alt allele
        #print(gt)
        if(gt=="0"):
            allele_depth=(col9[2].split(","))
            allele_depth=[int(s) for s in allele_depth]
            idx=allele_depth.index(max(allele_depth[1:]))
            AD=[col9[2].split(",")[0],col9[2].split(",")[idx]] #reference calls should have the max alt allele depth as second entry (or 0 if all reads agree on 0 gt)
            GL=[col9[3].split(",")[0], col9[3].split(",")[int(gt)]]
            col9[0]="0" #gt always 1
        else:
            AD=[col9[2].split(",")[0],col9[2].split(",")[int(gt)]] #if genotype is not 0 get the alt allele depth at this postion 
            GL=[col9[3].split(",")[0], col9[3].split(",")[int(gt)]]
            col9[0]="1" #gt always 1
            
            
        col9[2]=",".join(AD)
        col9[3]=",".join(GL)
        splitline[9]=":".join(col9)
        
        col7=splitline[7].split(";")
        var_info=[col7[3].split("=")[0],"=",col7[3].split("=")[1].split(",")[int(gt)-1]] #0=Vartype of first alt allle, 1=Vartype of second altallele
        var_info="".join(var_info)
        col7[3]=var_info
        splitline[7]=";".join(col7)
        #print(splitline)
        
    return splitline

#Method for writing valid variants in a VCF file 
def write_valid_vcf(outfile, mod_outfile, vg_valid, infile):
    valid_keys=[]
    for k in vg_valid.keys():
        if vg_valid[k][1][0]==1:
            valid_keys.append(k)
    
    with open (infile, "r") as vg_input:
        with open (outfile, "w") as valid_output: #output without modifying the vcf file (keep all alt alleles)
            with open (mod_outfile, "w") as mod_output: #output modified vcf file, only the called alt allele is kept and GT is set to 1 
                for line in vg_input:
                    if line.startswith("#"):
                        valid_output.write(line)
                        mod_output.write(line)
                    else:
                        splitline=line.strip().split("\t")
                        chrom = splitline[0]
                        pos = splitline[1]
                        variant=(chrom,pos)
                        if variant in valid_keys:
                            if ("," in splitline[4]):
                                splitline_mod=get_alt_allele(splitline, vg_valid[variant])
                                mod_line="\t".join(splitline_mod)
                                mod_output.write(mod_line)
                                mod_output.write("\n")
                            else:
                                mod_output.write(line)
                                
                            valid_output.write(line)
                            
def get_alt_allele(splitline, dict_entry):
    if ("," in splitline[4]):
        #modify gt column, keep only information for called gt and set it to 1 
        #remove information on not called gt alleles, keep information about reference and called gt
        splitline[4]=dict_entry[0][2]
        
        col9=splitline[9].split(":")
        gt=col9[0]
        refAD=col9[2].split(",")[0]
        gtAD=col9[2].split(",")[int(gt)]
        col9[2]=refAD + "," + gtAD
        refGL=col9[3].split(",")[0]
        gtGL=col9[3].split(",")[int(gt)]
        col9[3]=refGL + "," + gtGL
        col9[0]="1"
        splitline[9]=":".join(col9)
        
        '''
        #modify gt column
        col9=splitline[9].split(":")
        gt=col9[0] #which gt is called, information to extract the suitable information for this alt allele

        AD=[col9[2].split(",")[0],col9[2].split(",")[int(gt)]] #0=ref, 1= first alt allele ...
        col9[2]=",".join(AD)
        
        GL=[col9[3].split(",")[0], col9[3].split(",")[int(gt)]]
        col9[3]=",".join(GL)
        
        col9[0]="1" #gt always 1
        splitline[9]=":".join(col9)
       '''
       #modify Vartype information in Info column, keep only information for called gt
        col7=splitline[7].split(";")
        var_info=[col7[3].split("=")[0],"=",col7[3].split("=")[1].split(",")[int(gt)-1]]
        var_info="".join(var_info)
        col7[3]=var_info
        splitline[7]=";".join(col7)

    return splitline
                         
#Method for searching a specific variant in the dictionary        
def search_entry(chrom, position, vg_valid_dictionary):
    for key in vg_valid_dictionary.keys():
        if (chrom, position) == key:
            print(vg_valid_dictionary[key])


#input and output files 
vg_infile=sys.argv[1]
#"/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2015T_graph_Alignment/2015T_variants_MQ30_BQ20_vartype.vcf"
#
agp_file=sys.argv[2]
vg_outfile=sys.argv[3]
mod_vg_outfile=sys.argv[4]
vg_filter_stats=sys.argv[5]
vg_total_out=sys.argv[6]
excl_contig=sys.argv[7]
depth= int(sys.argv[8])
alfreq=int(sys.argv[9])
vg_dict=read_vcf(vg_infile)
unplaced_seq=read_agp(agp_file,excl_contig)
vg_dict_valid=validate_calls(vg_dict,vg_filter_stats ,unplaced_seq)
count_calls(vg_dict_valid, vg_filter_stats)
write_valid_vcf(vg_outfile, mod_vg_outfile, vg_dict_valid, vg_infile)
write_total(vg_total_out, vg_infile,vg_dict_valid)
#search_entry("JACEHA010000005.1", "600987", vg_dict_valid)
