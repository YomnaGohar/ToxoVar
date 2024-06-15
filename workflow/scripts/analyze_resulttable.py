#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 10:23:08 2024

@author: marie
"""
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
#script to produce the results from the resulttable: 
combined_table="/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/results/merged_vg_combined_table_placed_ref.txt"
df = pd.read_csv(combined_table, sep='\t', dtype="str")

def update_vartype(row):
    var=row["ALT"].split(",")
    ref=row["Ref"].split(",")

    if row["Vartype"] == "INS":

        if len(var[0]) < 50:

            return "sINS"

        else:

            return "INS"

    elif row["Vartype"] == "DEL":

        if len(ref[0]) < 50:

            return "sDEL"

        else:

            return "DEL"

    elif row["Vartype"] == "SNV":

        return "SNP"

    else:

        return row["Vartype"]



df["Vartype"] = df.apply(update_vartype, axis=1)



def exclude_unplaced_variants(df):

    total_rows=len(df)

    placed_sequences_df = df[~(df[['2015T_merged', '2020T_merged', '2000B_merged', '2015T_vg', '2020T_vg', '2000B_vg']] == '.').all(axis=1)]

    total_placed=len(placed_sequences_df)

    return total_rows, placed_sequences_df, total_placed



def count_total(df_placed):

    columns_to_analyze = ['2015T_merged', '2020T_merged', '2000B_merged', '2015T_vg', '2020T_vg', '2000B_vg']

    # Dictionary to store counts for each column

    total_counts = {}



    # Iterate over columns and calculate counts

    for column in columns_to_analyze:

        total_counts[column] = {

            '0': (df_placed[column] == '0').sum(),

            '.': ((df_placed[column] == '.') | (df_placed[column] == '-')).sum(),

            '1, 2, or 3': ((df_placed[column] == '1') | (df_placed[column] == '2') | (df_placed[column] == '3')).sum()

        }

    

    total_counts_df = pd.DataFrame.from_dict(total_counts, orient='index')

    return total_counts_df

    



def count_vartype(df_placed):

    counts_per_vartype = {}

    vartypes=["SNP", "sINS", "INS", "sDEL", "DEL"]

    columns_to_analyze = ['2015T_merged', '2020T_merged', '2000B_merged', '2015T_vg', '2020T_vg', '2000B_vg']

    # Iterate over vartypes and calculate counts for each column

    for vartype in vartypes:

        counts_per_vartype[vartype] = {}

        for column in columns_to_analyze:

            counts_per_vartype[vartype][column] = {

                '0': ((df_placed['Vartype'] == vartype) & (df_placed[column] == '0')).sum(),

                '.': ((df_placed['Vartype'] == vartype) & (df_placed[column] == '.')| (df_placed[column] == '-')).sum(),

                '1, 2, or 3': ((df_placed['Vartype'] == vartype) & (df_placed[column].isin(['1', '2', '3']))).sum()

                }



    counts_vartype_df = pd.DataFrame.from_dict({(vartype, column): counts_per_vartype[vartype][column] 

                                    for vartype in counts_per_vartype.keys() 

                                    for column in counts_per_vartype[vartype].keys()},

                                    orient='index')

    return(counts_vartype_df)

    

# Apply the function to update vartype column

df["Vartype"] = df.apply(update_vartype, axis=1)

total_rows, df_placed, total_placed=exclude_unplaced_variants(df) 

total_counts_df=count_total(df_placed)

total_vartype_counts_df=count_vartype(df_placed)



print("total_rows:", total_rows)

print("total_placed:", total_placed)

print(total_counts_df)

print(total_vartype_counts_df)



#####################################

idx_pairs=[(5,8),(6,9),(7,10)]

sample=["2015T", "2020T","2000B"]



#result dictionary for allele table

categories_table_allele = {

    'both_identical':{"total":0,"2015T": 0, "2020T": 0, "2000B":0}, 

    'both_disagree': {"total":0,"2015T": 0, "2020T": 0, "2000B":0},

    'graph_only_called_allele': {"total":0,"2015T": 0, "2020T": 0, "2000B":0},

    'linear_only_called_allele': {"total":0,"2015T": 0, "2020T": 0, "2000B":0},

    'both_not': {"total":0,"2015T": 0, "2020T": 0, "2000B":0}

}



#methods for the allele table 

def count_total_sample_wise(df, idx_pair, categories, sample):  

    for index, row in df.iterrows():

        merged_value = row[idx_pair[0]]

        graph_value = row[idx_pair[1]]

        # Check conditions and update categories dictionary accordingly

        if str(merged_value) == str(graph_value):

            if merged_value != '.':

                categories['both_identical']["total"] += 1

                categories['both_identical'][sample] += 1

                

            else:

                categories['both_not']["total"] += 1

                categories['both_not'][sample] += 1

                

        elif str(merged_value) != '.' and str(graph_value) != '.':

            categories['both_disagree']["total"] += 1

            categories['both_disagree'][sample] += 1

           

        elif str(merged_value) != '.' and str(graph_value) == '.':

            categories['linear_only_called_allele']["total"] += 1

            categories['linear_only_called_allele'][sample] += 1

           

            

            #print(merged_value, graph_value, "merged_notgraph")

        elif str(merged_value) == '.' and str(graph_value) != '.':

            categories['graph_only_called_allele']["total"] += 1

            categories['graph_only_called_allele'][sample] += 1

            

        #elif str(merged_value) == '.' and str(graph_value) == '.':

        #    categories['both_not'] += 1

        #    print(merged_value, graph_value, "both_not")

        

    return categories

###################################

def plot_pie_chart(ax, categories, sample_name):

    labels = categories.keys()

    sizes = [categories[category][sample_name] for category in labels]

    explode = (0, 0, 0, 0, 0)  # explode the 1st slice



    wedges, texts, autotexts = ax.pie(sizes, explode=explode, autopct='%1.1f%%', startangle=140)

    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    ax.set_title(f'Sample-wise allele comparison \n {sample_name}')

    plt.setp(autotexts, size=12, weight="bold")



    

for i in range(0,3):

    #print(sample[i], idx_pairs[i])

    categories_table_allele=count_total_sample_wise(df_placed, idx_pairs[i], categories_table_allele, sample[i])



print(categories_table_allele)

# Create subplots

fig, axs = plt.subplots(1, 3, figsize=(15, 5))



legend_artists = []

for i, (ax, s) in enumerate(zip(axs, sample)):

    plot_pie_chart(ax, categories_table_allele, s)

    # Collect artists for the legend from the first subplot

    if i == 0:

        legend_artists += ax.patches



# Create a legend using artists from the first subplot

plt.legend(legend_artists, categories_table_allele.keys(), title="Categories", loc="center left", bbox_to_anchor=(1, 0.5))



plt.tight_layout()

plt.show()

########################
#Analyze sdv: 

categories_table_sda = {
    'both_nosdA':0, 
    'both_sdA': 0,
    'graph_sdA':0,
    'merged_sdA':0
    }
rows_table_sda={
    'both_no_sdv':[], 
    'both_sdv':[],
    'graph_sdv':[],
    'merged_sdv':[]
}

def remove_dot(allele_list):
    allele_list = [item for item in allele_list if item != '.']
    return allele_list
        
def sdv_code_list(code_list, vg, merged):
    if (merged not in code_list):
        code_list.append(merged)
    if (vg not in code_list):
        code_list.append(vg)   
    return code_list
    
def count_total_sdA(df, result_dict, rows_dict):
    code_list=[]
    for index,row in df.iterrows():
        update_vartype(row)
        merged=[row[5], row[6], row[7]]        
        vg=[row[8], row[9], row[10]]      
        #print(merged)
        merged_l=list(set(merged))
        vg_l=list(set(vg))
        #print(merged, vg)        
        merged_r=remove_dot(merged_l)
        vg_r=remove_dot(vg_l)        
        if len(merged_r) <=1 and len(vg_r) <=1:
            result_dict["both_nosdA"] +=1
            rows_table_sda["both_no_sdv"].append((merged, vg,len(row[3]), row[4]))
        elif len(merged_r) >1 and len(vg_r) >1:
            result_dict["both_sdA"] += 1
            rows_table_sda["both_sdv"].append((merged, vg, len(row[3]),row[4]))
            sdv_code_list(code_list, vg, merged)
        elif (len(merged_r) <= 1 and len(vg_r) > 1):
            result_dict["graph_sdA"] +=1
            rows_table_sda["graph_sdv"].append((vg,len(row[3]), row[4]))
            if (vg not in code_list):
                code_list.append(vg)
            #print(row, "graph_sdA")
        elif (len(merged_r) > 1 and len(vg_r) <= 1):
            result_dict["merged_sdA"] += 1
            rows_table_sda["merged_sdv"].append((merged,len(row[3]), row[4]))
            if (merged not in code_list):
                code_list.append(merged)
        else:
            print("not assigned to any categorie")        
    return result_dict, rows_dict, code_list

def which_code(code_set):
    #code_joined="".join(code_set)
    sample_1=code_set[0]
    sample_2=code_set[1]
    sample_3=code_set[2]
    if  "." not in code_set and "-" not in code_set and len(set(code_set)) < 3:
        if sample_1 not in code_set[1:3]:
           return "2015T"
        elif sample_2 not in [sample_1,sample_3]:
             return "2020T"
        elif sample_3 not in [sample_1,sample_2]:
             return "2000B"
    elif sample_1 == "." or sample_1 == "-":
        return "2020T_or_2000B"
    elif sample_2 == "." or sample_2 == "-":
        return "2015T_or_2000B"
    elif sample_3 == "." or sample_3 == "-":
        return "2015T_or_2020T"
    else:
        return "sdv_all"
    #else:
        #print("missed key: ", code_joined)
def sdv_which_variant(rows_dict, code_list):
    sda_categorie_vartype_sample_codes={
        "both_sdv":{ 
            "both_sdv_agree": {
                "SNP":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                "sINS": {"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                "sDEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                "INS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                "DEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0},
                "Total":{"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}
                },
            "both_sdv_disagree":{ 
                "merged":{
                    "SNP":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "sINS": {"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "sDEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "INS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "DEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0},
                    "Total":{"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}
                    }, 
                "vg": {
                    "SNP":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "sINS": {"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "sDEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "INS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
                    "DEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0},
                    "Total":{"2015T":0, "2020T":0, "2000B":0,"2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}
                    }
                }
        }, 
        "graph_sdv":{ 
            "SNP":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "sINS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "sDEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "INS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "DEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0},
            "Total":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}
        }, 
        "merged_sdv":{ 
            "SNP":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "sINS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "sDEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "INS": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}, 
            "DEL": {"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0},
            "Total":{"2015T":0, "2020T":0, "2000B":0, "2015T_or_2020T":0, "2020T_or_2000B":0, "2015T_or_2000B":0, "sdv_all":0}
        }
    }
    
    #extract information from the dictionary, including the codes of merged and/or vg and the vartype, of all three categories with sdv:
        #both_sdv: include the variants were both approaches agree on the fact that this variant contains sda, but they do not necessaryly agree on the code or the sample/mix of sample containing sda
        #graph_sdv: include variants were only the graph calls the position as containing sdv
        #merged_sdv: include variants were only the merged called the positions as contiaing sdv
        
        
    for key in rows_table_sda.keys():
        if key == "both_sdv":      #results between merged and vg can either are on the code or disagree
            for i in range(0, len(rows_table_sda[key])):    #entries of the three categories
                if rows_table_sda[key][i][0] == rows_table_sda[key][i][1]:  #codes agree between merged and vg calling, categorized and counted
                
                    cat=which_code(rows_table_sda[key][i][0])
                    sda_categorie_vartype_sample_codes[key]["both_sdv_agree"][rows_table_sda[key][i][-1]][cat] +=1    #counting per vartype
                    sda_categorie_vartype_sample_codes[key]["both_sdv_agree"]["Total"][cat] +=1                       #coutnting in total 
                else:        #codes disagree between merged and vg calling
                    cat_m=which_code(rows_table_sda[key][i][0])                             #category of merged 
                    cat_vg= cat=which_code(rows_table_sda[key][i][1])                       #category of vg calling 
                    
                    #if categories disagree: categorize merged information:  
                    sda_categorie_vartype_sample_codes[key]["both_sdv_disagree"]["merged"][rows_table_sda[key][i][-1]][cat_m] +=1
                    sda_categorie_vartype_sample_codes[key]["both_sdv_disagree"]["merged"]["Total"][cat_m] +=1
                    
                    #if categories disagree: categorize vg information:  
                    sda_categorie_vartype_sample_codes[key]["both_sdv_disagree"]["vg"][rows_table_sda[key][i][-1]][cat_vg] +=1
                    sda_categorie_vartype_sample_codes[key]["both_sdv_disagree"]["vg"]["Total"][cat_vg] +=1
                
        elif key=="graph_sdv" or key=="merged_sdv":    #categorizing the entries of the other two categories, where only one approach called the variant as sdv
            for i in range(0, len(rows_table_sda[key])):
                cat=which_code(rows_table_sda[key][i][0])
                if cat: 
                    sda_categorie_vartype_sample_codes[key][rows_table_sda[key][i][-1]][cat] +=1
                    sda_categorie_vartype_sample_codes[key]["Total"][cat] +=1
  
    return sda_categorie_vartype_sample_codes   

def make_result_df_count_category(result_dict):
    data_df={}
    for k in result_dict.keys():
        if k=="both_sdv":
            for k2 in result_dict[k].keys():
                if k2 == "both_sdv_disagree":
                    for k3 in result_dict[k][k2].keys(): 
                        row_name=f"{k2}_{k3}"
                        row_t=result_dict[k][k2][k3]["Total"]
                        data_df[row_name]=row_t
                       # print(k2, k3, row_t)
                        
                else:
                    row_name=k2
                    row_t=result_dict[k][k2]["Total"]
                    data_df[row_name]=row_t
                    #print(k2, row_t)
        else:
            row_name=k
            row_t=result_dict[k]["Total"]
        #    print(k, row_t)
        
            data_df[row_name]=row_t
    #print(data_df)
    
    columns = ["categorie","2015T", "2020T", "2000B", "2015T_or_2020T", "2020T_or_2000B", "2015T_or_2000B", "sdv_all"]
 
    rows = [{"categorie": "codes", "2015T": "100 011 211" , "2020T": "010, 101, 121", "2000B": "001, 110", "2015T_or_2020T": "01., 10., 12., 12-, 01-" , "2020T_or_2000B": ".10, .01, .21", "2015T_or_2000B": "2.1, 0.1, 1.0", "sdv_all": "231, 120, 021, 201"}]
    
    # Populate the rows list with data
    for key, values in data_df.items():
        row = {"categorie": key}
        row.update(values)
        rows.append(row)
    
    # Create the DataFrame
    df_res = pd.DataFrame(rows, columns=columns)
    
    # Print the DataFrame
    #print(df_res)
    return df_res


# Define function to create stacked bar plot
def create_stacked_barplot(data):
    titles = {
       'both_sdv_agree': 'linear-reference-based and  \n graph-based approach agree on SDV alleles',
       'both_sdv_disagree_merged': 'SDV alleles disagreeing \n corresponding to linear-refernence-based calling',
       'both_sdv_disagree_vg': 'SDV alleles disagreeing \n corresponding to graph-based genotyping',
       'merged_sdv': 'SDV called only by \n linear-reference-based approach',
       'graph_sdv': 'SDV called only by \n graph-based approach'
   }
    def plot_variations(variations, title):
        variation_types = list(variations.keys())
        variation_types =variation_types[0:-1]
        # Extract all years from the nested dictionaries
        years = set()
        for variation in variations.values():
            for year in variation.keys():
                years.add(year)

        # Initialize plot
        fig, ax = plt.subplots()

        # Initialize bottom to 0 for stacking
        bottom = np.zeros(len(variation_types))
        bar_width = 0.35

        # Color map for different years
        colors = plt.cm.get_cmap("tab20", len(years))
        non_zero_years=[]
        for i, year in enumerate(sorted(years)):
            values = []
            for variation in variation_types:
                values.append(variations[variation].get(year, 0))

            if any(values):
                ax.bar(variation_types, values, bar_width, bottom=bottom, label=year, color=colors(i))
                non_zero_years.append(year)
            bottom += values


        ax.set_xlabel('Variation Types')
        ax.set_ylabel('Counts')
        ax.set_title(title)
        ax.legend(title= "SDV category according to \n allele combination", labels=[year for year in non_zero_years])

        plt.show()

    for key, variations in data.items():
        if key == 'both_sdv':
            for subkey, subvariations in variations.items():
                if(subkey == "both_sdv_disagree"):
                    for subsubkey, subsubvariations in subvariations.items():
                        plot_variations(subsubvariations, titles[f'{subkey}_{subsubkey}'])
                else: 
                    plot_variations(subvariations, titles[subkey])
        else:
            plot_variations(variations, titles[key])

# Call the function with sample data
#create_stacked_barplot(data)

# def create_stacked_barplot(inner_dic):
#     print(inner_dic)
#     for k2, variations in inner_dic.items():
#         variation_types = list(variations.keys())
#         #inner_values=set()
#         inner_cats= set(k for k in variations.keys())
        
#         # Initialize plot
#         fig, ax = plt.subplots()

#         # Initialize bottom to 0 for stacking
#         bottom = np.zeros(len(variation_types))
#         bar_width = 0.35

#         # Color map for different years
#         colors = plt.cm.get_cmap('tab20', len(inner_cats))

#         for i, inner_value in enumerate(sorted(inner_cats)):
#             values = [variations[v].get(inner_value, 0) for v in variation_types]
#             ax.bar(variation_types, values, bar_width, bottom=bottom, label=inner_value, color=colors(i))
#             bottom += values

#         ax.set_xlabel('Variation Types')
#         ax.set_ylabel('Counts')
#         ax.set_title(f'Stacked Bar Plot for {k2}')
#         ax.legend(title='Year')

#         plt.show()

# Call the function with sample data

categories_table_sda, rows_table_sda, code_list=count_total_sdA(df, categories_table_sda, rows_table_sda)
result_dict=sdv_which_variant(rows_table_sda, code_list)
df_res=make_result_df_count_category(result_dict)
create_stacked_barplot(result_dict)

#df_res.to_csv("/home/marie/masterthesis/results/benchmarking/exclusive_analysis/analysis_1804/manual_benchmarking/sda_analysis_categories.csv", index=False)


