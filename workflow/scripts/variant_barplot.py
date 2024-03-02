#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 22:03:22 2024

@author: yomna
"""
import pandas as pd
import matplotlib.pyplot as plt
import sys
import random
# Define a color palette with 100 colors
color_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                 '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
                 '#393b79', '#ff9c00', '#5254a3', '#d32f2f', '#669900', '#ffcc00', '#0099c6', '#dd4477', '#66aa00', '#b82e2e',
                 '#316395', '#994499', '#22aa99', '#aaaa11', '#6633cc', '#e67300', '#8b0707', '#651067', '#329262', '#5574a6',
                 '#3b3eac', '#b77322', '#16d620', '#b91383', '#f4359e', '#9c5935', '#a9c413', '#2a778d', '#668d1c', '#bea413',
                 '#0c5922', '#743411', '#31383b', '#505233', '#fc461d', '#73007d', '#f2e000', '#0e84b6', '#a35122', '#ec8025',
                 '#c983c8', '#d4adde', '#d3b683', '#f781bf', '#74c476', '#a695c7', '#c75d7b', '#afafaf', '#ad6f3b', '#dd3497',
                 '#7a6a83', '#cc33cc', '#669966', '#cc99cc', '#c0c0c0', '#ffd700', '#ff00ff', '#800080', '#ff6666', '#f78181',
                 '#81f781', '#39ff14', '#f6546a', '#b200b2', '#ba55d3', '#4682b4', '#006400', '#800000', '#663399', '#dc143c',
                 '#ffdab9', '#f08080', '#ff4500', '#da70d6', '#eee8aa', '#b0e0e6', '#4682b4', '#7b68ee', '#ff00ff', '#800080',
                 '#ffc0cb', '#87cefa', '#32cd32', '#66cdaa', '#20b2aa', '#778899']

# Shuffle the color palette to add some randomness
random.shuffle(color_palette)
def assign_variant_colors(variant_types):
    # Create a dictionary to store variant types and their assigned colors
    variant_colors = {}    
    # Assign colors to variant types
    for i, variant_type in enumerate(variant_types):
        # If we run out of colors, cycle through the palette
        if variant_type not in variant_colors:
           color_index = i % len(color_palette)
           variant_colors[variant_type] = color_palette[color_index]  
    return variant_colors

def combine_files_and_plot(file_names, output):
    # Define a colormap with more colors
    plotted_variants = {}

    # Plotting
    plt.figure(figsize=(10, 6))
    
    for file_name in file_names:
        cl = file_name.split("/")[-2].split("_")[0]  # Extract column label from file name
        # Read the file into a DataFrame
        df = pd.read_csv(file_name)
        df = df.rename(columns={"count": cl})
        df.set_index('Variant Type', inplace=True)
        # Convert counts to percentages
        df = df.div(df.sum(axis=0), axis=1) * 100
        # Sort variant types by the sum of counts across samples
        df_sorted = df.loc[df.sum(axis=1).sort_values(ascending=False).index]
        variant_colors = assign_variant_colors(df_sorted.index)
        for i, sample in enumerate(df_sorted.index):
            if sample not in plotted_variants:
                plt.bar(df_sorted.columns, df_sorted.loc[sample], label=sample, color=variant_colors.get(sample, 'black'))
                plotted_variants[sample] = True
            else:
                plt.bar(df_sorted.columns, df_sorted.loc[sample], color=variant_colors.get(sample, 'black'))

    plt.xlabel('Sample')
    plt.ylabel('Percentage')
    plt.title('Combined Bar Plot for Each Sample (Percentage)')
    plt.xticks(rotation=90)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save the plot as a PDF file
    plt.savefig(output)
    
    # Show the plot in a separate window
    #plt.show()

# Example usage
output = sys.argv[1]
#'/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/INS.variant_consequences.all.pdf'
file_names = sys.argv[2:]
#['/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2000B_graph_Alignment/INS.variant_consequences.all.txt', '/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2020T_graph_Alignment/INS.variant_consequences.all.txt', '/home/yomna/Desktop/PhD_Yomna_Gohar/graph_genome_project/ToxoVar/analysis/Graph/graph_construction/2015T_graph_Alignment/INS.variant_consequences.all.txt']
combine_files_and_plot(file_names, output)


