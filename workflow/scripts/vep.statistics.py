#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 00:05:58 2024

@author: yomna
"""

import pandas as pd
import sys

def parse_vep_output(filename,segment,sample_names):
    parse_data = False  
    consequences_all = {}  
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(segment):
                parse_data = True
                continue 
            elif parse_data and line.strip() == '':  # A new line break would suggest a segment close.
                parse_data = False
                continue
            elif parse_data:
                line_parts = line.strip().split('\t')  
                if len(line_parts) >= 2:  # Expecting each row to have at least one quantity.
                    # Categorizing or earmarking the primary form.
                    cons_type, count = line_parts[0], line_parts[1]
                    consequences_all[cons_type] = count
                # Here, a replication or reshaping could go on for your mentioned "Coding" section.

    # Here, we attempt the stratagem or resource for the feed-into your ask; again, it'll key to the pattern material.
    # Crafting the feed into DataFrame for 'Consequences (all)':
    df_consequences_all = pd.DataFrame(list(consequences_all.items()), columns=['Variant Type', sample_names])

    return df_consequences_all

# Employ the tooling with the defined access or road to the statement's response.
file_name = sys.argv[1]
sample_names=sys.argv[2]
df_consequences = parse_vep_output(file_name,'[Consequences (all)]',sample_names)
coding=parse_vep_output(file_name,'[Coding consequences]',sample_names)
df_consequences.to_csv(sys.argv[3], index=False)
coding.to_csv(sys.argv[4], index=False)