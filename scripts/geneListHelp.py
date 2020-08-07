#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 08:21:50 2020

@author: vivekmodi
"""

def geneListHelp(pwd,df):
    print('Creating gene list help file...')
    fhandle_genelist=open(f'{pwd}/static/geneListHelpFile.tab','w')
    fhandle_genelist.write('List of Organism, HGNC gene, uniprot id and protein names with a known structure and included in the database\n')
    uniprotid_list=list()
    df_sorted=df.sort_values(['Specie','Gene']).copy()
    for i in df_sorted.index:
        gene_name=df_sorted.at[i,'Gene']
        organism=df_sorted.at[i,'Specie']
        uniprotid=df_sorted.at[i,'UniprotID']
        protein_name=df_sorted.at[i,'Protein']
        
        if uniprotid not in uniprotid_list:
            fhandle_genelist.write(f'{organism}\t{gene_name}\t{uniprotid}\t{protein_name}\n')
            uniprotid_list.append(uniprotid)
        