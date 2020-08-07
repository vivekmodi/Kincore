#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:21:12 2020

@author: vivekmodi
"""

def gene_synonym(pwd,df):
    #pdb_gene_synonym_dict=dict()
    fhandle_synonym=open((pwd+'/Uniprot-kinase-acc-gene-list.txt'),'r')
    for i in df.index:
        df.at[i,'Synonym']='None'
        for lines in fhandle_synonym:
            lines=lines.strip()
            lines=lines.split(';')
            
            if df.at[i,'Gene']==lines[2] and lines[3]!='':
                df.at[i,'Synonym']=lines[3]
        fhandle_synonym.seek(0)
    fhandle_synonym.close()
    return df