#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:55:19 2020

@author: vivekmodi
"""

import sys, os, gzip, re
import pandas as pd

def identify_author_dfg(pwd,conf_df):
    #conf_df=pd.read_csv(filename,sep='\t',header='infer')
    
    for i in conf_df.index:
        pdbid=conf_df.at[i,'PDBid'];pdb=pdbid[0:4].lower()
        chainid=pdbid[4];uni_dfgnum=conf_df.at[i,'DFGnum']
        
        sifts_handle=gzip.open(f'{pwd}/kinasesifts/{pdb}.csv.gz','rt')
        sifts_handle.seek(0)
        sifts_handle.readline()
        
        for lines in sifts_handle:
            lines=lines.strip();lines=lines.split(',')
            sifts_dfgnum=lines[5]
            
            if chainid==lines[4]:
                if re.search('[A-Za-z]',sifts_dfgnum):     #To remove any insertion codes
                    sifts_dfgnum=sifts_dfgnum[:-1]
                    
                if int(uni_dfgnum)==int(sifts_dfgnum):
                    author_dfgnum=lines[2]     #When residue number is not present in pdb (e.g. in cases with residue not resolved), this value is printed as null
                    conf_df.at[i,'Author_DFGnum']=author_dfgnum
                    break
       
        
    return conf_df
    
#if __name__=='__main_':
#    filename=sys.argv[1]
#    pwd=os.getcwd()
#    conf_df=identify_author_dfg(pwd,filename)
#    conf_df.to_csv(f'{filename}-new',sep='\t',index=False)
    