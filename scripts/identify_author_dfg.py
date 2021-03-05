#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:55:19 2020

@author: vivekmodi
"""

import sys, os, gzip, re
import pandas as pd

def identify_author_dfg(pwd,conf_df):


    for i in conf_df.index:
     
        pdbid=conf_df.at[i,'PDBid'];pdb=pdbid[0:4].lower()
        chainid=pdbid[4:];uni_dfgnum=conf_df.at[i,'DFGnum'];uni_aspnum=int(uni_dfgnum)-1

        sifts_handle=gzip.open(f'{pwd}/kinasesifts/{pdb}.csv.gz','rt')
        sifts_handle.seek(0)
        sifts_handle.readline()

        for lines in sifts_handle:
            lines=lines.strip();lines=lines.split(',')
            sifts_aspnum=lines[5]

            if chainid==lines[4]:
                if re.search('[A-Za-z]',sifts_aspnum):     #To remove any insertion codes
                    sifts_aspnum=sifts_aspnum[:-1]

                if int(uni_aspnum)==int(sifts_aspnum):
                    author_aspnum=lines[2]     #When residue number is not present in pdb (e.g. in cases with residue not resolved), this value is printed as null
                    author_aspres=lines[3]
                    conf_df.at[i,'Author_Aspnum']=author_aspnum
                    conf_df.at[i,'Author_Aspres']=author_aspres
                    break


    return conf_df

if __name__=='__main__':
    filename=sys.argv[1]
    pwd=os.getcwd()
    conf_df=pd.read_csv(filename,sep='\t',header='infer')
    conf_df=identify_author_dfg(pwd,conf_df)
    conf_df.to_csv(f'{filename}-new',sep='\t',index=False)
