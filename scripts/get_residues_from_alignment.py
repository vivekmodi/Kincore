#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 13:10:08 2020

@author: vivekmodi
"""
import pandas as pd
pwd='/home/vivekmodi/Applications/Flask/Kinases'
df=pd.read_csv(f'{pwd}/Kinases_df-2020-06-25.csv',sep='\t')
def get_residues_from_alignment(pwd,df):
    print('Getting residue numbers from alignment...')
    df_align=pd.read_csv(f'{pwd}/All-organisms-alignment-residue-corresspondence.tab',sep=' ')  #make sure this file has a header
    
    gtk=425;hinge1=426;dfg_asp=1338
    
    for i in df.index:
        uniprotid=df.at[i,'UniprotID']
        domain=df.at[i,'Domain']
         
        try:
            df.at[i,'GTKnum']=int(df_align[(df_align.Domain==domain) & (df_align.Uniprot==uniprotid) & (df_align.AlignNum.astype(int)==gtk)].UniNum)
            df.at[i,'GTKres']=df_align[(df_align.Domain==domain) & (df_align.Uniprot==uniprotid) & (df_align.AlignNum.astype(int)==gtk)].ResType.to_string(index=False).strip()
            df.at[i,'Hinge1']=int(df_align[(df_align.Domain==domain) & (df_align.Uniprot==uniprotid) & (df_align.AlignNum.astype(int)==hinge1)].UniNum)
            df.at[i,'DFG_Aspres']=df_align[(df_align.Domain==domain) & (df_align.Uniprot==uniprotid) & (df_align.AlignNum.astype(int)==dfg_asp)].ResType.to_string(index=False).strip()
        except:
            print(df.at[i,'PDBid'],df.at[i,'UniprotID'])
            
     #       fhandle_alignment.seek(0)
     #       for lines in fhandle_alignment:
     #           lines=lines.strip();lines=lines.split()
        
     #           if domain==lines[0] and gtk==int(lines[3]):
     #               df.at[i,'GTKnum']=int(lines[2])
     #               df.at[i,'GTKres']=lines[1]
     #           if domain==lines[0] and hinge1==int(lines[3]):
     #               df.at[i,'Hinge1']=int(lines[2])
     #           if domain==lines[0] and dfg_asp==int(lines[3]):
     #               df.at[i,'DFG_Aspres']=lines[1]            #Asp res type, the same entry is also filled in motifs_for_nonhuman.py
     #               break
                    
    #fhandle_alignment.close()
    return df


                