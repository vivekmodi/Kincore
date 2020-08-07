#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:11:58 2020

@author: vivekmodi
"""
import pandas as pd

def identify_pseudokinases(df):
    print('Identifying pseudokinases...')
    #fhandle_list=open('Pseudokinase-list','r')
    

    for i in df.index:
        #domain=df.at[i,'Domain']
        #uniprotid=df.at[i,'UniprotID']
        
        #if 'HUMAN' in uniprotid: 
        #    pseudo=0
        #    fhandle_list.seek(0)
        #    for lines in fhandle_list:
        #        lines=lines.split()
                
        #        if lines[1]==domain:
        #            df.at[i,'Status']='Pseudo'
        #            pseudo=1
        #    if pseudo==0:
        #        df.at[i,'Status']='Typical'
        
        #else:
        alk=str(df.at[i,'ALKres']);hrd=str(df.at[i,'HRDres']);dfg_asp=str(df.at[i,'DFG_Aspres'])
            
        if alk!='K' or hrd!='D' or dfg_asp!='D':
            df.at[i,'Status']='Pseudo'
        else:
            df.at[i,'Status']='Typical'
    
    #fhandle_list.close()
    return df