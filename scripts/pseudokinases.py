#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:11:58 2020

@author: vivekmodi
"""

def identify_pseudokinases(df):
    print('Identifying pseudokinases...')
    
    for i in df.index:
       
        alk=str(df.at[i,'ALKres']);hrd=str(df.at[i,'HRDres']);dfg_asp=str(df.at[i,'DFG_Aspres'])
            
        if alk!='K' or hrd!='D' or dfg_asp!='D':
            df.at[i,'Status']='Pseudo'
        else:
            df.at[i,'Status']='Typical'
    
    
    return df