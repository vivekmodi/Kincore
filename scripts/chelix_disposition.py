#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:00:57 2020

@author: vivekmodi
"""

def chelix_disposition(pwd,df):
    print('Identifying C-helix conformation...')
    for i in df.index:
        if float(df.at[i,'Lys_Glu'])<=10.0:
            df.at[i,'C-helix']='in'
        if float(df.at[i,'Lys_Glu'])>10.0:
            df.at[i,'C-helix']='out'
        if float(df.at[i,'Lys_Glu'])==999:
            df.at[i,'C-helix']='None'
            
    return df