#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:00:43 2020

@author: vivekmodi
"""
from Bio import PDB
import gzip
import pandas as pd

def get_release_date(kinasecifs,df):
    print('Get release date from .cif files...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        filename=pdbs[0:4].lower()+'.cif.gz'
        handle=gzip.open(f'{kinasecifs}/{filename}','rt')
        parser=PDB.MMCIFParser(QUIET=True)
        structure=parser.get_structure(pdbs[0:4],handle)
        df.at[i,'Date_time']=pd.to_datetime(structure.header['deposition_date'])
    df['Date']=df['Date_time'].dt.date
    df.drop(['Date_time'],axis=1)
    return df
