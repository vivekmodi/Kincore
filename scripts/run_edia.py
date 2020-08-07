#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 13:34:31 2020

@author: vivekmodi
"""
import os, subprocess

def run_edia(pwd,df):
    print('Running EDIA...')
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4].lower()
        
        if not os.path.isfile(f'{pwd}/mtz_files/Edia_out/{pdb}atomscores.csv'):
            cmd=f'Ediascorer -t {pwd}/mtz_files/{pdb}.cif -d {pwd}/mtz_files/{pdb}_2mFo-DFc_map.ccp4 -o {pwd}/mtz_files/Edia_out/'
            subprocess.call(cmd,shell=True)