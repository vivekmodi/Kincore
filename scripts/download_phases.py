#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 20:12:22 2020

@author: vivekmodi
"""
import subprocess,os

def download_phases(pwd,df):
    print('Downloading phase files from PDB...')
    
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4].lower()
        if not os.path.isfile(f'{pwd}/mtz_files/{pdb}.mtz'):
            if not os.path.isfile(f'{pwd}/mtz_files/{pdb}.mtz.gz'):
                cmd=(f'wget -P {pwd}/mtz_files http://edmaps.rcsb.org/coefficients/{pdb}.mtz')    #Phenix requires unzipped mtz and .cif files
                subprocess.call(cmd,shell=True)
            
        if not os.path.isfile(f'{pwd}/mtz_files/{pdb}.cif'):    #copy cif files from kinasecifs folder to here
            if not os.path.isfile(f'{pwd}/mtz_files/{pdb}.cif.gz'):
                cmd=f'cp {pwd}/kinasecifs/{pdb}.cif.gz {pwd}/mtz_files/{pdb}.cif.gz;gunzip {pwd}/mtz_files/{pdb}.cif.gz'
                subprocess.call(cmd,shell=True)
