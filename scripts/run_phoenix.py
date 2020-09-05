#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 21:25:46 2020

@author: vivekmodi
"""
import subprocess,os

def run_phoenix(pwd,df):
    print('Running Phoenix...')
    log=open(f'{pwd}/kinasepml.log','a')
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4].lower()
        
        if not os.path.isfile(f'{pwd}/mtz_files/{pdb}.mtz'):
            log.write(f'run_phoenix: file not found {pdb}.mtz')
            continue
        else:    
            if not os.path.isfile(f'{pwd}/mtz_files/{pdb}_2mFo-DFc_map.ccp4'):
                if not os.path.isfile(f'{pwd}/mtz_files/{pdb}_2mFo-DFc_map.ccp4.gz'):
                    cmd=f'phenix.maps {pwd}/mtz_files/{pdb}.cif {pwd}/mtz_files/{pdb}.mtz'
                    subprocess.call(cmd,shell=True)
                    
    log.close()
                   
                
        