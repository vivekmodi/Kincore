#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 19:36:30 2020

@author: vivekmodi
"""
import gzip
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def read_header(pwd,kinasecifs,pdb_domain_dict):
    print('Reading header from MMCIF...')
    pdb_reso_dict=dict()
    fhandle_reso_read=open((pwd+'/Resolution.tab'),'r')
    fhandle_reso_append=open((pwd+'/Resolution.tab'),'a')       #Open the file again to append new values
    
    for pdbs in pdb_domain_dict:
        fhandle_reso_read.seek(0)
        pdb_present=0
        for lines in fhandle_reso_read:
            lines=lines.strip();lines=lines.split()
            if pdbs in lines[0]:
                pdb_reso_dict[pdbs]=lines[1]
                pdb_present=1
                break
        
        if pdb_present==0:
            handle=gzip.open((kinasecifs+'/'+pdbs[0:4].lower()+'.cif.gz'),'rt')
            mmcif_dict=MMCIF2Dict(handle)
            resolution=999
            if '_refine.ls_d_res_high' in mmcif_dict:
                try:
                    resolution=mmcif_dict['_refine.ls_d_res_high']      #Other items could be obtained in the same way
                    if resolution=='?':
                        resolution=999
                    if resolution=='.':
                        resolution=999
                except:
                    #if resolution=='?':
                    resolution=999
                #print(pdbs,resolution)
                if type(resolution)==list:      #Some cif files have two values for resolution, like 5ZN0A
                    resolution=float(resolution[0])
                resolution=float(resolution)
                pdb_reso_dict[pdbs]=np.round(resolution,2)
                fhandle_reso_append.write(f'{pdbs} {pdb_reso_dict[pdbs]}\n')            
            else:
                resolution=999
                pdb_reso_dict[pdbs]=resolution
                fhandle_reso_append.write(f'{pdbs} {pdb_reso_dict[pdbs]}\n')
            handle.close()
    fhandle_reso_append.close()
    fhandle_reso_read.close()
    return pdb_reso_dict

#if __name__ == '__main__':
#    input_dir=sys.argv[1]
#    pdb=sys.argv[2]
#    resolution=read_header(input_dir,pdb)
#    print(resolution)