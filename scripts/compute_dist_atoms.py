#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:03:40 2020

@author: vivekmodi
"""
#Usage from commandline python compute_dist_atoms.py /input_dir 1XYZA res1 atm1 res2 atm2
import gzip,sys
#import numpy as np
from Bio import PDB

def compute_distance(input_dir,filename,res1,atm1,res2,atm2):
    handle=gzip.open((input_dir+'/'+filename+'.cif.gz'),'rt')
    parser=PDB.MMCIFParser(QUIET=True)
    structure=parser.get_structure("PDB",handle)
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    atom_present=0
    for model in structure:
        for chain in model:
            for residue in chain:
                ignoremodified.seek(0)
                if (int(residue.id[1])==int(res1) and residue.get_id()[0]==' ') or (int(residue.id[1])==int(res1) and ((residue.id[0][2:]+'\n') in ignoremodified.readlines())):
                    if residue.has_id(atm1):
                        res1_object=residue
                        #residue1=chain[res1]
                        atom_present=atom_present+1
                ignoremodified.seek(0)
                if (int(residue.id[1])==int(res2) and residue.get_id()[0]==' ') or (int(residue.id[1])==int(res2) and ((residue.id[0][2:]+'\n') in ignoremodified.readlines())):
                    if residue.has_id(atm2):
                        res2_object=residue
                        #residue2=chain[res2]
                        atom_present=atom_present+1

    if atom_present==2:
        distance=round(float(res1_object[atm1]-res2_object[atm2]),2)    #round works only on type float, so first convert the number to float
        #print(distance)
        #print(f'{distance:.2f}')
        return distance
    else:
        return 999

if __name__ == '__main__':
    input_dir=sys.argv[1];filename=sys.argv[2];res1=sys.argv[3];atm1=sys.argv[4];res2=sys.argv[5];atm2=sys.argv[6]
    distance=compute_distance(input_dir,filename,res1,atm1,res2,atm2)
    print(distance)
