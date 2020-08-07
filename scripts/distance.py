#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
from Bio import PDB
#homedir='/home/vivekmodi/Applications/Flask/Kinases/'



def compute_distance(pwd,pdbfilename,res1,res2,res1_type,res2_type,setting):
    if setting=='rre4-phe':
        phe_atom_type=identify_phe_atom(res2_type)
        dis_phe_rre4=distance_atoms(pwd,pdbfilename,res1,res2,'CA',phe_atom_type)      #Change atom names for other residue types
        return dis_phe_rre4

    if setting=='lys-phe':
        phe_atom_type=identify_phe_atom(res2_type)
        dis_phe_lys=distance_atoms(pwd,pdbfilename,res1,res2,'CA',phe_atom_type)      #Change atom names for other residue types
        return dis_phe_lys

    if setting=='saltbridge':
        dis_sb=distance_atoms(pwd,pdbfilename,res1,res2,'CB','CB')      #Change atom names for other residue types
        return dis_sb

def identify_phe_atom(res2_type):
    if res2_type=='F':
        return 'CZ'
    if res2_type=='L' or res2_type=='P':
        return 'CG'
    if res2_type=='M':
        return 'CE'
    if res2_type=='S':
        return 'OG'
    if res2_type=='V':
        return 'CB'
    if res2_type=='W':
        return 'CZ3'
    if res2_type=='Y':
        return 'OH'
    if res2_type=='A':
        return 'CB'

def distance_atoms(pwd,pdbfilename,res1,res2,atm1,atm2):
    if '.pdb' in pdbfilename.lower():
        parser=PDB.PDBParser()
    if '.cif' in pdbfilename.lower():
        parser=PDB.MMCIFParser()
    structure=parser.get_structure('PDB',(pwd+'/server/uploads/'+pdbfilename))
    atom_present=0; res1=int(res1); res2=int(res2)
    print(res1,res2,atm1,atm2)
    for model in structure:
        for chain in model:
            for residue in chain:
                if int(residue.id[1])==int(res1) and residue.get_id()[0]==' ':
                    if residue.has_id(atm1):
                        #print(res1,residue.id,chain[res1])
                        residue1=chain[res1]
                        atom_present=atom_present+1
                if int(residue.id[1])==int(res2) and residue.get_id()[0]==' ':
                    if residue.has_id(atm2):
                        residue2=chain[res2]
                        atom_present=atom_present+1

    if atom_present==2:
        distance=np.round((residue1[atm1]-residue2[atm2]),2)
        #print(distance)
        #print(f'{distance:.2f}')
        return distance
    else:
        return 999

#Usage from commandline
#python distance.py 1GAGA.pdb 1078 1178 MET PHE rre4-phe
#python distance.py 1GAGA.pdb 1057 1074 LYS GLU saltbridge

#if __name__ == '__main__':
#    pdbfilename=sys.argv[1];res1=sys.argv[2];res2=sys.argv[3];res1_type=sys.argv[4];res2_type=sys.argv[5];setting=sys.argv[6]
#    distance=compute_distance(pwd,pdbfilename,res1,res2,res1_type,res2_type,setting)
#    print(distance)
