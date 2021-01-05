#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:57:22 2020

@author: vivekmodi
"""
import gzip, sys
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd

#Did not use this script because mmcif's record for _chem_comp.pdbx_synonyms is not consistent. In some files there are spaces between names, in some semi colon
#and in some files chemical name is included as common name. Therefore, parsing it is impossible.

def create_ligand_syn_dict(pwd,df):
    ligSynDict=dict();pdblist=list();synonym_list=list()
    for i in df.index:
        if 'No_ligand' in df.at[i,'Ligand']:
            continue
        
        pdbs=df.at[i,'PDBid']
        if pdbs[0:4].lower() in pdblist:
            continue
        pdblist.append(pdbs[0:4].lower())
        filename=pdbs[0:4].lower()+'.cif.gz'
        handle=gzip.open(f'{pwd}/kinasecifs/{filename}','rt')
        mmcif_dict = MMCIF2Dict(handle)
        
        ligandlist=df.at[i,'Ligand'].split(',')
        for items in ligandlist:
            ligandname=items.split(':')[0]
            
            
            for count,compound in enumerate(mmcif_dict['_chem_comp.id']):
                if compound==ligandname:
                    if mmcif_dict['_chem_comp.pdbx_synonyms'][count]!='?':
                        synonym_list=mmcif_dict['_chem_comp.pdbx_synonyms'][count].split(',')
                        for items in synonym_list:
                            ligSynDict[items]=ligandname
                            print(pdbs,items,ligandname)
        if i==30:
            break
                    
                    
if __name__ == '__main__':
    pwd='/home/vivekmodi/Applications/Flask/Kinases'
    filename=sys.argv[1]
    df=pd.read_csv(filename,sep='\t',header='infer')
    create_ligand_syn_dict(pwd,df)