#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:13:52 2020

@author: vivekmodi
"""
import gzip
from Bio import PDB

def extract_ligands_webserver(pwd,pdbfilename,index,df,structure):
    omitligands=open(f'{pwd}/ListOrganicMoleculesNotLigands.txt','r')
    modified_aa=open(f'{pwd}/List_modified_aminoacid.txt','r')
    
    model_id=df.at[index,'Model_id']
    chain_id=df.at[index,'Chain_id'] 
    #domain_begin=df.at[index,'DomainBegin']    
    #domain_end=df.at[index,'DomainEnd']
        
    ligandname='';ligandpresent=0;ligandlist=list();
            
    for model in structure:
        for chain in model:
            if int(model.id)==int(model_id) and chain.id==chain_id:
                for residue in chain:
                        if 'H_' in residue.id[0]:
                            ligand_id=residue.id[1]
                            ligandname=str(residue.id[0][2:].strip(" "))
                            omitligands.seek(0)
                            if (ligandname+'\n') in omitligands.readlines():
                                continue                 #Modified residues are also included in this list
                            modified_aa.seek(0)
                            if (ligandname+'\n') in modified_aa.readlines():
                                continue
                            
                                
                            else:
                                ligandcount=0
                                for atom in residue:       #iterate over ligand atoms
                                    if ligandcount==0:     #if the same ligand is not counted before
                                        for residue2 in chain:   #iterate over protein residues
                                            if PDB.is_aa(residue2):
                                                for atom2 in residue2:
                                                    distance=residue[atom.fullname.strip()]-residue2[atom2.fullname.strip()]       #if ligand atom makes contact with residue within kinase domain; Some PDB files have spaces around atom name so use strip()
                                                    if distance<=4 and ligandcount==0:                             #to make sure atoms for the same ligand are not counted twice          
                                                        ligandpresent=1
                                                        ligandcount=1
                                                        ligandlist.append(f'{ligandname}:{ligand_id}')     #format ['ATP:300']
                                                        df.at[index,'Ligand']=','.join(ligandlist)
                                                        
                                                                        
                                print(df.at[index,'Ligand'])
                                #print(ligandlist)
                                #return ligandname
    if ligandpresent==0:
        ligandname='No_ligand'
        ligandlist.append(ligandname)
        df.at[index,'Ligand']=','.join(ligandlist)
            
    return df