#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:13:52 2020

@author: vivekmodi
"""
import gzip,re
from Bio import PDB

def extract_ligands(pwd,df):
    #kinasechains=f'{pwd}/kinasechains'
    fhandle_read_ligands=open((pwd+'/Ligands.tab'),'r')
    fhandle_append_ligands=open((pwd+'/Ligands.tab'),'a')
    omitligands=open(f'{pwd}/ListOrganicMoleculesNotLigands.txt','r')
    modified_aa=open(f'{pwd}/List_modified_aminoacid.txt','r')
    
    print('Extracting ligands...')
    #pdb_ligand_dict=dict()
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        domain_begin=df.at[i,'DomainBegin']    
        domain_end=df.at[i,'DomainEnd']
        fhandle_read_ligands.seek(0)
        pdb_present=0
        #for lines in fhandle_read_ligands:
        #    lines=lines.strip();lines=lines.split()
        #    if pdbs in lines[0]:
        #        df.at[i,'Ligand']=lines[1]
        #        pdb_present=1
        #        break
        
        if pdb_present==0:
            ligandname='';ligandpresent=0;ligandlist=list();
            parser=PDB.MMCIFParser(QUIET=True)
            handle=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz','rt')
            structure=parser.get_structure("PDB",handle)
            
            parser2=PDB.MMCIFParser(QUIET=True)
            handle2=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz','rt')
            structure2=parser2.get_structure("PDB",handle2)
    
            for model in structure:
                for chain in model:
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
                                        for model2 in structure2:      #iterate our protein residues
                                            for chain2 in model2:
                                                for residue2 in chain2:
                                                    if PDB.is_aa(residue2):
                                                        if int(residue2.id[1])>=int(domain_begin) & int(residue2.id[1])<=int(domain_end):
                                                            for atom2 in residue2:
                                                                distance=residue[atom.fullname]-residue2[atom2.fullname]       #if ligand atom makes contact with residue within kinase domain
                                                                if distance<=4 and ligandcount==0:                             #to make sure atoms for the same ligand are not counted twice          
                                                                    
                                                                    ligandpresent=1
                                                                    ligandcount=1
                                                                    ligandlist.append(f'{ligandname}:{ligand_id}')     #format ['ATP:300']
                                                                    if i==0:              #if its 0th element then pandas does not print list brackets []
                                                                        df.at[i,'Ligand']=','.join(ligandlist)
                                                                    else:
                                                                        df.at[i,'Ligand']=','.join(ligandlist)
                                                                    
                                #print(df.at[i,'Ligand'])
                                #print(ligandlist)
                                #return ligandname
            if ligandpresent==0:
                ligandname='No_ligand'
                ligandlist.append(ligandname)
                df.at[i,'Ligand']=','.join(ligandlist)
            
            #df.at[i,'Ligand']=sorted(set(df.loc[i,'Ligand']),key = str,reverse=True)      #To do alphanumeric sorting
            #pdb_ligand_dict[pdbs]=sorted(pdb_ligand_dict[pdbs])
            #df.at[i,'Ligand']=','.join(df.loc[i,'Ligand'])
            fhandle_append_ligands.write(f"{pdbs} {df.at[i,'Ligand']}\n")
            #print(pdbs,pdb_ligand_dict[pdbs])
    
    fhandle_append_ligands.close()
    fhandle_read_ligands.close()    
    return df