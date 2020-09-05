#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:55:35 2020

@author: vivekmodi
"""

import gzip
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
    
    
def compute_all_distances(pwd,df):
    #Compute distance between residues
    kinasechains_renumber_uniprot=f'{pwd}/kinasechains_renumber_uniprot'
    
    print('Computing distance between residue pairs...')

    for i in df.index:
        pdbs=df.at[i,'PDBid']
        pdb_present=0
    #    fhandle_read_distances.seek(0)
        #for lines in fhandle_read_distances:
        #    lines=lines.strip();lines=lines.split()
        #    if pdbs in lines[0]:
        #        df.at[i,'Lys_Glu']=float(lines[1]);df.at[i,'Phe_Glu4']=float(lines[2]);df.at[i,'Phe_Lys']=float(lines[3])
        #        pdb_present=1
        #        break

        if pdb_present==0:
            #pdbid=pdbs[0:4]
            #chainid=pdbs[4]
            #domain_name=df.at[i,'Domain']
            df.at[i,'Lys_Glu']=999.0;df.at[i,'Phe_Glu4']=999.0;df.at[i,'Phe_Lys']=999.0

            df.at[i,'Lys_Glu']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'ALKnum']),'CB',int(df.at[i,'RREnum']),'CB')

            if df.at[i,'DFGres']=='F' or df.at[i,'DFGres']=='R':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CZ',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CZ',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='L' or df.at[i,'DFGres']=='P' or df.at[i,'DFGres']=='N':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CG',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CG',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='M':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CE',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CE',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='S':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'OG',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'OG',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='H':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'NE2',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'NE2',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='V' or df.at[i,'DFGres']=='A':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CB',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CB',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='W':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CZ3',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CZ3',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='Y':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'OH',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'OH',int(df.at[i,'ALKnum']),'CA')
            if df.at[i,'DFGres']=='G':
                df.at[i,'Phe_Glu4']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CA',int(df.at[i,'RREnum'])+4,'CA')
                df.at[i,'Phe_Lys']=compute_distance(kinasechains_renumber_uniprot,pdbs,int(df.at[i,'DFGnum']),'CA',int(df.at[i,'ALKnum']),'CA')

            #fhandle_append_distances.write(f"{pdbs} {df.at[i,'Lys_Glu']} {df.at[i,'Phe_Glu4']} {df.at[i,'Phe_Lys']}\n")


    #fhandle_read_distances.close()
    #fhandle_append_distances.close()
    return df