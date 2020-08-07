#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 14:40:41 2020

@author: vivekmodi
"""
from Bio import PDB
import gzip

def min_distance_from_ligand(pwd,df):
    print('Computing minimum distance from ligand...')
    
    fhandle_output=open('Ligands-All.tab','w')
    fhandle_output.write('UniprotID\tPDBid\tSpatial\tDihedral\tLigand\tLigandID\tLig_RRE4\tLig_Hinge\tPhe_RRE4\tLig_RRE4CA\tPhe_Lys\tLig_label\n')
    
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        uniprotid=df.at[i,'UniprotID']
        spatial=df.at[i,'Spatial']
        dihedral=df.at[i,'Dihedral']
        phe_lys=df.at[i,'Phe_Lys']
        chelix=df.at[i,'C-helix']
        ligandlist=df.at[i,'Ligand']
        
        rre4num=df.at[i,'RREnum']+4
        hinge1num=df.at[i,'Hinge1']
        phe_rre4=df.at[i,'Phe_Glu4']
        
        ligand_rre4=list()
        ligand_hinge=list()      #minimum distances stored in a list for the structures with more than one ligand
        ligand_label=list()
        ligand_rre4CA=list()
        
        if 'NO_LIGAND' in ligandlist:
            ligand_rre4.append(999)
            ligand_hinge.append(999)
            df.at[i,'Lig_RRE4']=ligand_rre4
            df.at[i,'Lig_Hinge']=ligand_hinge
            df.at[i,'Lig_RRE4CA']=ligand_rre4CA.append(999)
            df.at[i,'Ligand_label']=ligand_label.append('None')
            continue
        
        handle=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz','rt')        
        parser=PDB.MMCIFParser()
        structure=parser.get_structure(pdbs,handle)
        
        handle2=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz','rt')        
        parser2=PDB.MMCIFParser()
        structure2=parser2.get_structure(pdbs,handle2)
        
        for ligands in ligandlist:
            ligand_name=ligands.split(':')[0]
            ligand_id=ligands.split(':')[1]
            
            min_rre4=999
            min_hinge=999
            min_rre4CA=999
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        #print(f'{ligand_name} {residue.id[0]}')
                        if residue.id[0]==('H_'+ligand_name) and residue.id[1]==int(ligand_id):
                            #print(f'{pdbs} {ligand_name}')
                            for atom in residue:
                            
                                for model2 in structure2:
                                    for chain2 in model2:
                                        for residue2 in chain2:
                                            if residue2.id[1]==rre4num:
                                                
                                                if residue2.has_id('CA'):
                                                    distance_rre4CA=round(float(residue[atom.fullname]-residue2['CA']),2)
                                                    if distance_rre4CA<min_rre4CA:
                                                        min_rre4CA=distance_rre4CA
                                                
                                                for atom2 in residue2:
                                                    distance_rre4=round(float(residue[atom.fullname]-residue2[atom2.fullname]),2)
                                                    
                                                    
                                                    if distance_rre4<min_rre4:
                                                        min_rre4=distance_rre4
                                                    
                                                        
                                
                                            if residue2.id[1]==hinge1num or residue2.id[1]==hinge1num+1 or residue2.id[1]==hinge1num+2:
                                                for atom2 in residue2:
                                                    if atom2.fullname=='O' or atom2.fullname=='N':
                                                        distance_hinge=round(float(residue[atom.fullname]-residue2[atom2.fullname]),2)
                                                    
                                                        if distance_hinge<min_hinge:
                                                            min_hinge=distance_hinge
                                                            
            ligand_rre4.append(str(min_rre4))
            ligand_rre4CA.append(str(min_rre4CA))
            ligand_hinge.append(str(min_hinge))      #converted into str because .join in the next step does not work in float
            
            #assign ligand label
            if min_rre4!=999 or min_hinge!=999:
                if min_rre4>=6.5 and min_hinge>=6.5:     #think about 5ORLA
                    ligand_label.append('Allosteric')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tAllosteric\n')
                    continue
                if min_rre4<=5 and min_hinge>=5:
                    ligand_label.append('Type3')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tType3\n')
                    continue
                if min_rre4<=4.5 and phe_rre4<=11.5:
                    ligand_label.append('Type1.5')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tType1.5\n')
                    continue
                if min_rre4>4.5 and phe_rre4<=11.5:
                    ligand_label.append('Type1')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tType1\n')
                    continue
                if min_rre4<=4.5 and phe_rre4>11.5:
                    ligand_label.append('Type2')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tType2\n')
                    continue
                if min_rre4>4.5 and phe_rre4>11.5:
                    ligand_label.append('Type1')
                    fhandle_output.write(f'{uniprotid}\t{pdbs}\t{spatial}\t{dihedral}\t{ligand_name}\t{ligand_id}\t{min_rre4}\t{min_hinge}\t{phe_rre4}\t{min_rre4CA}\t{phe_lys}\t{chelix}\tType1\n')
                    continue
                
            
            
        
        
        ligand_rre4=','.join(ligand_rre4)
        ligand_rre4CA=','.join(ligand_rre4CA)
        ligand_hinge=','.join(ligand_hinge)
        ligand_label=','.join(ligand_label)
        #print(f'{pdbs} {ligands} {ligand_rre4}')    
        df.at[i,'Lig_RRE4']=ligand_rre4
        df.at[i,'Lig_RRE4CA']=ligand_rre4CA
        df.at[i,'Lig_Hinge']=ligand_hinge
        df.at[i,'Ligand_label']=ligand_label
        

    fhandle_output.close()    
    return df