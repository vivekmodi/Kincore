#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 08:02:41 2020

@author: vivekmodi
"""
from Bio import PDB
import gzip
import sys
import pandas as pd

def compute_distance_from_rre4(structure,ligandname,ligandid,rre4num):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    min_rre4=999
    for model in structure:
        for chain in model:
            for residue1 in chain:
                if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                    for atom1 in residue1:
                        if atom1.element!='H':

                            for residue2 in chain:
                                if residue2.get_id()[1]==rre4num:    #Only rre4num
                                    ignoremodified.seek(0)
                                    if residue2.get_id()[0]==' ' or ((residue2.id[0][2:]+'\n') in ignoremodified.readlines()):  #Only protein atoms or modified residues
                                        for atom2 in residue2:
                                            if atom2.element!='H' and atom2.fullname not in ('CA','O','N','C'):   #Side chain contact
                                                distance=residue1[atom1.fullname]-residue2[atom2.fullname]
                                                if min_rre4>distance:
                                                    min_rre4=distance
    return min_rre4

def compute_distance_from_hinge(structure,ligandname,ligandid,hinge1):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    min_hinge=999;
    for model in structure:
        for chain in model:
            for residue1 in chain:
                if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                    for atom1 in residue1:
                        if atom1.element!='H':

                            for residue2 in chain:
                                if residue2.get_id()[1]>=hinge1 and residue2.get_id()[1]<=hinge1+2:        #compute distance from hinge
                                    ignoremodified.seek(0)
                                    if residue2.get_id()[0]==' ' or ((residue2.id[0][2:]+'\n') in ignoremodified.readlines()):
                                        for atom2 in residue2:
                                            if atom2.element!='H' and (atom2.fullname=='O' or atom2.fullname=='N'):    #Main chain contact
                                                distance=residue1[atom1.fullname]-residue2[atom2.fullname]
                                                if min_hinge>distance:
                                                    min_hinge=distance
    return min_hinge

def compute_distance_from_pocket_residues(structure,ligandname,ligandid):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    dfgcontact=0;dfgoutcontact=0;distance=dict()
    contact_list=list();contact_list_front=list()
    backpocket_count=dict();frontpocket_count=dict()
    backpocket_count[ligandname+':'+ligandid]=0
    frontpocket_count[ligandname+':'+ligandid]=0

    for model in structure:
        for chain in model:
            for residue1 in chain:
                if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                    for atom1 in residue1:
                        if atom1.element!='H':
                            for residue2 in chain:
                                ignoremodified.seek(0)
                                if residue2.get_id()[0]==' ' or ((residue2.id[0][2:]+'\n') in ignoremodified.readlines()):    #Only protein atoms
                                    res_id=residue2.get_id()[1]
                                    if res_id in range(106,185) or res_id in range(187,196)  or res_id in range(420,423) or res_id in range(1337,1340) or res_id==1011 or res_id==959 :
                                        for atom2 in residue2:
                                            if atom2.element!='H':
                                                distance[res_id]=round(float((residue1[atom1.fullname]-residue2[atom2.fullname])),2)
                                                
                                                #Contact with Type2 pocket present
                                                if res_id in (149,153,959,1011) and distance[res_id]<=4.5:  
                                                        dfgoutcontact+=1
                                                
                                                #Contact with X-D mainchain
                                                if distance[res_id]<=4 and res_id in range(1337,1339) and (atom2.fullname=='O' or atom2.fullname=='N') and dfgcontact==0 and res_id not in contact_list: 
                                                    dfgcontact=1
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)
                                                
                                                #Contact with Phe sidechain
                                                elif distance[res_id]<=4 and res_id==1339 and res_id not in contact_list and (atom2.fullname!='O' and atom2.fullname!='N' and atom2.fullname!='CA'):
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)
                                                
                                                #Contact with non-XDF backpocket residues
                                                elif distance[res_id]<=4 and res_id not in range(1337,1340) and res_id not in contact_list:
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)

                                                    
                                                    #Contact with N-ter of C-helix (Frontpocket)
                                                    if res_id in range(106,144) and res_id not in contact_list_front:
                                                        frontpocket_count[ligandname+':'+ligandid]+=1
                                                        contact_list_front.append(res_id)
                                                        
    return frontpocket_count, backpocket_count, dfgoutcontact, distance

def correct_chain_diff_in_ligand_type_labels(df):   #If two chains in the same PDB have Type1 and Type1.5 labels then keep only Type1.5 for both the chains
    for i in df.index:
        pdb1=df.at[i,'PDBid'][0:4]
        chain1=df.at[i,'PDBid'][4:]
        ligand_label1=df.at[i,'Ligand_label']
        ligand_name1=df.at[i,'Ligand']
        
        
        for j in df.index:
            pdb2=df.at[j,'PDBid'][0:4]
            chain2=df.at[j,'PDBid'][4:]
            ligand_label2=df.at[j,'Ligand_label']
            ligand_name2=df.at[j,'Ligand']
            
            if pdb1==pdb2 and chain1!=chain2:
                if ',' in ligand_name1:
                    for position1,ligand_n1 in enumerate(ligand_name1.split(',')):
                        for position2,ligand_n2 in enumerate(ligand_name2.split(',')):
                            if ligand_n1.split(':')[0]==ligand_n2.split(':')[0]:
                                if 'Type1'==ligand_label1.split(',')[position1] and 'Type1.5' in ligand_label2.split(',')[position2]:
                                    df.at[i,'Ligand_label']=ligand_label2
                else:
                    if ligand_name1.split(':')[0]==ligand_name2.split(':')[0]:
                        if 'Type1'==ligand_label1 and 'Type1.5' in ligand_label2:
                            df.at[i,'Ligand_label']=ligand_label2
    return df
                    
                    
def classify_ligands(pwd,df):
    fhandle_output=open(f'Ligand_classification.tab','w')
    print('Classifying different kinds of ligands...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        uniprotid=df.at[i,'UniprotID']
        spatial=df.at[i,'Spatial']
        dihedral=df.at[i,'Dihedral']
        rre4num=150;hinge1=426
        ligand_label=list()

        if 'No_ligand' in df.at[i,'Ligand']:
            df.at[i,'Ligand_label']='No_ligand'
            continue
        else:
            df.at[i,'Ligand_label']='None'    #make default label None

        try:
            handle=gzip.open(f'{pwd}/kinasechains_renumber_alignment/{pdbs}.cif.gz','rt')
        except:
            log=open(f'{pwd}/kinasepml.log','a')
            log.write(f'classify_ligands: File not found {pdbs}.cif.gz\n')
            continue

        parser=PDB.MMCIFParser()
        structure=parser.get_structure(pdbs,handle)

        ligandlist=df.at[i,'Ligand'].split(',')
        for items in ligandlist:
            ligandname=items.split(':')[0]
            ligandid=items.split(':')[1]



            #Identify allosteric ligands
            min_rre4=compute_distance_from_rre4(structure,ligandname,ligandid,rre4num)
            min_hinge=compute_distance_from_hinge(structure,ligandname,ligandid,hinge1)

            #Contacts with pocket residues
            (frontpocket_count, backpocket_count, dfgoutcontact,distance)=compute_distance_from_pocket_residues(structure,ligandname,ligandid)
           

            if min_rre4!=999 or min_hinge!=999:
                if min_rre4>=6.5 and min_hinge>=6.5:
                    ligand_label.append('Allosteric')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tAllosteric\n')
                elif min_hinge>=6 and backpocket_count[ligandname+':'+ligandid]>=3:
                    ligand_label.append('Type3')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tType3\n')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and frontpocket_count[ligandname+':'+ligandid]==0 and dfgoutcontact==0:
                    ligand_label.append('Type1.5_Back')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tType1.5_Back\n')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and frontpocket_count[ligandname+':'+ligandid]>=1 and dfgoutcontact==0:
                    ligand_label.append('Type1.5_Front')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tType1.5_Front\n')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and dfgoutcontact>=1 and spatial=='DFGout':
                    ligand_label.append('Type2')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tType2\n')
                else:
                    ligand_label.append('Type1')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tType1\n')


        df.at[i,'Ligand_label']=','.join(ligand_label)
        
        
    fhandle_output.close()      #This file does not contain updated labels
    
    df=correct_chain_diff_in_ligand_type_labels(df)
    return df

if __name__=='__main__':
    pwd='/home/vivekmodi/Applications/Flask/Kinases'
    filename=sys.argv[1]
    df=pd.read_csv(filename,sep='\t',header='infer')
    classify_ligands(pwd,df)
    