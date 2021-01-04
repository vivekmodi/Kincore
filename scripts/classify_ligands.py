#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 08:02:41 2020

@author: vivekmodi
"""
from Bio import PDB
import gzip

def compute_distance_from_rre4(structure,ligandname,ligandid,rre4num):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    min_rre4=999
    for chain1 in structure:
        for model1 in chain1:
            for residue1 in model1:
                if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                    for atom1 in residue1:
                        if atom1.element!='H':

                            for chain2 in structure:
                                for model2 in chain2:
                                    for residue2 in model2:
                                        #res_id=residue2.get_id()[1]
                                        ignoremodified.seek(0)
                                        if ((residue2.get_id()[0]==' ') or ((residue2.id[0][2:]+'\n') in ignoremodified.readlines())):
                                            if residue2.get_id()[1]==rre4num:    #compute distance from rre4
                                                for atom2 in residue2:
                                                    if atom2.element!='H':
                                                        if atom2.fullname!='CA' and atom2.fullname!='O' and atom2.fullname!='N' and atom2.fullname!='C':
                                                            distance=residue1[atom1.fullname]-residue2[atom2.fullname]
                                                            if min_rre4>distance:
                                                                min_rre4=distance
    return min_rre4

def compute_distance_from_hinge(structure,ligandname,ligandid,hinge1):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    min_hinge=999;
    for chain1 in structure:
                for model1 in chain1:
                    for residue1 in model1:
                        if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                            for atom1 in residue1:
                                if atom1.element!='H':

                                    for chain2 in structure:
                                        for model2 in chain2:
                                            for residue2 in model2:
                                                #res_id=residue2.get_id()[1]
                                                ignoremodified.seek(0)
                                                if residue2.get_id()[1]>=hinge1 and residue2.get_id()[1]<=hinge1+2:        #compute distance from hinge
                                                    for atom2 in residue2:
                                                        if atom2.element!='H':
                                                            if atom2.fullname=='O' or atom2.fullname=='N':
                                                                distance=residue1[atom1.fullname]-residue2[atom2.fullname]
                                                                if min_hinge>distance:
                                                                    min_hinge=distance
    return min_hinge

def compute_distance_from_pocket_residues(structure,ligandname,ligandid):
    ignoremodified=open(f'List_modified_aminoacid.txt','r')
    dfgcontact=0;dfgoutcontact=0;
    contact_list=list();contact_list_front=list()
    backpocket_count=dict();frontpocket_count=dict()
    backpocket_count[ligandname+':'+ligandid]=0
    frontpocket_count[ligandname+':'+ligandid]=0
    
    for chain1 in structure:
                for model1 in chain1:
                    for residue1 in model1:
                        if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                            for atom1 in residue1:
                                if atom1.element!='H':

                                    for chain2 in structure:
                                        for model2 in chain2:
                                            for residue2 in model2:
                                                res_id=residue2.get_id()[1]
                                                ignoremodified.seek(0)
                                                if (res_id>=106 and res_id<=184) or (res_id>=187 and res_id<=195)  or (res_id>=420 and res_id<=422) or \
                                                        res_id==1011 or res_id==959 or (res_id>=1337 and res_id<=1339):
                                                          
                                                            for atom2 in residue2:
                                                                if atom2.element!='H':
                                                                    if (res_id>=1337 and res_id<=1338) and (atom2.fullname!='O' and atom2.fullname!='N'):
                                                                        continue #only polar backbone contacts with XD
                                                                    if (res_id>=1337 and res_id<=1339) and (atom2.fullname=='O' or atom2.fullname=='N') and dfgcontact==1:
                                                                        continue

                                                                    distance=residue1[atom1.fullname]-residue2[atom2.fullname]
                                                                    if distance<=4 and res_id not in contact_list:
                                                                        distance=round(float(distance),2)
                                                                        backpocket_count[ligandname+':'+ligandid]+=1
                                                                        if res_id==1339 and (atom2.fullname=='O' or atom2.fullname=='N'):   #Always update contact list except for phe backbone
                                                                            pass
                                                                        else:
                                                                            contact_list.append(res_id)

                                                                        if (res_id>=1337 and res_id<=1339) and (atom2.fullname=='O' or atom2.fullname=='N'):  #contact with XD backbone present
                                                                            dfgcontact=1

                                                                        if ((res_id>=106 and res_id<=143)) and res_id not in contact_list_front:     #contact with front pocket present
                                                                            frontpocket_count[ligandname+':'+ligandid]+=1
                                                                            contact_list_front.append(res_id)

                                                                    if (res_id==1011 or res_id==959 or res_id==153 or res_id==149) and distance<=4.5:  #contact with Type2 pocket present
                                                                        dfgoutcontact+=1
                                                                        
    return frontpocket_count, backpocket_count, dfgoutcontact

def classify_ligands(pwd,df):
    fhandle_output=open(f'Ligand_classification.tab','w')
    print('Classifying different kinds of ligands...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        uniprotid=df.at[i,'UniprotID']
        spatial=df.at[i,'Spatial']
        dihedral=df.at[i,'Dihedral']
        rre4num=150;hinge1=426

        #ligand_rre4=list()
        #ligand_hinge=list()      #minimum distances stored in a list for the structures with more than one ligand
        ligand_label=list()

        ligand_label=list()
        


        if 'No_ligand' in df.at[i,'Ligand']:
            df.at[i,'Ligand_label']='No_ligand'
            continue
        else:
            df.at[i,'Ligand_label']='None'    #make default label None
            
        handle=gzip.open(f'{pwd}/kinasechains_renumber_alignment/{pdbs}.cif.gz','rt')
        parser=PDB.MMCIFParser()
        structure=parser.get_structure(pdbs,handle)

        #handle2=gzip.open(f'{pwd}/kinasechains_renumber_alignment/{pdbs}.cif.gz','rt')
        #parser2=PDB.MMCIFParser()
        #structure2=parser2.get_structure(pdbs,handle2)

        ligandlist=df.at[i,'Ligand'].split(',')
        for items in ligandlist:
            ligandname=items.split(':')[0]
            ligandid=items.split(':')[1]

            
            
            #Identify allosteric ligands
            min_rre4=compute_distance_from_rre4(structure,ligandname,ligandid,rre4num)
            min_hinge=compute_distance_from_hinge(structure,ligandname,ligandid,hinge1)
            
            #Contacts with pocket residues
            (frontpocket_count, backpocket_count, dfgoutcontact)=compute_distance_from_pocket_residues(structure,ligandname,ligandid)
            
            if min_rre4!=999 or min_hinge!=999:
                if min_rre4>=6.5 and min_hinge>=6.5:
                    ligand_label.append('Allosteric')
                    fhandle_output.write(f'{pdbs}\t{uniprotid}\t{spatial}\t{dihedral}\t{ligandname}\t{ligandid}\tAllosteric\n')
                elif min_hinge>=6 and backpocket_count[ligandname+':'+ligandid]>=3:     #min_hinge changed from 5 to 6 to correct TypeIII classification in DFGout, e.g 3II5A
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
    fhandle_output.close()
    return df
