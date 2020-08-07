#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 00:46:37 2020

@author: vivekmodi
"""

def least_atoms_missing(loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,domain_name,spatial_label,dihedral_label):
    pdb_reso_sorted=sorted(pdb_reso_dict.items(),key=lambda x: x[1], reverse=True)      #Trick to sort dictionary in descending order
    min_missing_residues=999;min_missing_pdb='NA';min_reso=999
    
#    for pdbs,reso in pdb_reso_sorted:
#        if pdb_domain_dict[pdbs]==domain_name:
#             if pdb_dihedral_dict[pdbs]==dihedral_label:
#                 if loop_break_dict[pdbs]<=min_missing_residues:
#                     min_missing_residues=loop_break_dict[pdbs]
#                     min_missing_pdb=pdbs
#                     min_reso=reso
#    return (min_missing_pdb,min_missing_residues,min_reso)
    #dihedral_list=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus')
    for pdbs,reso in pdb_reso_sorted:
        if pdb_domain_dict[pdbs]==domain_name:
            if pdb_spatial_dict[pdbs]==spatial_label:
                    if pdb_dihedral_dict[pdbs]==dihedral_label:
                            if loop_break_dict[pdbs]<=min_missing_residues:
                                min_missing_residues=loop_break_dict[pdbs]
                                min_missing_pdb=pdbs
                                min_reso=reso
    return (min_missing_pdb,min_missing_residues,min_reso)


def representative_structure(pwd,loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict):
    print('Finding representative structure for each domain...')
    domain_list=sorted(set(pdb_domain_dict.values()))
    #dihedral_list=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus')
    spatial_list=('DFGin','DFGinter','DFGout')
    fhandle_output=open((pwd+'/Representative_structures.tab'),'w')
    for domain_name in domain_list:
        for spatial_label in spatial_list:
            dfgin=0;dfginter=0;dfgout=0
            if spatial_label=='DFGin':
                for dihedral_label in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing(loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,domain_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfgin=dfgin+1
                if dfgin==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                    
            if spatial_label=='DFGinter':
                for dihedral_label in ('BABtrans','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing(loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,domain_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfginter=dfginter+1
                if dfginter==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                    
            if spatial_label=='DFGout':
                for dihedral_label in ('BBAminus','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing(loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,domain_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfgout=dfgout+1
                if dfgout==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{domain_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
    fhandle_output.close()
    

def least_atoms_missing_ligand(loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,ligand_name,spatial_label,dihedral_label):
    pdb_reso_sorted=sorted(pdb_reso_dict.items(),key=lambda x: x[1], reverse=True)      #Trick to sort dictionary in descending order
    min_missing_residues=999;min_missing_pdb='NA';min_reso=999
    
#    for pdbs,reso in pdb_reso_sorted:
#        if pdb_domain_dict[pdbs]==domain_name:
#             if pdb_dihedral_dict[pdbs]==dihedral_label:
#                 if loop_break_dict[pdbs]<=min_missing_residues:
#                     min_missing_residues=loop_break_dict[pdbs]
#                     min_missing_pdb=pdbs
#                     min_reso=reso
#    return (min_missing_pdb,min_missing_residues,min_reso)
    #dihedral_list=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus')
    for pdbs,reso in pdb_reso_sorted:
        if ligand_name in pdb_ligand_dict[pdbs]:
            if pdb_spatial_dict[pdbs]==spatial_label:
                    if pdb_dihedral_dict[pdbs]==dihedral_label:
                            if loop_break_dict[pdbs]<=min_missing_residues:
                                min_missing_residues=loop_break_dict[pdbs]
                                min_missing_pdb=pdbs
                                min_reso=reso
    return (min_missing_pdb,min_missing_residues,min_reso)

def representative_structure_ligand(pwd,loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict):
    print('Finding representative structure for each ligand...')
    ligand_list=list()
    for ligands in sorted(set(pdb_ligand_dict.values())):
        ligands=ligands.split(',')
        for items in ligands:
            ligand_list.append(items)
    ligand_list=sorted(set(ligand_list))
    #ligand_list=sorted(set(pdb_ligand_dict.values()))
    spatial_list=('DFGin','DFGinter','DFGout')
    fhandle_output=open((pwd+'/Representative_structures_ligands.tab'),'w')
    for ligand_name in ligand_list:
        for spatial_label in spatial_list:
            dfgin=0;dfginter=0;dfgout=0
            if spatial_label=='DFGin':
                for dihedral_label in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing_ligand(loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,ligand_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfgin=dfgin+1
                if dfgin==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                    
            if spatial_label=='DFGinter':
                for dihedral_label in ('BABtrans','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing_ligand(loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,ligand_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfginter=dfginter+1
                if dfginter==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                    
            if spatial_label=='DFGout':
                for dihedral_label in ('BBAminus','NA'):
                    (min_missing_pdb,min_missing_residues,min_reso)=least_atoms_missing_ligand(loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict,ligand_name,spatial_label,dihedral_label)
                    if min_missing_pdb!='NA' and dihedral_label!='NA':
                        fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
                        dfgout=dfgout+1
                if dfgout==0 and min_missing_pdb!='NA':      #If their is no cluster assigned to any pdb then print the best noise point
                    fhandle_output.write(f'{ligand_name} {spatial_label} {dihedral_label} {min_missing_pdb} {min_reso} {min_missing_residues}\n')
    fhandle_output.close()
    
    
    