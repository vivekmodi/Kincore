#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:04:12 2020

@author: vivekmodi
"""
from datetime import datetime
from collections import defaultdict
import subprocess

def statistics(pwd,pdb_domain_dict,pdb_gene_dict,pdb_group_dict,pdb_uniprot_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Calculating statistics...')
    today=str(datetime.now())[0:10].strip()
    fhandle_statistics=open((pwd+'/Statistics-'+today+'.tab'),'w')
    fhandle_domain=open((pwd+'/Domain-count-'+today+'.tab'),'w')
    fhandle_group=open((pwd+'/Group-count-'+today+'.tab'),'w')
    fhandle_spatial=open((pwd+'/Label-count-spatial-'+today+'.tab'),'w')
    fhandle_dihedral=open((pwd+'/Label-count-dihedral-'+today+'.tab'),'w')
    pdbname=list();count_spatial=dict();count_dihedral=dict();count_domain=dict();count_group=dict();group_total=dict()
    count_group_spatial=defaultdict(dict)       #Initialize two dimensional dictionary
    count_group_dihedral=defaultdict(dict)
    chain_num=0;pdb_num=0;uni_num=0;gene_num=0;
    spatial_list=('DFGin','DFGinter','DFGout','NA')
    dihedral_list=('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','BABtrans','BBAminus','NA')
    group_list=('AGC','CAMK','CK1','CMGC','OTHER','NEK','RGC','STE','TKL','TYR')
   
    for pdbs in pdb_domain_dict:            
        #count_spatial[pdb_spatial_dict[pdbs]]=0
        #count_dihedral[pdb_dihedral_dict[pdbs]]=0
        count_domain[pdb_domain_dict[pdbs]]=0
        #count_group[pdb_group_dict[pdbs]]=0
    for spatials in spatial_list:
        count_spatial[spatials]=0
    for dihedrals in dihedral_list:
        count_dihedral[dihedrals]=0
    for groups in group_list:
        count_group[groups]=0
        group_total[groups]=0
        
             
    for spatials in spatial_list:
        for groups in group_list:
            count_group_spatial[spatials][groups]=0
    for dihedrals in dihedral_list:
        for groups in group_list:
            count_group_dihedral[dihedrals][groups]=0       #datastructure {DFGin:{AGC:0}}
            #print(f'{groups} {dihedrals} {count_group_dihedral[groups][dihedrals]}\n')
    
    for pdbs in pdb_domain_dict:
        pdbname.append(pdbs[0:4])
        count_spatial[pdb_spatial_dict[pdbs]]=count_spatial[pdb_spatial_dict[pdbs]]+1
        count_dihedral[pdb_dihedral_dict[pdbs]]=count_dihedral[pdb_dihedral_dict[pdbs]]+1
        count_domain[pdb_domain_dict[pdbs]]=count_domain[pdb_domain_dict[pdbs]]+1
        count_group[pdb_group_dict[pdbs]]=count_group[pdb_group_dict[pdbs]]+1
        count_group_spatial[pdb_spatial_dict[pdbs]][pdb_group_dict[pdbs]]=count_group_spatial[pdb_spatial_dict[pdbs]][pdb_group_dict[pdbs]]+1
        count_group_dihedral[pdb_dihedral_dict[pdbs]][pdb_group_dict[pdbs]]=count_group_dihedral[pdb_dihedral_dict[pdbs]][pdb_group_dict[pdbs]]+1
        
    chain_num=len(pdb_domain_dict)
    pdb_num=len(set(pdbname))
    uni_num=len(set(pdb_uniprot_dict.values()))
    gene_num=len(set(pdb_gene_dict.values()))
    sum_chain_spatial=0;sum_chain_dihedral=0
    fhandle_statistics.write(f"PDBs {pdb_num}\nChains {chain_num}\nGene {gene_num}\nUniprot {uni_num}\n")
    
    fhandle_statistics.write(f"Spatial DFGin DFGinter DFGout NA Total\nChains ")
    for spatials in spatial_list:
    #for keys in count_spatial:
        fhandle_statistics.write(f"{count_spatial[spatials]} ")
        sum_chain_spatial=sum_chain_spatial+count_spatial[spatials]
    fhandle_statistics.write(f"{sum_chain_spatial}\n")
    fhandle_statistics.write(f"Dihedral DFGin-BLAminus DFGin-BLAplus DFGin-ABAminus DFGin-BLBminus DFGin-BLBplus DFGin-BLBtrans DFGinter-BABtrans DFGout-BBAminus NA Total\nChains ")
    for dihedrals in dihedral_list:
    #for keys in count_dihedral:
        fhandle_statistics.write(f"{count_dihedral[dihedrals]} ")
        sum_chain_dihedral=sum_chain_dihedral+count_dihedral[dihedrals]
    fhandle_statistics.write(f"{sum_chain_dihedral}\n")
    fhandle_statistics.close()
    
    for keys in count_domain:
        fhandle_domain.write(f"{keys} {count_domain[keys]}\n")
    fhandle_domain.close()
    
    for groups in group_list:
    #for keys in count_group:
        fhandle_group.write(f"{groups} {count_group[groups]}\n")
    fhandle_group.close()
    
    fhandle_spatial.write(f'Spatial_label AGC CAMK CK1 CMGC OTHER NEK RGC STE TKL TYR Total\n')
    for spatial_labels in spatial_list:
        fhandle_spatial.write(f'{spatial_labels} ')
        for group_names in group_list:
            fhandle_spatial.write(f'{count_group_spatial[spatial_labels][group_names]} ')
        fhandle_spatial.write(f'{count_spatial[spatial_labels]}\n')
    fhandle_spatial.write(f"Total {count_group['AGC']} {count_group['CAMK']} {count_group['CK1']} {count_group['CMGC']} {count_group['OTHER']} {count_group['NEK']} {count_group['RGC']} {count_group['STE']} {count_group['TKL']} {count_group['TYR']} {chain_num}\n")
    fhandle_spatial.close()
    
    fhandle_dihedral.write(f'Dihedral_label AGC CAMK CK1 CMGC OTHER NEK RGC STE TKL TYR Total\n')
    for dihedral_labels in dihedral_list:
        fhandle_dihedral.write(f'{dihedral_labels} ')
        for group_names in group_list:
            fhandle_dihedral.write(f'{count_group_dihedral[dihedral_labels][group_names]} ')
        fhandle_dihedral.write(f'{count_dihedral[dihedral_labels]}\n')
    fhandle_dihedral.write(f"Total {count_group['AGC']} {count_group['CAMK']} {count_group['CK1']} {count_group['CMGC']} {count_group['OTHER']} {count_group['NEK']} {count_group['RGC']} {count_group['STE']} {count_group['TKL']} {count_group['TYR']} {chain_num}\n")
    fhandle_dihedral.close()
    
    cmd=('column -t '+'Label-count-spatial-'+today+'.tab > temp\nmv temp '+'Label-count-spatial-'+today+'.tab')
    subprocess.call(cmd,shell=True)
    cmd=('column -t '+'Label-count-dihedral-'+today+'.tab > temp\nmv temp '+'Label-count-dihedral-'+today+'.tab')
    subprocess.call(cmd,shell=True)