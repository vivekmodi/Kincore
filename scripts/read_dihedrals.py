#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:31:55 2020

@author: vivekmodi
"""
import numpy as np

def read_dihedrals(pwd,df):
    kinasechains_dihedrals=f'{pwd}/kinasechains_dihedrals'
    print('Reading dihedrals...')
    # pdb_xphi_dict=dict();pdb_xpsi_dict=dict();pdb_aspphi_dict=dict();pdb_asppsi_dict=dict();pdb_phephi_dict=dict();
    # pdb_aspchi1_dict=dict();pdb_aspchi2_dict=dict();pdb_phepsi_dict=dict();pdb_phechi1_dict=dict();pdb_phechi2_dict=dict();
    # pdb_glyphi_dict=dict();pdb_glypsi_dict=dict();
    
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        #domain_name=df.at[i,'Domain']
        fhandle_dihedral=open((kinasechains_dihedrals+'/'+pdbs+'.dih'),"r")
        # pdb_xphi_dict[pdbs]=999;pdb_xpsi_dict[pdbs]=999;pdb_aspphi_dict[pdbs]=999;pdb_asppsi_dict[pdbs]=999;pdb_aspchi1_dict[pdbs]=999;
        # pdb_aspchi2_dict[pdbs]=999;pdb_phephi_dict[pdbs]=999;pdb_phepsi_dict[pdbs]=999;pdb_phechi1_dict[pdbs]=999;
        # pdb_phechi2_dict[pdbs]=999;pdb_glyphi_dict[pdbs]=999;pdb_glypsi_dict[pdbs]=999;
        pdb_aspchi1=999;pdb_aspchi2=999;pdb_phechi1=999;pdb_phechi2=999
        df.at[i,'X_Phi']=999;df.at[i,'X_Psi']=999;df.at[i,'Asp_Phi']=999;df.at[i,'Asp_Psi']=999;df.at[i,'Asp_Chi1']=999;df.at[i,'Asp_Chi2']=999;
        df.at[i,'Phe_Phi']=999;df.at[i,'Phe_Psi']=999;df.at[i,'Phe_Chi1']=999;df.at[i,'Phe_Chi2']=999;df.at[i,'Gly_Phi']=999;df.at[i,'Gly_Psi']=999;
        
        #if 'HUMAN' not in df.at[i,'UniprotID']:     #skip non-human for now
        #    continue
        for lines in fhandle_dihedral:
            lines=lines.strip("\n");lines=lines.split(" ")
            #print(pdbs)
            if int(lines[3])==int(int(df.at[i,'DFGnum'])-2):     #Match X-DFG
                df.at[i,'X_Phi']=float(lines[5]);df.at[i,'X_Psi']=float(lines[6])
                
            if int(lines[3])==int(int(df.at[i,'DFGnum'])-1):     #Match DFG-Asp
                df.at[i,'Asp_Phi']=float(lines[5]);df.at[i,'Asp_Psi']=float(lines[6]);
                pdb_aspchi1=float(lines[8]);pdb_aspchi2=float(lines[9])
                df.at[i,'Asp_Chi1']=pdb_aspchi1;df.at[i,'Asp_Chi2']=pdb_aspchi2
                
                if pdb_aspchi2<(-90) and pdb_aspchi2>=(-180):
                    df.at[i,'Asp_Chi2']=np.round((pdb_aspchi2+180),2)
                if pdb_aspchi2>90 and pdb_aspchi2<=180:
                    df.at[i,'Asp_Chi2']=np.round((pdb_aspchi2-180),2)
                if pdb_aspchi1<0:
                    df.at[i,'Asp_Chi1']=np.round((pdb_aspchi1+360),2)

            if int(lines[3])==int(df.at[i,'DFGnum']):           #Match DFG-Phe
                df.at[i,'Phe_Phi']=float(lines[5]);df.at[i,'Phe_Psi']=float(lines[6]);
                pdb_phechi1=float(lines[8]);pdb_phechi2=float(lines[9])
                df.at[i,'Phe_Chi1']=pdb_phechi1;df.at[i,'Phe_Chi2']=pdb_phechi2
            
                if pdb_phechi1<0:
                    df.at[i,'Phe_Chi1']=np.round((pdb_phechi1+360),2)
                if pdb_phechi2<0:
                    df.at[i,'Phe_Chi2']=np.round((pdb_phechi2+180),2)
        
            if int(lines[3])==int(int(df.at[i,'DFGnum'])+1):           #Match DFG-Gly
                df.at[i,'Gly_Phi']=float(lines[5]);df.at[i,'Gly_Psi']=float(lines[6])
                break

        
        fhandle_dihedral.close()
    
    return (df)