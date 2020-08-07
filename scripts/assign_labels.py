#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:13:13 2020

@author: vivekmodi
"""

import math
import numpy as np

def spatial_labels(pwd,df):
    print('Assigning spatial labels...')
    #pdb_spatial_dict=dict()
    for i in df.index:
        #pdbs=df.at[i,'PDBid']
        #pdb_spatial_dict[pdbs]='NA'
        df.at[i,'Spatial']='None'
        dis_phe_glu=float(df.at[i,'Phe_Glu4']);dis_phe_lys=float(df.at[i,'Phe_Lys'])
        if dis_phe_glu<=11 and dis_phe_lys>=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGin'
        elif dis_phe_glu>11 and dis_phe_lys<=14 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGout'
        elif dis_phe_glu<=11 and dis_phe_lys<=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGinter'
        
            
        #print(f'{pdbs} {pdb_pheglu4_dict[pdbs]} {pdb_phelys_dict[pdbs]} {pdb_spatial_dict[pdbs]}')
    return df


def dihedral_labels(pwd,df):
    fhandle_cosine_with_chi=open(f'{pwd}/cosine_with_chi','w')
    print('Assigning dihedral labels...')
    fhandle_dfgin=open(f"{pwd}/Cluster-centroids-DFGin-filtered-2.25","r")
    fhandle_dfginter=open(f"{pwd}/Cluster-centroids-DFGinter-filtered-2.25","r")
    fhandle_dfgout=open(f"{pwd}/Cluster-centroids-DFGout-filtered-2.25","r")
    
    #pdb_dihedral_dict=dict()
    
    for i in df.index:
        min_cosine=999
        dist_label=df.at[i,'Spatial']
        #pdb_dihedral_dict[pdbs]='NA'
        #if dist_label!='':
            #pdb_dihedral_dict[pdbs]='NA'
        #else:
        x_dfg_phi=float(df.at[i,'X_Phi']);x_dfg_psi=float(df.at[i,'X_Psi']);dfg_asp_phi=float(df.at[i,'Asp_Phi']);
        dfg_asp_psi=float(df.at[i,'Asp_Psi']);dfg_phe_phi=float(df.at[i,'Phe_Phi'])
        dfg_phe_psi=float(df.at[i,'Phe_Psi']);dfg_phe_chi1=float(df.at[i,'Phe_Chi1'])
    
        #if x_dfg_phi==999 or x_dfg_psi==999 or dfg_asp_phi==999 or dfg_asp_psi==999 or dfg_phe_phi==999 or  dfg_phe_psi==999 or dfg_phe_chi1==999:
        df.at[i,'Dihedral']='None'
         #   continue
    
        if dist_label=='DFGin':
            min_cluster=0;
            for lines in fhandle_dfgin:
                lines=lines.strip("\n");lines=lines.split(" ")
    
                cosine_dis=(2/7)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                            (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))
               
                if cosine_dis<=min_cosine:
                    min_cosine=cosine_dis
                if cosine_dis<0.3:
                    min_cluster=int(lines[1])
                #print(f'{pdbs} {cosine_dis} {dfg_phe_chi1}')
                                        
            #if min_cluster==0:
            #    pdb_dihedral_dict[pdbs]='NA'
            #if min_cluster>0:
            if min_cluster==1:
                df.at[i,'Dihedral']='BLAminus'
            if min_cluster==2:
                df.at[i,'Dihedral']='BLAplus'
            if min_cluster==3:
                df.at[i,'Dihedral']='ABAminus'
            if min_cluster==4:
                df.at[i,'Dihedral']='BLBminus'
            if min_cluster==5:
                df.at[i,'Dihedral']='BLBplus'
            if min_cluster==6:
                df.at[i,'Dihedral']='BLBtrans'
          
#                else:
#                    pdb_dihedral_dict[pdbs]='NA'
            fhandle_cosine_with_chi.write(f'{df.at[i,"PDBid"]} {dist_label} {dfg_phe_chi1} {min_cosine}\n')        
        fhandle_dfgin.seek(0)

        if dist_label=='DFGinter':
            min_cluster=0
            for lines in fhandle_dfginter:
                lines=lines.strip("\n");lines=lines.split(" ")
                cosine_dis=(2/7)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                            (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))
                
                if cosine_dis<0.3:
                    df.at[i,'Dihedral']='BABtrans'
                #else:
                #    df.at[i,'Dihedral']=''
        fhandle_dfginter.seek(0)

        if dist_label=='DFGout':
            min_cluster=0
            for lines in fhandle_dfgout:
                lines=lines.strip("\n");lines=lines.split(" ")
                cosine_dis=(2/7)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                            (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))               
               
                if cosine_dis<0.3:
                    df.at[i,'Dihedral']='BBAminus'
                #else:
                #    df.at[i,'Dihedral']=''
        fhandle_dfgout.seek(0)
         #   print(f'{pdbs} {pdb_dihedral_dict[pdbs]}')   
    fhandle_cosine_with_chi.close()             
    return df
            

def dihedral_labels_no_chi(pwd,df):
    print('Assigning dihedral labels...')
    fhandle_dfgin=open(f"{pwd}/Cluster-centroids-DFGin-filtered-2.25","r")
    fhandle_dfginter=open(f"{pwd}/Cluster-centroids-DFGinter-filtered-2.25","r")
    fhandle_dfgout=open(f"{pwd}/Cluster-centroids-DFGout-filtered-2.25","r")
    fhandle_cosine_no_chi=open(f'{pwd}/cosine_no_chi','w')
    #pdb_dihedral_dict=dict()
    for i in df.index:
        
        dist_label=df.at[i,'Spatial']
        #pdb_dihedral_dict[pdbs]='NA'
        #if dist_label!='':
            #pdb_dihedral_dict[pdbs]='NA'
        #else:
        x_dfg_phi=float(df.at[i,'X_Phi']);x_dfg_psi=float(df.at[i,'X_Psi']);dfg_asp_phi=float(df.at[i,'Asp_Phi']);
        dfg_asp_psi=float(df.at[i,'Asp_Psi']);dfg_phe_phi=float(df.at[i,'Phe_Phi'])
        dfg_phe_psi=float(df.at[i,'Phe_Psi']);dfg_phe_chi1=float(df.at[i,'Phe_Chi1'])
    
        #if x_dfg_phi==999 or x_dfg_psi==999 or dfg_asp_phi==999 or dfg_asp_psi==999 or dfg_phe_phi==999 or  dfg_phe_psi==999 or dfg_phe_chi1==999:
        #    df.at[i,'Dihedral_NoChi1']=''
        #    continue
        df.at[i,'Dihedral_NoChi1']='None'
        if dist_label=='DFGin':
            min_cluster=0;min_cosine=999
            for lines in fhandle_dfgin:
                lines=lines.strip("\n");lines=lines.split(" ")
    
                cosine_dis=(2/6)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7])))))
                
                #if (1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))<=0.14 and (1-math.cos(math.radians(x_dfg_psi-float(lines[3])))) <=0.14 and \
                #            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4])))) <=0.14 and  (1-math.cos(math.radians(dfg_asp_psi-float(lines[5])))) <=0.14 and \
                #            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6])))) <=0.14 and (1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))<=0.14:
                                        
                if dfg_phe_chi1>240 and dfg_phe_chi1<=360:
                    if cosine_dis<0.3:
                        if cosine_dis<=min_cosine:
                            min_cosine=cosine_dis
                            min_cluster=int(lines[1])
                        
                            if min_cluster==1:
                                df.at[i,'Dihedral_NoChi1']='BLAminus'
                            if min_cluster==3:
                                df.at[i,'Dihedral_NoChi1']='ABAminus'
                            if min_cluster==4:
                                df.at[i,'Dihedral_NoChi1']='BLBminus'
                        
                if dfg_phe_chi1>0 and dfg_phe_chi1<=120:
                    if cosine_dis<0.3:
                        if cosine_dis<=min_cosine:
                            min_cosine=cosine_dis
                    
                            min_cluster=int(lines[1])
                            if min_cluster==2:
                                df.at[i,'Dihedral_NoChi1']='BLAplus'
                            if min_cluster==5:
                                df.at[i,'Dihedral_NoChi1']='BLBplus'
                        
                if dfg_phe_chi1>120 and dfg_phe_chi1<=240:
                    if cosine_dis<0.3:
                        if cosine_dis<=min_cosine:
                            min_cosine=cosine_dis
                    
                            min_cluster=int(lines[1])
                            if min_cluster==6:
                                df.at[i,'Dihedral_NoChi1']='BLBtrans'
          
#                else:
#                    pdb_dihedral_dict[pdbs]='NA'
            fhandle_cosine_no_chi.write(f'{df.at[i,"PDBid"]} {dist_label} {dfg_phe_chi1} {min_cluster} {df.at[i,"Dihedral_NoChi1"]} {min_cosine}\n')        
        fhandle_dfgin.seek(0)

        if dist_label=='DFGinter':
            min_cluster=0;min_cosine=999
            for lines in fhandle_dfginter:
                lines=lines.strip("\n");lines=lines.split(" ")
                cosine_dis=(2/7)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                            (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))
                
                if dfg_phe_chi1>120 and dfg_phe_chi1<=240:
                    if cosine_dis<0.3:
                        if cosine_dis<=min_cosine:
                            min_cosine=cosine_dis
                    
                            min_cluster=int(lines[1])
                            if min_cluster==1:
                                df.at[i,'Dihedral_NoChi1']='BABtrans'
            fhandle_cosine_no_chi.write(f'{df.at[i,"PDBid"]} {dist_label} {dfg_phe_chi1} {min_cluster} {df.at[i,"Dihedral_NoChi1"]} {min_cosine}\n')
        fhandle_dfginter.seek(0)

        if dist_label=='DFGout':
            min_cluster=0;min_cosine=999
            for lines in fhandle_dfgout:
                lines=lines.strip("\n");lines=lines.split(" ")
                cosine_dis=(2/6)*((1-math.cos(math.radians(x_dfg_phi-float(lines[2]))))+(1-math.cos(math.radians(x_dfg_psi-float(lines[3]))))+\
                            (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                            (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7])))))
                
                if dfg_phe_chi1>240 and dfg_phe_chi1<=360:
                    if cosine_dis<0.3:
                        if cosine_dis<=min_cosine:
                            min_cosine=cosine_dis
                    
                            min_cluster=int(lines[1])
                            if min_cluster==1:
                                df.at[i,"Dihedral_NoChi1"]='BBAminus'
            fhandle_cosine_no_chi.write(f'{df.at[i,"PDBid"]} {dist_label} {dfg_phe_chi1} {min_cluster} {df.at[i,"Dihedral_NoChi1"]} {min_cosine}\n')
        fhandle_dfgout.seek(0)
         #   print(f'{pdbs} {pdb_dihedral_dict[pdbs]}')                
    fhandle_cosine_no_chi.close()
    return df