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
    
    for i in df.index:
      
        df.at[i,'Spatial']='None'
        dis_phe_glu=float(df.at[i,'Phe_Glu4']);dis_phe_lys=float(df.at[i,'Phe_Lys'])
        if dis_phe_glu<=11 and dis_phe_lys>=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGin'
        elif dis_phe_glu>11 and dis_phe_lys<=14 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGout'
        elif dis_phe_glu<=11 and dis_phe_lys<=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial']='DFGinter'

    return df

#def cosine_dis_with_chi1(df,i,spatial,cutoff):
#    x_dfg_phi=float(df.at[i,'X_Phi']);x_dfg_psi=float(df.at[i,'X_Psi']);dfg_asp_phi=float(df.at[i,'Asp_Phi']);
#    dfg_asp_psi=float(df.at[i,'Asp_Psi']);dfg_phe_phi=float(df.at[i,'Phe_Phi'])
#    dfg_phe_psi=float(df.at[i,'Phe_Psi']);dfg_phe_chi1=float(df.at[i,'Phe_Chi1'])
#    min_spatial=999
#    
#    for clusters in spatial:    
#        cosine_dis=(2/7)*((1-math.cos(math.radians(x_dfg_phi-float(spatial[clusters][0]))))+(1-math.cos(math.radians(x_dfg_psi-float(spatial[clusters][1]))))+\
#                (1-math.cos(math.radians(dfg_asp_phi-float(spatial[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(spatial[clusters][3]))))+\
#                (1-math.cos(math.radians(dfg_phe_phi-float(spatial[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(spatial[clusters][5]))))+\
#                (1-math.cos(math.radians(dfg_phe_chi1-float(spatial[clusters][6])))))
#        
#        if cosine_dis<=min_spatial:
#            df.at[i,'Dihedral_dis']=np.round(cosine_dis,2)
#            min_spatial=cosine_dis
#            
#            if cosine_dis<=cutoff:
#                df.at[i,'Dihedral']=clusters
#        
#    return df

def cosine_dis_without_chi1(df,i,spatial,cutoff):
    x_dfg_phi=float(df.at[i,'X_Phi']);x_dfg_psi=float(df.at[i,'X_Psi']);dfg_asp_phi=float(df.at[i,'Asp_Phi']);
    dfg_asp_psi=float(df.at[i,'Asp_Psi']);dfg_phe_phi=float(df.at[i,'Phe_Phi'])
    dfg_phe_psi=float(df.at[i,'Phe_Psi'])
    min_spatial=999
    
    for clusters in spatial:    
        cosine_dis=(2/6)*((1-math.cos(math.radians(x_dfg_phi-float(spatial[clusters][0]))))+(1-math.cos(math.radians(x_dfg_psi-float(spatial[clusters][1]))))+\
                (1-math.cos(math.radians(dfg_asp_phi-float(spatial[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(spatial[clusters][3]))))+\
                (1-math.cos(math.radians(dfg_phe_phi-float(spatial[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(spatial[clusters][5])))))
        
        if cosine_dis<=min_spatial:
            df.at[i,'Dihedral_dis_NoChi1']=np.round(cosine_dis,2)
            min_spatial=cosine_dis
            
            if cosine_dis<=float(cutoff):
                df.at[i,'Dihedral']=clusters      #Only Dihedral column name is used for final labeling without chi1
    return df

    
def dihedral_labels(df,cutoff):
    print('Assigning dihedral labels...')
    
    dfginter={'BABtrans':(-80.20,128.22,-117.47,23.76,-85.16,133.21,181.42)}
    dfgout={'BBAminus':(-138.56,-176.12,-144.35,103.66,-82.59,-9.03,290.59)}
    
    dfgin_minus={'BLAminus':(-128.64,178.67,61.15,81.21,-96.89,20.53,289.12),\
       'ABAminus':(-111.82,-7.64,-141.55,148.01,-127.79,23.32,296.17),\
       'BLBminus':(-134.79,175.48,60.44,65.35,-79.44,145.34,287.56)}
    
    dfgin_plus={'BLAplus':(-119.24,167.71,58.94,34.08,-89.42,-8.54,55.63),\
       'BLBplus':(-125.28,172.53,59.98,32.92,-85.51,145.28,49.01)}
    
    dfgin_trans={'BLBtrans':(-106.16,157.24,69.37,21.33,-61.73,134.56,215.23)}
    
    
    for i in df.index:
        x_dfg_phi=float(df.at[i,'X_Phi']);x_dfg_psi=float(df.at[i,'X_Psi']);dfg_asp_phi=float(df.at[i,'Asp_Phi']);
        dfg_asp_psi=float(df.at[i,'Asp_Psi']);dfg_phe_phi=float(df.at[i,'Phe_Phi'])
        dfg_phe_psi=float(df.at[i,'Phe_Psi']);dfg_phe_chi1=float(df.at[i,'Phe_Chi1'])
        
        df.at[i,'Dihedral']='None'    #Default label is None
        df.at[i,'Dihedral_dis_NoChi1']=999
        
#        if df.at[i,'Spatial']=='Unassigned':
#            df.at[i,'Dihedral']='Unassigned'
#            df.at[i,'Dihedral_dis_NoChi1']=999
        
        if x_dfg_phi==999 or x_dfg_psi==999 or dfg_asp_phi==999 or dfg_asp_psi==999 or dfg_phe_phi==999 or  \
        dfg_phe_psi==999 or dfg_phe_chi1==999:
            df.at[i,'Dihedral']='None'
            df.at[i,'Dihedral_dis_NoChi1']=999
            continue
        
        if df.at[i,'Spatial']=='DFGin':     
            if float(df.at[i,'Phe_Chi1'])>240 and float(df.at[i,'Phe_Chi1'])<=360:
                df=cosine_dis_without_chi1(df,i,dfgin_minus,cutoff)
            if float(df.at[i,'Phe_Chi1'])>0 and float(df.at[i,'Phe_Chi1'])<=120:
                df=cosine_dis_without_chi1(df,i,dfgin_plus,cutoff)
            if float(df.at[i,'Phe_Chi1'])>120 and float(df.at[i,'Phe_Chi1'])<=240:
                df=cosine_dis_without_chi1(df,i,dfgin_trans,cutoff)
            
        if df.at[i,'Spatial']=='DFGinter': 
            if float(df.at[i,'Phe_Chi1'])>120 and float(df.at[i,'Phe_Chi1'])<=240:
                df=cosine_dis_without_chi1(df,i,dfginter,cutoff)
        
        if df.at[i,'Spatial']=='DFGout':        
            if float(df.at[i,'Phe_Chi1'])>240 and float(df.at[i,'Phe_Chi1'])<=360:
                df=cosine_dis_without_chi1(df,i,dfgout,cutoff)

    return df