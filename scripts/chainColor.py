#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:42:41 2020

@author: vivekmodi
"""

def chain_color(df):
    for i in df.index:    
        if df.at[i,'Dihedral']=='BLAminus':
            df.at[i,'Color']=str('(179,215,229)')         #Cyan - ##b3d7e5 
        if df.at[i,'Dihedral']=='BLAplus':
            df.at[i,'Color']=str('(254, 214, 154)')         #Yellowish
        if df.at[i,'Dihedral']=='ABAminus':
            df.at[i,'Color']=str('(243,140,184)')         #Magenta
        if df.at[i,'Dihedral']=='BLBminus':
            df.at[i,'Color']=str('(250, 128, 114)')         #Salmon
        if df.at[i,'Dihedral']=='BLBplus':              
            df.at[i,'Color']=str('(199,233,175)')          #LightGreen -  #c7e9af
        if df.at[i,'Dihedral']=='BLBtrans':
            df.at[i,'Color']=str('(253,156,104)')         #Orange #fd9c68
        if df.at[i,'Dihedral']=='BABtrans':
            df.at[i,'Color']=str('(2,71,254)')          #Darkblue #0247FE
        if df.at[i,'Dihedral']=='BBAminus':
            df.at[i,'Color']=str('(154,180,254)')         #LightBlue - #9ab4fe similar to skyblue
        if df.at[i,'Spatial']=='DFGin' and df.at[i,'Dihedral']=='None':
            df.at[i,'Color']=str('(230,230,230)')         #LightGrey - #e6e6e6 
        if df.at[i,'Spatial']=='DFGinter' and df.at[i,'Dihedral']=='None':
            df.at[i,'Color']=str('(230,230,230)')         #LightGrey - #e6e6e6 
        if df.at[i,'Spatial']=='DFGout' and df.at[i,'Dihedral']=='None':
            df.at[i,'Color']=str('(230,230,230)')         #LightGrey - #e6e6e6 
        if df.at[i,'Spatial']=='None' and df.at[i,'Dihedral']=='None':
            df.at[i,'Color']=str('(230,230,230)')         #LightGrey - #e6e6e6 
    
    return df
        
#    set_color mycyan, (179,215,229)
#    set_color myyellow, (254, 214, 154)
#    set_color mymagenta, (243,140,184)
#    set_color mysalmon, (250, 128, 114)
#    set_color mygreen, (199,233,175)
#    set_color myorange, (253,156,104)
#    set_color myblue, (2,71,254)
#    set_color mylightblue, (154,180,254)