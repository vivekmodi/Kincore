#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 12:51:03 2020

@author: vivekmodi
"""

import pandas as pd

def gene_dict(pwd,df):
    print('Creating gene dictionary...')
    df_motifs=pd.read_csv(f'{pwd}/motifs_all.tsv',sep=' ')
    fhandle_newuniprots=open(f'{pwd}/New_uniprots.txt','w')



    for i in df.index:
        #pdb=df.at[i,'PDBid']
        uniprot=df.at[i,'UniprotID']
        res1=df.at[i,'StrBegin']
        res2=df.at[i,'StrEnd']
        df.at[i,'UniSeq']=''             #initialize this here otherwise in uniprotseq function it gives an error that key does not exist

        if df.at[i,'PDBid']=='4OLIA': #Both domains present, use only the first one
            res2=875

        if uniprot in ('KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN','OBSCN_HUMAN','SPEG_HUMAN','E2AK4_HUMAN',\
                        'JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN','TYK2_HUMAN','KS6A3_MOUSE','JAK2_MOUSE','TYK2_MOUSE'):


                df.at[i,'Gene']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Gene.to_string(index=False).strip()
                df.at[i,'Domain']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Domain.to_string(index=False).strip()
                df.at[i,'Group']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Group.to_string(index=False).strip()
                df.at[i,'UniAcc']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Uni_Acc.to_string(index=False).strip()
                df.at[i,'DomainBegin']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Domain_start)
                df.at[i,'DomainEnd']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Domain_end)
                df.at[i,'ALKres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].ALK_res.to_string(index=False).strip()
                df.at[i,'ALKnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].ALK_num)
                df.at[i,'RREres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].RRE_res.to_string(index=False).strip()
                df.at[i,'RREnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].RRE_num)
                df.at[i,'HRDres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].HRD_res.to_string(index=False).strip()
                df.at[i,'HRDnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].HRD_num)
                df.at[i,'DFGres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].DFGPhe_res.to_string(index=False).strip()
                df.at[i,'DFGnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].DFGPhe_num)
                df.at[i,'DFG_Aspres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].DFGAsp_res.to_string(index=False).strip()
                df.at[i,'APEres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].APE_res.to_string(index=False).strip()
                df.at[i,'APEnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].APE_num)
                df.at[i,'GTKres']=df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].GTK_res.to_string(index=False).strip()
                df.at[i,'GTKnum']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].GTK_num)
                df.at[i,'Hinge1']=int(df_motifs[(df_motifs.Uniprot==uniprot) & (df_motifs.Domain_start<=(res1+res2)/2) & (df_motifs.Domain_end>=(res1+res2)/2)].Hinge1_num)



        else:
            try:
                #print(uniprot,res1,res2)
                df.at[i,'Gene']=df_motifs[df_motifs.Uniprot==uniprot].Gene.to_string(index=False).strip()
                df.at[i,'Domain']=df_motifs[df_motifs.Uniprot==uniprot].Domain.to_string(index=False).strip()
                df.at[i,'Group']=df_motifs[df_motifs.Uniprot==uniprot].Group.to_string(index=False).strip()
                df.at[i,'UniAcc']=df_motifs[df_motifs.Uniprot==uniprot].Uni_Acc.to_string(index=False).strip()
                df.at[i,'DomainBegin']=int(df_motifs[df_motifs.Uniprot==uniprot].Domain_start)
                df.at[i,'DomainEnd']=int(df_motifs[df_motifs.Uniprot==uniprot].Domain_end)
                df.at[i,'ALKres']=df_motifs[df_motifs.Uniprot==uniprot].ALK_res.to_string(index=False).strip()
                df.at[i,'ALKnum']=int(df_motifs[df_motifs.Uniprot==uniprot].ALK_num)
                df.at[i,'RREres']=df_motifs[df_motifs.Uniprot==uniprot].RRE_res.to_string(index=False).strip()
                df.at[i,'RREnum']=int(df_motifs[df_motifs.Uniprot==uniprot].RRE_num)
                df.at[i,'HRDres']=df_motifs[df_motifs.Uniprot==uniprot].HRD_res.to_string(index=False).strip()
                df.at[i,'HRDnum']=int(df_motifs[df_motifs.Uniprot==uniprot].HRD_num)
                df.at[i,'DFGres']=df_motifs[df_motifs.Uniprot==uniprot].DFGPhe_res.to_string(index=False).strip()
                df.at[i,'DFGnum']=int(df_motifs[df_motifs.Uniprot==uniprot].DFGPhe_num)
                df.at[i,'DFG_Aspres']=df_motifs[df_motifs.Uniprot==uniprot].DFGAsp_res.to_string(index=False).strip()
                df.at[i,'APEres']=df_motifs[df_motifs.Uniprot==uniprot].APE_res.to_string(index=False).strip()
                df.at[i,'APEnum']=int(df_motifs[df_motifs.Uniprot==uniprot].APE_num)
                df.at[i,'GTKres']=df_motifs[(df_motifs.Uniprot==uniprot)].GTK_res.to_string(index=False).strip()
                df.at[i,'GTKnum']=int(df_motifs[(df_motifs.Uniprot==uniprot)].GTK_num)
                df.at[i,'Hinge1']=int(df_motifs[(df_motifs.Uniprot==uniprot)].Hinge1_num)
                #df.at[i,'LRLres']='X';df.at[i,'LRLnum']=int(0);df.at[i,'DLWres']='X';df.at[i,'DLWnum']=int(0);df.at[i,'ISRres']='X';df.at[i,'ISRnum']=int(0);  these motifs are not included for now
            except:
                fhandle_newuniprots.write(f"Missing {uniprot} {df.at[i,'PDBid']}\n")
                df.at[i,'Gene']='XXX'
                df.at[i,'Domain']='XXX'
                df.at[i,'Group']='XXX'
                df.at[i,'UniAcc']='XXX'
                df.at[i,'DomainBegin']=int(999)
                df.at[i,'DomainEnd']=int(999)
                df.at[i,'ALKres']='XXX'
                df.at[i,'ALKnum']=int(999)
                df.at[i,'RREres']='XXX'
                df.at[i,'RREnum']=int(999)
                df.at[i,'HRDres']='XXX'
                df.at[i,'HRDnum']=int(999)
                df.at[i,'DFGres']='XXX'
                df.at[i,'DFGnum']=int(999)
                df.at[i,'DFG_Aspres']='XXX'
                df.at[i,'APEres']='XXX'
                df.at[i,'APEnum']=int(999)
                df.at[i,'GTKres']='XXX'
                df.at[i,'GTKnum']=int(999)
                df.at[i,'Hinge1']=int(999)
    fhandle_newuniprots.close()
    return(df)
