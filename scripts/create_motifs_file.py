#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 14:21:58 2020

@author: vivekmodi
"""
import pandas as pd

#pwd='/home/vivekmodi/Applications/Flask/Kinases'
def create_motifs_file(pwd):
    print('Creating motifs file...')
    df_align=pd.read_csv(f'{pwd}/All-organisms-alignment-residue-corresspondence.tab',sep=' ')  #make sure this file has a header Org Uniprot Domain ResType UniNum AlignNum
    fhandle_notfound=open(f'{pwd}/Not_found_in-alignment.txt','w')
    fhandle_motifs=open(f'{pwd}/motifs_all.tsv','w')
    fhandle_motifs.write('Gene Domain Group Uniprot Uni_Acc Domain_start Domain_end ALK_res ALK_num RRE_res RRE_num HRD_res HRD_num DFGPhe_res DFGPhe_num DFGAsp_res DFGAsp_num APE_res APE_num DLW_res DLW_num GTK_res GTK_num Hinge1_res Hinge1_num\n')

    count=dict()
    lys=100;glu=146;his_asp=1013;asp=1338;phe=1339;ape=1916;
    gtk=425;hinge1=426;dlw_asp=1961

    for i in df_align.index:
        uniprotid=df_align.at[i,'Uniprot']
        domain=df_align.at[i,'Domain']
        count[uniprotid,domain]=0

    for i in df_align.index:
        uniprotid=df_align.at[i,'Uniprot']
        domain=df_align.at[i,'Domain']
        uniprotacc=df_align.at[i,'UniprotAcc']
        group=df_align.at[i,'Group']
        domain_begin=df_align.at[i,'DomainBegin']
        domain_end=df_align.at[i,'DomainEnd']
        
        alk_num=rre_num=hrd_num=dfg_num=dfgAsp_num=ape_num=hinge1_num=gtk_num=9999
        alk_res=rre_res=hrd_res=dfg_res=dfgAsp_res=ape_res=dlw_res=hinge1_res=gtk_res='X'
        
        if '_1' in domain:
            gene=domain[0:-2]
        elif '_2' in domain:
            gene=domain[0:-2]
        else:
            gene=domain
            
        if count[uniprotid,domain]==0:
            try:
                alk_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==lys)].UniNum)
                alk_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==lys)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} ALKnum\n')
            try:                
                rre_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==glu)].UniNum)
                rre_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==glu)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} RREnum\n')
            try:                
                hrd_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==his_asp)].UniNum)
                hrd_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==his_asp)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} HRDnum\n')
            try:
                dfg_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==phe)].UniNum)
                dfg_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==phe)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} DFGnum\n')
            try:
                dfgAsp_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==asp)].UniNum)
                dfgAsp_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==asp)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} DFGAspnum\n')
            try:
                ape_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==ape)].UniNum)
                ape_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==ape)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} APEnum\n')
            try:
                dlw_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==dlw_asp)].UniNum)
                dlw_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==dlw_asp)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} DLWnum\n')
            try:
                hinge1_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==hinge1)].UniNum)
                hinge1_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==hinge1)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} Hinge1num\n')
            try:
                gtk_num=int(df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==gtk)].UniNum)
                gtk_res=df_align[(df_align.Uniprot==uniprotid) & (df_align.Domain==domain) & (df_align.AlignNum.astype(int)==gtk)].ResType.to_string(index=False).strip()
            except:
                fhandle_notfound.write(f'Not found in alignment {uniprotid} {domain} GTKnum\n')

            fhandle_motifs.write(f'{gene} {domain} {group} {uniprotid} {uniprotacc} {domain_begin} {domain_end} {alk_res} {alk_num} {rre_res} {rre_num} {hrd_res} {hrd_num} {dfg_res} {dfg_num} {dfgAsp_res} {dfgAsp_num} {ape_res} {ape_num} {dlw_res} {dlw_num} {gtk_res} {gtk_num} {hinge1_res} {hinge1_num}\n')

            count[uniprotid,domain]=+1

            
    fhandle_motifs.close()
    fhandle_notfound.close()
