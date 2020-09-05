#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 13:57:40 2020

@author: vivekmodi
"""
import pandas as pd

def read_edia(pwd,df):
    print('Reading EDIA...')
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4]
        pdb_lower=pdb.lower()
        chain=df.at[i,'PDBid'][4]
        uniprotacc=df.at[i,'UniAcc']
        uni_phenum=int(df.at[i,'DFGnum'])
        uni_xnum=uni_phenum-2
        uni_aspnum=uni_phenum-1
        uni_glynum=uni_phenum+1
        
        #Assign default values
        #df.at[i,'X_N_Edia']=999;
        df.at[i,'X_O_Edia']=999;#df.at[i,'Asp_N_Edia']=999;
        df.at[i,'Asp_O_Edia']=999;#df.at[i,'Phe_N_Edia']=999;
        df.at[i,'Phe_O_Edia']=999;#df.at[i,'Gly_N_Edia']=999;
        df.at[i,'Gly_O_Edia']=999
        
        try:               #Try if Edia exists otherwise continue
            edia_df=pd.read_csv(f'{pwd}/mtz_files/Edia_out/{pdb_lower}atomscores.csv',sep=',')
        except:
            continue
        
        sifts_df=pd.read_csv(f'{pwd}/kinasesifts/{pdb_lower}.csv.gz',sep=',')
        
        try:
            pdb_x=int(sifts_df[(sifts_df['Uniprot_resnum'].astype(int)==uni_xnum) & (sifts_df['Uniprot_accession']==uniprotacc) & (sifts_df['PDB_Chain']==chain)].PDB_resnum)
        except:
            pdb_x=999
        try:
            pdb_asp=int(sifts_df[(sifts_df['Uniprot_resnum'].astype(int)==uni_aspnum) & (sifts_df['Uniprot_accession']==uniprotacc) & (sifts_df['PDB_Chain']==chain)].PDB_resnum)
        except:
            pdb_asp=999
        try:
            pdb_phe=int(sifts_df[(sifts_df['Uniprot_resnum'].astype(int)==uni_phenum) & (sifts_df['Uniprot_accession']==uniprotacc) & (sifts_df['PDB_Chain']==chain)].PDB_resnum)
        except:
            pdb_phe=999
        try:
            pdb_gly=int(sifts_df[((sifts_df['Uniprot_resnum']).astype(int)==uni_glynum) & (sifts_df['Uniprot_accession']==uniprotacc) & (sifts_df['PDB_Chain']==chain)].PDB_resnum)
        except:          #if any of the above values is null
            pdb_gly=999
       
               
        
        if pdb_x!=999:
            try:
                #df.at[i,'X_N_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_x)) & (edia_df['Atom name']=='N') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'X_CA_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_x)) & (edia_df['Atom name']=='CA') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'X_C_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_x)) & (edia_df['Atom name']=='C') & (edia_df['Chain']==chain)].EDIA)
                df.at[i,'X_O_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_x)) & (edia_df['Atom name']=='O') & (edia_df['Chain']==chain)].EDIA)
            except:
                pass
            
        if pdb_asp!=999:        
            try:
                #df.at[i,'Asp_N_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_x)) & (edia_df['Atom name']=='N') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Asp_CA_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_asp)) & (edia_df['Atom name']=='CA') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Asp_C_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_asp)) & (edia_df['Atom name']=='C') & (edia_df['Chain']==chain)].EDIA)
                df.at[i,'Asp_O_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_asp)) & (edia_df['Atom name']=='O') & (edia_df['Chain']==chain)].EDIA)
            except:
                pass
                
        if pdb_phe!=999:
            try:
                #df.at[i,'Phe_N_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_phe)) & (edia_df['Atom name']=='N') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Phe_CA_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_phe)) & (edia_df['Atom name']=='CA') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Phe_C_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_phe)) & (edia_df['Atom name']=='C') & (edia_df['Chain']==chain)].EDIA)
                df.at[i,'Phe_O_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_phe)) & (edia_df['Atom name']=='O') & (edia_df['Chain']==chain)].EDIA)
            except:
                pass
                
        if pdb_gly!=999:
            try:
                #df.at[i,'Gly_N_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_gly)) & (edia_df['Atom name']=='N') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Gly_CA_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_gly)) & (edia_df['Atom name']=='CA') & (edia_df['Chain']==chain)].EDIA)
                #df.at[i,'Gly_C_Edia']=float(edia_df[(edia_df['Substructure id']==(pdb_gly)) & (edia_df['Atom name']=='C') & (edia_df['Chain']==chain)].EDIA)
                df.at[i,'Gly_O_Edia']=float(edia_df[(edia_df['Substructure id'].astype(int)==(pdb_gly)) & (edia_df['Atom name']=='O') & (edia_df['Chain']==chain)].EDIA)
            except:
                pass
            
    return df