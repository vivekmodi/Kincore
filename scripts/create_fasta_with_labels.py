#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:58:30 2021

@author: vivekmodi
"""
import sys
import pandas as pd
from Bio import SeqIO

    
def create_fasta_with_labels(pwd,df):
    for i in df.index:
        for record in SeqIO.parse(f"{pwd}/pdbaa_psiblast_dir/pdbaa","fasta"):
            if df.at[i,'PDBid']==record.id.split()[0]:
                df.at[i,'PDBAA_Seq']=str(record.seq)
                break
    
    fhandle_fasta=open(f'{pwd}/static/downloads/fasta_with_labels/PK_labels_PDB.fasta','w')
    for i in df.index:
        pdbs=df.at[i,'PDBid'];uniprotid=df.at[i,'UniprotID'];gene=df.at[i,'Gene'];spatial=df.at[i,'Spatial'];dihedral=df.at[i,'Dihedral']
        seq=df.at[i,'PDBAA_Seq']
        fhandle_fasta.write(f'>{pdbs}\t{gene}\t{uniprotid}\t{spatial}\t{dihedral}\n')
        fhandle_fasta.write(f'{seq}\n')
    fhandle_fasta.close()
    
    return df
    

if __name__ == '__main__':
    filename=sys.argv[1]     #Input .csv file
    pwd='/home/vivekmodi/Applications/Flask/Kinases'
    df=pd.read_csv(filename,sep='\t',header='infer')
    df=create_fasta_with_labels(pwd,df)