#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 14:15:45 2020

@author: vivekmodi
"""
import gzip, os
from Bio import SeqIO

def genes_from_trembl(pwd,df):
    print('Identifying genes from Trembl...')
    #fhandleUni=gzip.open('/mnt/Data/Databases/Uniprot/uniprot_trembl.fasta.gz','rt')
    fhandle_output_read=open(f'{pwd}/Trembl_genenames.txt','r')    
    fhandle_output_append=open(f'{pwd}/Trembl_genenames.txt','a')
    
    for i in df.index:
        if df.at[i,'Gene']=='X':
           uniname=df.at[i,'UniprotID']
           found=0
           
           fhandle_output_read.seek(0)
           for lines in fhandle_output_read:
               lines=lines.strip();lines=lines.split('\t')
               if lines[0]==uniname:
                   
                   df.at[i,'UniAcc']=lines[1].strip()
                   df.at[i,'Gene']=lines[2].strip().upper()
                   df.at[i,'Domain']=lines[3].strip().upper()
                   df.at[i,'UniSeq']=lines[4].strip()
                   found=1
                   break
           
           
           if found==0:
               description=os.popen(f'grep -w {uniname} /mnt/Data/Databases/Uniprot/TremblIDGeneProteinMapping.csv').read()               
               if description:
                   description=description.split('\t')
                   
                   df.at[i,'UniAcc']=description[0].strip()
                   df.at[i,'Gene']=description[2].strip().upper()
                   df.at[i,'Domain']=description[2].strip().upper()
                   df.at[i,'UniSeq']=description[5].strip()
                   fhandle_output_append.write(f'{uniname}\t{description[0]}\t{description[2]}\t{description[2]}\t{description[5]}\n')
    
    fhandle_output_append.close()
    fhandle_output_read.close()
    return df
                
                