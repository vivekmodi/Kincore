#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:41:28 2020

@author: vivekmodi
"""
import os
import pandas as pd

def uniprotseq(pwd,df):
    print('Reading Uniprot and Trembl sequences...')

    df_seq=pd.read_csv('Uniprot_sequences.txt',sep='\t')     
    #fhandle_uniprot_read=open(f'{pwd}/Uniprot_sequences.txt','r')    
    #fhandle_uniprot_append=open(f'{pwd}/Uniprot_sequences.txt','a')
    #uniprot_seq=dict()
    
    for i in df.index:
        uniname=df.at[i,'UniprotID']
        found=0
        
        if df_seq[df_seq.Uniprot==uniname].empty:    #If empty then found=0
            found=0            
        else:
            df.at[i,'UniSeq']=df_seq[df_seq.Uniprot==uniname].Sequence.to_string(index=False).strip()
            found=1
            
        if found==0:
            description=os.popen(f'grep -w {uniname} {pwd}/SwissProtIDGeneProteinMapping.csv').read()               
            if description:
                description=description.strip();description=description.split('\t')
                df.at[i,'UniSeq']=str(description[5])
                df_seq=df_seq.append({'Uniprot':uniname,'Sequence':str(description[5])},ignore_index=True)
                found=1
                
        if found==0:
            description=os.popen(f'grep -w {uniname} /mnt/Data/Databases/Uniprot/TremblIDGeneProteinMapping.csv').read()               
            if description:
                description=description.strip()
                description=description.split('\t')
                df.at[i,'UniSeq']=description[5].strip()
                df_seq=df_seq.append({'Uniprot':uniname,'Sequence':str(description[5])},ignore_index=True)
                found=1
        
    #    fhandle_uniprot_read.seek(0)
    #    for lines in fhandle_uniprot_read:
    #           lines=lines.strip();lines=lines.split('\t')
    #           if lines[0]==uniname:
#        if uniname in uniprot_seq.keys():
                   #df.at[i,'UniSeq']=lines[4].strip()
#                   df.at[i,'UniSeq']=uniprot_seq[uniname]
#                   found=1
#                   break
        
#        if found==0:
#            description=os.popen(f'grep -w {uniname} {pwd}/SwissProtIDGeneProteinMapping.csv').read()               
#            if description:
#                description=description.strip()
#                description=description.split('\t')
#                df.at[i,'UniSeq']=str(description[5])
#                uniprot_seq[uniname]=str(description[5])
#                fhandle_uniprot_append.write(f"{uniname}\t{df.at[i,'UniAcc']}\t{df.at[i,'Gene']}\t{df.at[i,'Domain']}\t{df.at[i,'UniSeq']}\n")
#                
#            else:
#                description=os.popen(f'grep -w {uniname} /mnt/Data/Databases/Uniprot/TremblIDGeneProteinMapping.csv').read()               
#                if description:
#                   description=description.strip()
#                   description=description.split('\t')
#                   df.at[i,'UniSeq']=description[5].strip()
#                   uniprot_seq[uniname]=str(description[5])
#                   fhandle_uniprot_append.write(f"{uniname}\t{df.at[i,'UniAcc']}\t{df.at[i,'Gene']}\t{df.at[i,'Domain']}\t{df.at[i,'UniSeq']}\n")
#                
#    
#    
#    fhandle_uniprot_append.close()
    #fhandle_uniprot_read.close()
    #print(df_seq.head())
    df_seq.to_csv('Uniprot_sequences.txt',sep='\t',index=False)
    return df