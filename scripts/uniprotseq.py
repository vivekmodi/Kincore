#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:41:28 2020

@author: vivekmodi
"""
import os

def uniprotseq(pwd,df):
    print('Reading Uniprot and Trembl sequences...')
     
 #   fhandle_uniprot_read=open(f'{pwd}/Uniprot_sequences.txt','r')    
    fhandle_uniprot_append=open(f'{pwd}/Uniprot_sequences.txt','a')
    uniprot_seq=dict()
    
    for i in df.index:
        uniname=df.at[i,'UniprotID']
        
        found=0
        
    #    fhandle_uniprot_read.seek(0)
    #    for lines in fhandle_uniprot_read:
    #           lines=lines.strip();lines=lines.split('\t')
    #           if lines[0]==uniname:
        if uniname in uniprot_seq.keys():
                   #df.at[i,'UniSeq']=lines[4].strip()
                   df.at[i,'UniSeq']=uniprot_seq[uniname]
                   found=1
                   break
        
        if found==0:
            description=os.popen(f'grep -w {uniname} {pwd}/SwissProtIDGeneProteinMapping.csv').read()               
            if description:
                description=description.strip()
                description=description.split('\t')
                df.at[i,'UniSeq']=str(description[5])
                uniprot_seq[uniname]=str(description[5])
                fhandle_uniprot_append.write(f"{uniname}\t{df.at[i,'UniAcc']}\t{df.at[i,'Gene']}\t{df.at[i,'Domain']}\t{df.at[i,'UniSeq']}\n")
                
            else:
                description=os.popen(f'grep -w {uniname} /mnt/Data/Databases/Uniprot/TremblIDGeneProteinMapping.csv').read()               
                if description:
                   description=description.strip()
                   description=description.split('\t')
                   df.at[i,'UniSeq']=description[5].strip()
                   uniprot_seq[uniname]=str(description[5])
                   fhandle_uniprot_append.write(f"{uniname}\t{df.at[i,'UniAcc']}\t{df.at[i,'Gene']}\t{df.at[i,'Domain']}\t{df.at[i,'UniSeq']}\n")
                
    
    
    fhandle_uniprot_append.close()
    #fhandle_uniprot_read.close()
    
    return df