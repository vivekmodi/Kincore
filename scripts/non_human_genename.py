#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:24:11 2020

@author: vivekmodi
"""
import os

def non_human_genename(pwd,df):
    print('Identifying gene names for non human proteins...')
    #fhandle_mapping=open(f'{pwd}/SwissProtIDGeneProteinMapping.csv','r')
    fhandle_output_read=open(f'{pwd}/Non_Human_genenames.txt','r')    
    fhandle_output_append=open(f'{pwd}/Non_Human_genenames.txt','a')
    for i in df.index:
        uniname=df.at[i,'UniprotID']
        if 'HUMAN' in uniname:
            continue
        found=0
        
        fhandle_output_read.seek(0)
        for lines in fhandle_output_read:
            lines=lines.strip();lines=lines.split()
            if lines[0]==uniname:
                df.at[i,'UniAcc']=lines[1]
                df.at[i,'Gene']=lines[2].upper()
                df.at[i,'Domain']=lines[3].upper()
                found=1
                break
        
        if found==0:
                
            description=os.popen(f'grep -w {uniname} {pwd}/SwissProtIDGeneProteinMapping.csv').read()               
            if description:
                description=description.split('\t')
                df.at[i,'UniAcc']=description[0]
                df.at[i,'Gene']=description[2].upper()
                df.at[i,'Domain']=description[2].upper()
                fhandle_output_append.write(f'{uniname}\t{description[0]}\t{description[2]}\t{description[2]}\n')
            
     #   fhandle_mapping.seek(0)
     #   for lines in fhandle_mapping:
     #       lines=lines.strip();lines=lines.split('\t')
     #       if lines[1]==uniname:
     #           df.at[i,'UniAcc']=lines[0]
     #           df.at[i,'Gene']=lines[2]
     #           df.at[i,'Domain']=lines[2]    #Change this in future for proteins with two domains
     #           break
    fhandle_output_append.close()
    fhandle_output_read.close()
    return df