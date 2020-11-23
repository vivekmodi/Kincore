#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 12:17:49 2020

@author: vivekmodi
"""

import os
from ftplib import FTP


def transfer_to_pdbe(pwd):
    print('Transferring JSON files to pdbe...')
    ftp = FTP('10.132.8.198')
    ftp.login(user='vivek', passwd='vivek2')
    ftp.cwd('/home/addit/vivek/check-ftp/')     #remote directory
    
    
    os.chdir(f'{pwd}/JSON/')
    for dirname in filter(os.path.isdir, os.listdir(os.getcwd())):      #Create directories in the remote machine
        try:
            ftp.mkd(f'{dirname}')
        except:
            print(f'Directory {dirname} probably already exists')
    
    for dirname in filter(os.path.isdir, os.listdir(os.getcwd())):     #iterate over local directories
        os.chdir(f'{pwd}/JSON/{dirname}')
        ftp.cwd(f'/home/addit/vivek/check-ftp/{dirname}')              #change directory in the remote machine
        
        for filename in os.listdir():           #iterate over files in individual directories like 'ga','fu'
            myfile=open(f'{filename}', 'rb')    #open the file to upload
            ftp.storlines('STOR '+filename,myfile)     #filename here refers to the name in remote machine
        
        os.chdir(f'{pwd}/JSON/')
    
    
    ftp.quit()
    os.chdir(f'{pwd}')
    