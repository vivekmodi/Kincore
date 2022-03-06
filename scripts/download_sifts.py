#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:01:50 2020

@author: vivekmodi
"""


import os, subprocess

def download_sifts(downld_dir, df):
    print('Downloading files from Sifts...')
    for i in df.index:

        pdbs=df.at[i,'PDBid']
        filename=pdbs[0:4]

        if not os.path.isfile(downld_dir+"/"+str(filename).lower()+".xml.gz"):
            print("Downloading Sifts "+filename+"...")
            cmd=("wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/"+str(filename).lower()+".xml.gz -P "+downld_dir)
            subprocess.call(cmd,shell=True)
