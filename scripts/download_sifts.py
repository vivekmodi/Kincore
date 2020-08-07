#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:01:50 2020

@author: vivekmodi
"""

#Usage from commandline python download_sifts.py /download/sifts 1xyz
#Usage from another script download_sifts(downld_dir, filename)

import os, subprocess

def download_sifts(downld_dir, df):
    print('Downloading files from Sifts...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        filename=pdbs[0:4]
   # cmd=("rm ./kinasesifts/*.gz*")            #Before downloading new files delete previous xml and .gz files otherwise it creates new files with extension like .gz1 or .gz2.
   # subprocess.call(cmd,shell=True)
   # cmd=("rm ./kinasesifts/*.xml*")
   # subprocess.call(cmd,shell=True)
        if not os.path.isfile(downld_dir+"/"+str(filename).lower()+".xml.gz"):
            
                print("Downloading Sifts "+filename+"...")
                cmd=("wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/"+str(filename).lower()+".xml.gz -P "+downld_dir)
                subprocess.call(cmd,shell=True)
           
           
#if __name__ == '__main__':
#    downld_dir=sys.argv[1]
#    filename=sys.argv[2]
#    download_sifts(downld_dir, filename)