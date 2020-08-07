#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:48:50 2020

@author: vivekmodi
"""


import os, subprocess
from Bio import PDB

def download_cifs(downld_dir,df):
    print('Downloading MMCIFs...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        filename=pdbs[0:4]
        #if not os.path.isfile(downld_dir+"/"+str(filename).lower()+".cif"):   #Download if file does not exist
        if not os.path.isfile(downld_dir+"/"+str(filename).lower()+".cif.gz"):
                #print("Downloading "+filename+"...")
                pdb1=PDB.PDBList()
                pdb1.retrieve_pdb_file(filename,pdir=downld_dir)
    
        if os.path.isfile(downld_dir+"/"+str(filename).lower()+".cif"):     #gzip files
            #print("gzip "+filename)
            cmd=("gzip -f "+downld_dir+"/"+str(filename).lower()+".cif")
            subprocess.call(cmd, shell=True)

#if __name__ == '__main__':
#    downld_dir=sys.argv[1]
#    filename=sys.argv[2]
#    download_cifs(downld_dir,filename)