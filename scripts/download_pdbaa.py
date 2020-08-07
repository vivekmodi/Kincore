#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:18:55 2020

@author: vivekmodi
"""
#Usage from commandline python download_pdbaa.py /path
#Usage from another script download_pdbaa(path)

import subprocess, sys

def download_pdbaa(path):
    print("Downloading pdbaa from dunbrack.fccc.edu...")
    cmd=("rm "+path+"/pdbaa.gz")
    subprocess.call(cmd, shell=True)
    cmd=("wget -P "+path+" http://dunbrack.fccc.edu/Guoli/culledpdb_hh/pdbaa.gz")
    subprocess.call(cmd,shell=True)
    cmd=("gunzip -f "+path+"/pdbaa.gz")
    subprocess.call(cmd,shell=True)
    
if __name__ == '__main__':
    path=sys.argv[1]
    download_pdbaa(path)