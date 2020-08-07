#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:36:06 2020

@author: vivekmodi
"""
#Usage from commandline python create-blastdb.py pdbaa path
#Usage from another script create_blastdb(pdbaa, path)

import sys, subprocess

def create_blastdb(filename, path):
    print("Creating blast database...")
    cmd=('rm '+path+'/'+filename+'.db.*')
    print(cmd)
    subprocess.call(cmd,shell=True)
    cmd=('makeblastdb -in '+path+'/'+filename+' -out '+path+'/'+filename+'.db -dbtype prot -title '+path+'/'+filename+'-database')
    subprocess.call(cmd,shell=True)
    
if __name__ == '__main__':
    filename=sys.argv[1]
    path=sys.argv[2]
    create_blastdb(filename, path)