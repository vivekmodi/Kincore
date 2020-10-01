#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:36:19 2020

@author: vivekmodi
"""
import subprocess, os

def copy_ngl_files(pwd,df):
        print('Copying .pdb files for NGL')
        cmd=f'cp {pwd}/kinasechains_renumber_uniprot/*.pdb.gz {pwd}/static/kinasechainsNGL/'
        subprocess.call(cmd,shell=True)
        os.chdir(f'{pwd}/static/kinasechainsNGL')
        cmd='gunzip -f *.pdb.gz'
        subprocess.call(cmd,shell=True)
        os.chdir(f'{pwd}')
