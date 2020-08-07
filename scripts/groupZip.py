#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 10:04:57 2020

@author: vivekmodi
"""
import subprocess, os

def groupZip(pwd,pdb_group_dict,pdb_domain_dict):
    print('Creating compressed directories for files of each group...')
    groupList=set(pdb_group_dict.values())
    
    for groups in groupList:
        if os.path.isfile(f'{pwd}/static/downloads/groupZip/{groups}.zip'):
            continue
        if not os.path.isdir(f'{pwd}/static/downloads/groupZip/{groups}'):
            cmd=(f'mkdir {pwd}/static/downloads/groupZip/{groups}')
            subprocess.call(cmd,shell=True)
        for pdbs in pdb_domain_dict:
            if groups==pdb_group_dict[pdbs]:
                cmd=(f'mv {pwd}/static/downloads/groupZip/{pdbs}.cif {pwd}/static/downloads/groupZip/{groups}/{pdb_domain_dict[pdbs]}-{pdbs}.cif')
                subprocess.call(cmd,shell=True)
        
        os.chdir(f'{pwd}/static/downloads/groupZip/')
        cmd=(f'zip -r {groups}.zip {groups}')
        subprocess.call(cmd,shell=True)
        os.chdir(f'{pwd}')
        
        