#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:12:50 2020

@author: vivekmodi
"""
import subprocess

def representativeSuperpositionGene(pwd,pdb_domain_dict):
    print('Creating superpositions for representative structures of genes...')
    fhandleFile=open(f'{pwd}/Representative_structures.tab','r')
    domainList=sorted(set(pdb_domain_dict.values()))
    
    for domainName in domainList:
        fhandlePymol=open(f'{pwd}/static/downloads/representativeSuperpositionDomain/{domainName}.pml','w')
        fhandleFile.seek(0)
        for lines in fhandleFile:
            lines=lines.strip();lines=lines.split()
            
            if lines[0]==domainName:
                pdbs=lines[3]
                fhandlePymol.write(f'load {pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz\n')
        fhandlePymol.write(f'alignto {pdbs}\n')
        
        fhandleFile.seek(0)
        for lines in fhandleFile:
            lines=lines.strip();lines=lines.split()
            
            if lines[0]==domainName:
                pdbs=lines[3]
                fhandlePymol.write(f'save {pwd}/static/downloads/representativeSuperpositionDomain/{pdbs}.pdb, {pdbs}\n')
        fhandlePymol.close()
        
        cmd=(f'pymol -c {pwd}/static/downloads/representativeSuperpositionDomain/{domainName}.pml')        
        subprocess.call(cmd,shell=True)
    
    fhandleFile.close()
    
        
    
     
        
    
    