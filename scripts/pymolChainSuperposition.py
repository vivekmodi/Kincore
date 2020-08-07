#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 20:12:51 2020

@author: vivekmodi
"""
import subprocess

def pymolChainSuperposition(pwd,pdb_gene_dict):
    print('Creating chain superpositions...')
    pdbChains=dict()
    for pdbs in pdb_gene_dict:
        pdbChains[pdbs[0:4]]=list()
    for pdbs in pdb_gene_dict:
        pdbChains[pdbs[0:4]].append(pdbs)
        
    for items in pdbChains:
        fhandle_pymol=open(pwd+'/pymolChainSuperposition/'+items+'.pml','w')
        for chains in sorted(pdbChains[items],reverse=True):
            fhandle_pymol.write(f'load {pwd}/kinasechains_renumber_uniprot/{chains}.pdb.gz\n')
        fhandle_pymol.write(f'alignto {chains}\n')
        fhandle_pymol.write(f'save {pwd}/pymolChainSuperposition/{items}.pdb')
        fhandle_pymol.close()

    for items in pdbChains:
        cmd=(f'pymol -c {pwd}/pymolChainSuperposition/{items}.pml')        
        subprocess.call(cmd,shell=True)