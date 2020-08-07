#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:44:28 2020

@author: vivekmodi
"""
import os, subprocess

def zip_coordinate_file(pwd,pdb_gene_dict,pdb_domain_dict):
    print('Compressing coordinate files...')
    os.chdir(f'{pwd}/static/downloads/uniprot-numbered')

    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if os.path.isfile(f'{pdbID}_uniprot.zip'):
           continue
        if not os.path.isdir(f'{pdbID}_uniprot'):
            cmd=(f'mkdir {pdbID}_uniprot')
            subprocess.call(cmd,shell=True)
        
       
        cmd=(f'mv {pdbs}.pdb {pdbID}_uniprot')
        subprocess.call(cmd,shell=True)
        cmd=(f'mv {pdbs}.cif {pdbID}_uniprot')
        subprocess.call(cmd,shell=True)
        
    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if not os.path.isfile(f'{pdbID}_uniprot.zip'):
            cmd=(f'zip -r {pdbID}_uniprot.zip {pdbID}_uniprot')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)
    
    #Alignment numbered directory
    os.chdir(f'{pwd}/static/downloads/alignment-numbered')

    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if os.path.isfile(f'{pdbID}_align.zip'):
           continue
        if not os.path.isdir(f'{pdbID}_align'):
            cmd=(f'mkdir {pdbID}_align')
            subprocess.call(cmd,shell=True)
        
       
        cmd=(f'mv {pdbs}.pdb {pdbID}_align')
        subprocess.call(cmd,shell=True)
        cmd=(f'mv {pdbs}.cif {pdbID}_align')
        subprocess.call(cmd,shell=True)
        
    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if not os.path.isfile(f'{pdbID}_align.zip'):
            cmd=(f'zip -r {pdbID}_align.zip {pdbID}_align')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)
    
    
    #PDB numbered directory
    os.chdir(f'{pwd}/static/downloads/pdb-numbered')

    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if os.path.isfile(f'{pdbID}_pdb.zip'):
           continue
        if not os.path.isdir(f'{pdbID}_pdb'):
            cmd=(f'mkdir {pdbID}_pdb')
            subprocess.call(cmd,shell=True)
        
       
        cmd=(f'mv {pdbs}.pdb {pdbID}_pdb')
        subprocess.call(cmd,shell=True)
        cmd=(f'mv {pdbs}.cif {pdbID}_pdb')
        subprocess.call(cmd,shell=True)
        
    for pdbs in pdb_gene_dict:
        pdbID=pdbs[0:4]
        
        if not os.path.isfile(f'{pdbID}_pdb.zip'):
            cmd=(f'zip -r {pdbID}_pdb.zip {pdbID}_pdb')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)
    
    #Domain directory - PDB numbered
    os.chdir(f'{pwd}/static/downloads/domain-pdb-numbered')

    domainList=set(pdb_domain_dict.values())
    
    for domainName in domainList:
        for pdbs in pdb_gene_dict:
            if pdb_domain_dict[pdbs]==domainName:
                pdbID=pdbs[0:4]
                
                
                if os.path.isfile(f'{domainName}_pdb.zip'):
                   continue
                if not os.path.isdir(f'{domainName}_pdb'):
                    cmd=(f'mkdir {domainName}_pdb')
                    subprocess.call(cmd,shell=True)
                
               
                cmd=(f'mv {pdbs}.pdb {domainName}_pdb')
                subprocess.call(cmd,shell=True)
                cmd=(f'mv {pdbs}.cif {domainName}_pdb')
                subprocess.call(cmd,shell=True)
                
    for domainName in domainList:
        
        if not os.path.isfile(f'{domainName}_pdb.zip'):
            cmd=(f'zip -r {domainName}_pdb.zip {domainName}_pdb')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)
    
    #Domain directory - Uniprot numbered
    os.chdir(f'{pwd}/static/downloads/domain-uniprot-numbered')

    domainList=set(pdb_domain_dict.values())
    
    for domainName in domainList:
        for pdbs in pdb_gene_dict:
            if pdb_domain_dict[pdbs]==domainName:
                pdbID=pdbs[0:4]
                
                
                if os.path.isfile(f'{domainName}_uniprot.zip'):
                   continue
                if not os.path.isdir(f'{domainName}_uniprot'):
                    cmd=(f'mkdir {domainName}_uniprot')
                    subprocess.call(cmd,shell=True)
                
               
                cmd=(f'mv {pdbs}.pdb {domainName}_uniprot')
                subprocess.call(cmd,shell=True)
                cmd=(f'mv {pdbs}.cif {domainName}_uniprot')
                subprocess.call(cmd,shell=True)
                
    for domainName in domainList:
        
        if not os.path.isfile(f'{domainName}_uniprot.zip'):
            cmd=(f'zip -r {domainName}_uniprot.zip {domainName}_uniprot')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)
    
    #Domain directory - Alignment numbered
    os.chdir(f'{pwd}/static/downloads/domain-alignment-numbered')

    domainList=set(pdb_domain_dict.values())
    
    for domainName in domainList:
        for pdbs in pdb_gene_dict:
            if pdb_domain_dict[pdbs]==domainName:
                pdbID=pdbs[0:4]
                
                
                if os.path.isfile(f'{domainName}_align.zip'):
                   continue
                if not os.path.isdir(f'{domainName}_align'):
                    cmd=(f'mkdir {domainName}_align')
                    subprocess.call(cmd,shell=True)
                
               
                cmd=(f'mv {pdbs}.pdb {domainName}_align')
                subprocess.call(cmd,shell=True)
                cmd=(f'mv {pdbs}.cif {domainName}_align')
                subprocess.call(cmd,shell=True)
                
    for domainName in domainList:
        
        if not os.path.isfile(f'{domainName}_align.zip'):
            cmd=(f'zip -r {domainName}_align.zip {domainName}_align')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
    
    os.chdir(pwd)