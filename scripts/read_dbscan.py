#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:43:03 2020

@author: vivekmodi
"""

def read_dbscan(pdb_gene_dict):
    fhandle_dfgin=open('dbscan-DFGin-clusters.tab','r')
    fhandle_dfginter=open('dbscan-DFGinter-clusters.tab','r')
    fhandle_dfgout=open('dbscan-DFGout-clusters.tab','r')
    pdb_dbscan_dict=dict()
    for pdbs in pdb_gene_dict:
        fhandle_dfgin.seek(0)
        fhandle_dfginter.seek(0)
        fhandle_dfgout.seek(0)
        pdb_dbscan_dict[pdbs]='NA'
        
        for lines in fhandle_dfgin:
            lines=lines.strip();lines=lines.split()
            if pdbs==lines[1]:
                
                if int(lines[32])==1:
                    pdb_dbscan_dict[pdbs]='BLAplus'
                    continue
                if int(lines[32])==2:
                    pdb_dbscan_dict[pdbs]='BLBplus'
                    continue
                if int(lines[32])==3:
                    pdb_dbscan_dict[pdbs]='BLBtrans'
                    continue
                if int(lines[32])==4:
                    pdb_dbscan_dict[pdbs]='BLAminus'
                    continue
                if int(lines[32])==5:
                    pdb_dbscan_dict[pdbs]='BLBminus'
                    continue
                if int(lines[32])==6:
                    pdb_dbscan_dict[pdbs]='ABAminus'
                    continue
                
        for lines in fhandle_dfginter:
            lines=lines.strip();lines=lines.split()
            if pdbs==lines[1]:
                if int(lines[32])==1:
                    pdb_dbscan_dict[pdbs]='BABtrans'
                    
        for lines in fhandle_dfgout:
            lines=lines.strip();lines=lines.split()
            if pdbs==lines[1]:
                if int(lines[32])==1 or int(lines[32])==2:
                    pdb_dbscan_dict[pdbs]='BBAminus'
    
    fhandle_dfgin.close()
    fhandle_dfginter.close()
    fhandle_dfgout.close()
    return pdb_dbscan_dict