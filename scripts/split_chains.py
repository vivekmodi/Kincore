#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:17:07 2020

@author: vivekmodi
"""

#Usage from commandline python split_chains.py kinasecifs, kinasechains, 1xyzA    (prints the chain id provided in a separate file)
#Usage from another script split_chains(kinasecifs, kinasechains, pdbs)

import os, gzip, subprocess
from Bio import PDB

############# FIX (The files which have disordered atoms give an error while copying while printing separate chains. This block is a fix copied from https://github.com/biopython/biopython/issues/455)

def get_unpacked_list(self):
#     """
#     Returns all atoms from the residue,
#     in case of disordered, keep only first alt loc and remove the alt-loc tag
#     """
    atom_list = self.get_list()
    undisordered_atom_list = []
    for atom in atom_list:
        if atom.is_disordered():
            atom.altloc=" "
            undisordered_atom_list.append(atom)
        else:
            undisordered_atom_list.append(atom)
    return undisordered_atom_list

PDB.Residue.Residue.get_unpacked_list = get_unpacked_list

###########################


def split_chains(pwd, df):
    kinasecifs=f'{pwd}/kinasecifs'
    kinasechains=f'{pwd}/kinasechains'
    print('Splitting pdbs and printing chains in separate files...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        filename=(kinasecifs+"/"+pdbs[0:4].lower()+".cif.gz")
        if not os.path.isfile(filename):
            print("Error: Function split_chains: file does not exist:"+filename+"\n")
            return
        
        
    
        chain_filename=(kinasechains+"/"+pdbs+".pdb.gz")
        if not os.path.isfile(chain_filename):
            #print("Splitting pdbs and printing chains in separate files:"+pdbs)
            handle=gzip.open(filename,"rt")
            parser=PDB.MMCIFParser(QUIET=True)
            structure=parser.get_structure(pdbs[0:4],handle)
            df.at[i,'Date']=structure.header['deposition_date']
            io=PDB.PDBIO()
            chain_id=pdbs[4]    #Take the chain id from pdbs
            for model in structure:
                for chain in model:
                    if(chain.id==chain_id):
                        io.set_structure(chain)
                        chain_filename=(kinasechains+"/"+structure.id+chain.id+".pdb")             #Output in PDB format
                        io.save(chain_filename)        
                        #cmd=(f'cp {kinasechains}/{pdbs[0:4]}{pdbs[4]}.pdb {pwd}/static/downloads/pdb-numbered')
                        #print(cmd)
                        #subprocess.call(cmd,shell=True)
                        cmd=('gzip -f '+chain_filename)
                        subprocess.call(cmd, shell=True)
            #handle.close()
            
        chain_filename=(kinasechains+"/"+pdbs+".cif.gz")
        if not os.path.isfile(chain_filename):
            handle=gzip.open(filename,"rt")
            parser=PDB.MMCIFParser(QUIET=True)
            structure=parser.get_structure(pdbs[0:4],handle)
            io=PDB.MMCIFIO()
            chain_id=pdbs[4]    #Take the chain id from pdbs    
            for model in structure:
                for chain in model:
                    if(chain.id==chain_id):
                        io.set_structure(chain)
                        chain_filename=(kinasechains+"/"+structure.id+chain.id+".cif")             #Output in MMCIF format
                        io.save(chain_filename)     
                        #cmd=(f'cp {kinasechains}/{pdbs[0:4]}{pdbs[4]}.cif {pwd}/static/downloads/pdb-numbered')
                        #print(cmd)
                        #subprocess.call(cmd,shell=True)
                        cmd=('gzip -f '+chain_filename)
                        subprocess.call(cmd, shell=True)
            #handle.close()
#if __name__ == '__main__':
#    kinasecifs=sys.argv[1]
#    kinasechains=sys.argv[2]
#    pdbs=sys.argv[4]
#    split_chains(kinasecifs, kinasechains, pdbs)