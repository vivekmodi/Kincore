#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 15:30:26 2020

@author: vivekmodi
"""
import gzip
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def get_seq_from_cif(pwd,df):
    print('Reading sequence from Uniprot numbered .cif file...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        handle=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz','rt')    #will contain only the Uniprot residues, not 
        for record in SeqIO.parse(handle, "cif-atom"):
            sequence=list(str(record.seq))
            for letters in sequence[::-1]:    #fix for a probable bug in Biopython, it introduces a trail  of 'X' in some cases at the end of the sequence
                if letters=='X':
                    sequence.pop()
                else:
                    break
            
            df.at[i,'StrSeq']=''.join(sequence)
        
        #The following part is to get SeqRes record sequence, this is the sequence of the construct
        #It should be done using .pdb files. If needed, download the files in future and uncomment the following code (not tested)
        
        #print('Reading sequence from original pdb file to get SeqRes...')
        df.at[i,'SeqRes']='XXXX'
        #pdb=pdbs[0:4].lower()
        #chainid=pdbs[4]
        #handle=gzip.open(f'{pwd}/kinasecifs/{pdb}.pdb.gz','rt')    #This sequence should come from original files not renumbered
        #for record in SeqIO.parse(handle,'pdb-seqres'):   #Biopython uses different chain id here, not from the same column which we use everywhere else, fix in future
        #    if record.annotations['chain']==chainid:
        #        df.at[i,'SeqRes']=str(record.seq)
            
    return df
        
        