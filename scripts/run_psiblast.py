#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:50:36 2020

@author: vivekmodi
"""
#Usage from commandline python run_psiblast.py AURKA.fasta pdbaa.db AURKA.check AURKA.pdbaa.xml
#Usage from another script run_psiblast(query, database, out_pssm, out)

import sys
from Bio.Blast.Applications import NcbipsiblastCommandline

def run_psiblast(query, database, out_pssm, out):
    print("Running psiblast...")
    psiblast_cline=NcbipsiblastCommandline(query=query,db=database,out_pssm=out_pssm, out=out,num_iterations=3,num_alignments=10000,num_descriptions=10000)
    psiblast_cline()

    psiblast_cline2=NcbipsiblastCommandline(in_pssm=out_pssm, db=database, out=out, num_alignments=20000,num_descriptions=20000,outfmt=5)
    psiblast_cline2()
    
if __name__ == '__main__':
    query=sys.argv[1]
    database=sys.argv[2]
    out_pssm=sys.argv[3]
    out=sys.argv[4]
    run_psiblast(query, database, out_pssm, out)