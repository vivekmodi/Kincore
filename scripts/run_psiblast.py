#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 09:59:38 2020

@author: vivekmodi
"""

# Usage from commandline $python3 run_psiblast.py /home/psiblast_dir pdbaa AurkaPsiblastIter6PSSM.asn AURKA.pdbaa.xml

import sys
from Bio.Blast.Applications import NcbipsiblastCommandline

def run_psiblast(path,database,in_pssm,out_xml):
    psiblast_cline=NcbipsiblastCommandline(db=f'{path}/{database}.db',in_pssm=f'{path}/{in_pssm}', out=f'{path}/{out_xml}',max_target_seqs=10000,outfmt=5)
    psiblast_cline()   
    
    
if __name__ == '__main__':
    path=sys.argv[1]
    database=sys.argv[2]
    in_pssm=sys.argv[3]
    out_xml=sys.argv[4]
    run_psiblast(path,database,in_pssm,out_xml)