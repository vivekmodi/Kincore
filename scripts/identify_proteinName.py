#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 12:49:53 2020

@author: vivekmodi
"""

def identify_proteinName(domain_uniprot_dict):
    uniprot_proteinName_dict=dict()
    fhandle=open('Kinases-uniprot-protein-name.tab','r')
    for uniprot in domain_uniprot_dict:
        uniprotacc=domain_uniprot_dict[uniprot]
        fhandle.seek(0)
        for lines in fhandle:
            lines=lines.strip();lines=lines.split(';')
            if lines[1]==uniprotacc:
                uniprot_proteinName_dict[uniprotacc]=lines[2]
    fhandle.close()
    return uniprot_proteinName_dict
        
        