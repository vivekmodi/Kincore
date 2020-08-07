#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:44:11 2020

@author: vivekmodi
"""
from datetime import datetime
import os, subprocess

def create_dirs(pwd):
         
    #today=str(datetime.now())[0:10].strip()
    #print("Date:"+today)
    print("Creating new directories...")
    #update=(pwd+"/update-"+today)
    #if not os.path.isdir(update):
    #    cmd=("mkdir "+update)
    #    subprocess.call(cmd,shell=True)
    
    pdbaa_psiblast_dir=(pwd+"/"+"pdbaa_psiblast_dir")
    kinasecifs=(pwd+"/"+"kinasecifs")
    #kinasebioassm=(pwd+"/"+"kinasebioassm")
    kinasesifts=(pwd+"/"+"kinasesifts")
    kinasechains=(pwd+"/"+"kinasechains")
    kinasechains_renumber_uniprot=(pwd+"/"+"kinasechains_renumber_uniprot")
    kinasechains_renumber_alignment=(pwd+"/"+"kinasechains_renumber_alignment")
    pymol_sessions=(pwd+"/"+"static/downloads/pymol-uniprot")
    kinasechains_dihedrals=(pwd+'/'+'kinasechains_dihedrals')
    formattedSeq=(pwd+'/'+'formattedSeq')
    pymolChainSuperposition=(pwd+'/'+'pymolChainSuperposition')
    
    # print(pdbaa_psiblast_dir)
    # print(kinasecifs)
    # #print(kinasebioassm)
    # print(kinasesifts)
    # print(kinasechains)
    # print(kinasechains_renumber_uniprot)
    # print(kinasechains_renumber_alignment)
    # print(pymol_sessions)
    # print(kinasechains_dihedrals)
    # print(formattedSeq)

    if not os.path.isdir(pdbaa_psiblast_dir):
        cmd=("mkdir "+pdbaa_psiblast_dir)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasecifs):
        cmd=("mkdir "+kinasecifs)
        subprocess.call(cmd,shell=True)
    #if not os.path.isdir(kinasebioassm):
    #    cmd=("mkdir "+kinasebioassm)
    #    subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasesifts):
        cmd=("mkdir "+kinasesifts)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasechains):
        cmd=("mkdir "+kinasechains)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasechains_renumber_uniprot):
        cmd=("mkdir "+kinasechains_renumber_uniprot)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasechains_renumber_alignment):
        cmd=("mkdir "+kinasechains_renumber_alignment)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(pymol_sessions):
        cmd=("mkdir "+pymol_sessions)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasechains_dihedrals):
        cmd=("mkdir "+kinasechains_dihedrals)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(formattedSeq):
        cmd=("mkdir "+formattedSeq)
        subprocess.call(cmd,shell=True)
    #if not os.path.isdir(pymolChainSuperposition):
    #    cmd=("mkdir "+pymolChainSuperposition)
    #    subprocess.call(cmd,shell=True)
    #return(pdbaa_psiblast_dir,kinasecifs,kinasesifts,kinasechains,kinasechains_renumber_uniprot,kinasechains_renumber_alignment,pymol_sessions,kinasechains_dihedrals)
    return