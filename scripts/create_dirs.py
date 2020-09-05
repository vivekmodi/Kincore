#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:44:11 2020

@author: vivekmodi
"""

import os, subprocess

def create_dirs(pwd):
    print("Creating new directories...")
    pdbaa_psiblast_dir=(pwd+"/"+"pdbaa_psiblast_dir")
    kinasecifs=(pwd+"/"+"kinasecifs")
    kinasesifts=(pwd+"/"+"kinasesifts")
    kinasechains=(pwd+"/"+"kinasechains")
    kinasechains_renumber_uniprot=(pwd+"/"+"kinasechains_renumber_uniprot")
    kinasechains_renumber_alignment=(pwd+"/"+"kinasechains_renumber_alignment")
    kinasechains_dihedrals=(pwd+'/'+'kinasechains_dihedrals')
    formattedSeq=(pwd+'/'+'formattedSeq')
    mtzfiles=(pwd+'/'+'mtz_files')
    edia=(pwd+'/'+'mtz_files/Edia_out')
    static=(pwd+'/'+'static')
    ngl=(pwd+'/static/kinasechainsNGL')

    if not os.path.isdir(pdbaa_psiblast_dir):
        cmd=("mkdir "+pdbaa_psiblast_dir)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(kinasecifs):
        cmd=("mkdir "+kinasecifs)
        subprocess.call(cmd,shell=True)
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
    if not os.path.isdir(kinasechains_dihedrals):
        cmd=("mkdir "+kinasechains_dihedrals)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(formattedSeq):
        cmd=("mkdir "+formattedSeq)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(mtzfiles):
        cmd=("mkdir "+mtzfiles)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(edia):
        cmd=("mkdir "+edia)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(static):
        cmd=("mkdir "+static)
        subprocess.call(cmd,shell=True)
    if not os.path.isdir(ngl):
        cmd=("mkdir "+ngl)
        subprocess.call(cmd,shell=True)
    return