#! /usr/bin/python

import math
import os
import numpy as np
from Bio import PDB
import gzip


def compute_phi(structure_,model_,chain_,prev_residue_,curr_residue_):          #The original variable names can not be assigned again, therefore copies of variables are created here
    phi=999.00
    if prev_residue_.has_id('C') and curr_residue_.has_id('N') and curr_residue_.has_id('CA') and curr_residue_.has_id('C') and (curr_residue_.id[1]==(prev_residue_.id[1]+1)):     #The last condition is required to make sure that the residues are consecutive
        prev_c=structure_[model_.id][chain_.id][prev_residue_.id]['C'].get_vector()
        curr_n=structure_[model_.id][chain_.id][curr_residue_.id]['N'].get_vector()
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        curr_c=structure_[model_.id][chain_.id][curr_residue_.id]['C'].get_vector()
        phi=round(math.degrees(PDB.calc_dihedral(prev_c,curr_n,curr_ca,curr_c)),2)
    return phi

def compute_psi(structure_,model_,chain_,curr_residue_,next_residue_):
    psi=999.00
    if curr_residue_.has_id('N') and curr_residue_.has_id('CA') and curr_residue_.has_id('C') and next_residue_.has_id('N') and (next_residue_.id[1]==(curr_residue_.id[1]+1)):     #The last condition is required to make sure that the residues are consecutive
        curr_n=structure_[model_.id][chain_.id][curr_residue_.id]['N'].get_vector()
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        curr_c=structure_[model_.id][chain_.id][curr_residue_.id]['C'].get_vector()
        next_n=structure_[model_.id][chain_.id][next_residue_.id]['N'].get_vector()
        psi=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_c,next_n)),2)
    return psi

def compute_omega(structure_,model_,chain_,prev_residue_,curr_residue_):
    omega=999.00
    if prev_residue_.has_id('CA') and prev_residue_.has_id('C') and curr_residue_.has_id('N') and curr_residue_.has_id('CA') and (curr_residue_.id[1]==(prev_residue_.id[1]+1)):     #The last condition is required to make sure that the residues are consecutive
        prev_ca=structure_[model_.id][chain_.id][prev_residue_.id]['CA'].get_vector()
        prev_c=structure_[model_.id][chain_.id][prev_residue_.id]['C'].get_vector()
        curr_n=structure_[model_.id][chain_.id][curr_residue_.id]['N'].get_vector()
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        omega=round(math.degrees(PDB.calc_dihedral(prev_ca,prev_c,curr_n,curr_ca)),2)
    return omega

def compute_chi1(structure_,model_,chain_,curr_residue_):
    chi1=999.00
    if curr_residue_.has_id('N') and curr_residue_.has_id('CA') and curr_residue_.has_id('CB'):
        curr_n=structure_[model_.id][chain_.id][curr_residue_.id]['N'].get_vector()
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        curr_cb=structure_[model_.id][chain_.id][curr_residue_.id]['CB'].get_vector()

        if  curr_residue_.has_id('CG') and (curr_residue_.resname=='ARG' or curr_residue_.resname=='ASN' or curr_residue_.resname=='ASP' or curr_residue_.resname=='GLN' or curr_residue_.resname=='GLU' or curr_residue_.resname=='HIS' or curr_residue_.resname=='LEU' or curr_residue_.resname=='LYS' or curr_residue_.resname=='MET' or curr_residue_.resname=='PHE' or curr_residue_.resname=='PRO' or curr_residue_.resname=='TRP' or curr_residue_.resname=='TYR'):
            curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()
            chi1=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_cb,curr_cg)),2)

        if curr_residue_.has_id('SG') and curr_residue_.resname=='CYS':
            curr_sg=structure_[model_.id][chain_.id][curr_residue_.id]['SG'].get_vector()
            chi1=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_cb,curr_sg)),2)

        if curr_residue_.has_id('CG1') and (curr_residue_.resname=='VAL' or curr_residue_.resname=='ILE'):
            curr_cg1=structure_[model_.id][chain_.id][curr_residue_.id]['CG1'].get_vector()
            chi1=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_cb,curr_cg1)),2)

        if  curr_residue_.has_id('OG') and curr_residue_.resname=='SER':
            curr_og=structure_[model_.id][chain_.id][curr_residue_.id]['OG'].get_vector()
            chi1=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_cb,curr_og)),2)

        if  curr_residue_.has_id('OG1') and curr_residue_.resname=='THR':
            curr_og1=structure_[model_.id][chain_.id][curr_residue_.id]['OG1'].get_vector()
            chi1=round(math.degrees(PDB.calc_dihedral(curr_n,curr_ca,curr_cb,curr_og1)),2)
    return chi1

def compute_chi2(structure_,model_,chain_,curr_residue_):
    chi2=999.00
    if curr_residue_.has_id('CA') and curr_residue_.has_id('CB') and curr_residue_.has_id('CG'):
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        curr_cb=structure_[model_.id][chain_.id][curr_residue_.id]['CB'].get_vector()
        curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()

        if  curr_residue_.has_id('CD') and (curr_residue_.resname=='ARG' or curr_residue_.resname=='GLN' or curr_residue_.resname=='GLU' or curr_residue_.resname=='LYS' or curr_residue_.resname=='PRO'):
            curr_cd=structure_[model_.id][chain_.id][curr_residue_.id]['CD'].get_vector()
            chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg,curr_cd)),2)

        if curr_residue_.has_id('OD1') and (curr_residue_.resname=='ASN' or curr_residue_.resname=='ASP'):
            curr_od1=structure_[model_.id][chain_.id][curr_residue_.id]['OD1'].get_vector()
            chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg,curr_od1)),2)

        if curr_residue_.has_id('ND1') and curr_residue_.resname=='HIS':
            curr_nd1=structure_[model_.id][chain_.id][curr_residue_.id]['ND1'].get_vector()
            chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg,curr_nd1)),2)

        if curr_residue_.has_id('CD1') and (curr_residue_.resname=='LEU' or curr_residue_.resname=='PHE' or curr_residue_.resname=='TRP' or curr_residue_.resname=='TYR'):
            curr_cd1=structure_[model_.id][chain_.id][curr_residue_.id]['CD1'].get_vector()
            chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg,curr_cd1)),2)

        if curr_residue_.has_id('SD') and curr_residue_.resname=='MET':
            curr_sd=structure_[model_.id][chain_.id][curr_residue_.id]['SD'].get_vector()
            chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg,curr_sd)),2)

    if curr_residue_.has_id('CA') and curr_residue_.has_id('CB') and curr_residue_.has_id('CG1') and curr_residue_.has_id('CD1') and curr_residue_.resname=='ILE':
        curr_ca=structure_[model_.id][chain_.id][curr_residue_.id]['CA'].get_vector()
        curr_cb=structure_[model_.id][chain_.id][curr_residue_.id]['CB'].get_vector()
        curr_cg1=structure_[model_.id][chain_.id][curr_residue_.id]['CG1'].get_vector()
        curr_cd1=structure_[model_.id][chain_.id][curr_residue_.id]['CD1'].get_vector()
        chi2=round(math.degrees(PDB.calc_dihedral(curr_ca,curr_cb,curr_cg1,curr_cd1)),2)
    return chi2

def compute_chi3(structure_,model_,chain_,curr_residue_):
    chi3=999.00
    if curr_residue_.has_id('CB') and curr_residue_.has_id('CG') and curr_residue_.has_id('CD'):
        curr_cb=structure_[model_.id][chain_.id][curr_residue_.id]['CB'].get_vector()
        curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()
        curr_cd=structure_[model_.id][chain_.id][curr_residue_.id]['CD'].get_vector()

        if curr_residue_.has_id('NE') and curr_residue_.resname=='ARG':
            curr_ne=structure_[model_.id][chain_.id][curr_residue_.id]['NE'].get_vector()
            chi3=round(math.degrees(PDB.calc_dihedral(curr_cb,curr_cg,curr_cd,curr_ne)),2)

        if curr_residue_.has_id('OE1') and (curr_residue_.resname=='GLN' or curr_residue_.resname=='GLU'):
            curr_oe1=structure_[model_.id][chain_.id][curr_residue_.id]['OE1'].get_vector()
            chi3=round(math.degrees(PDB.calc_dihedral(curr_cb,curr_cg,curr_cd,curr_oe1)),2)

        if curr_residue_.has_id('CE') and curr_residue_.resname=='LYS':
            curr_ce=structure_[model_.id][chain_.id][curr_residue_.id]['CE'].get_vector()
            chi3=round(math.degrees(PDB.calc_dihedral(curr_cb,curr_cg,curr_cd,curr_ce)),2)

    if curr_residue_.has_id('CB') and curr_residue_.has_id('CG') and curr_residue_.has_id('SD') and curr_residue_.has_id('CE') and curr_residue_.resname=='MET':
        curr_cb=structure_[model_.id][chain_.id][curr_residue_.id]['CB'].get_vector()
        curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()
        curr_sd=structure_[model_.id][chain_.id][curr_residue_.id]['SD'].get_vector()
        curr_ce=structure_[model_.id][chain_.id][curr_residue_.id]['CE'].get_vector()
        chi3=round(math.degrees(PDB.calc_dihedral(curr_cb,curr_cg,curr_sd,curr_ce)),2)
    return chi3

def compute_chi4(structure_,model_,chain_,curr_residue_):
    chi4=999.00
    if curr_residue_.has_id('CG') and curr_residue_.has_id('CD') and curr_residue_.has_id('NE'):
        curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()
        curr_cd=structure_[model_.id][chain_.id][curr_residue_.id]['CD'].get_vector()
        curr_ne=structure_[model_.id][chain_.id][curr_residue_.id]['NE'].get_vector()

        if curr_residue_.has_id('CZ') and curr_residue_.resname=='ARG':
            curr_cz=structure_[model_.id][chain_.id][curr_residue_.id]['CZ'].get_vector()
            chi4=round(math.degrees(PDB.calc_dihedral(curr_cg,curr_cd,curr_ne,curr_cz)),2)

    if curr_residue_.has_id('CG') and curr_residue_.has_id('CD') and curr_residue_.has_id('CE'):
        curr_cg=structure_[model_.id][chain_.id][curr_residue_.id]['CG'].get_vector()
        curr_cd=structure_[model_.id][chain_.id][curr_residue_.id]['CD'].get_vector()
        curr_ce=structure_[model_.id][chain_.id][curr_residue_.id]['CE'].get_vector()

        if curr_residue_.has_id('NZ') and curr_residue_.resname=='LYS':
            curr_nz=structure_[model_.id][chain_.id][curr_residue_.id]['NZ'].get_vector()
            chi4=round(math.degrees(PDB.calc_dihedral(curr_cg,curr_cd,curr_ce,curr_nz)),2)

    return chi4

def compute_dihedrals(pwd,df):
    ignoremodified=open(f'{pwd}/List_modified_aminoacid.txt','r')
    log=open(f'{pwd}/kinasepml.log','a')
    
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        pdbfilename=pdbs
        if os.path.isfile(f'{pwd}/kinasechains_dihedrals/{pdbfilename}.dih'):
            continue
        try:
            handle=gzip.open(f'{pwd}/kinasechains_renumber_uniprot/{pdbfilename}.cif.gz','rt')
        except:
            log.write(f'dihedrals: File does not exist {pdbfilename}.cif.gz\n')
            continue
        
        fhandle_output=open((pwd+'/kinasechains_dihedrals/'+pdbfilename+'.dih'),'w')

        parser=PDB.MMCIFParser()
       
        
        structure=parser.get_structure('PDB',handle)
        for model in structure:
            for chain in model:
                first=1
                for residue in chain:
                    ignoremodified.seek(0)
                    if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():
                        #if residue.id[0]!=' ' and residue.id[0]!='H_PTR' and residue.id[0]!='H_TPO' and residue.id[0]!='H_SEP':         #If residue is not protein ' ' or PTR, TPO, SEP then skip
                        #    continue
        
                        if first==1:        #The 'first' blocks are required to assign first and second residue to variables
                            prev_residue=residue
                            first=2
                            continue
        
                        if first==2:        #This block computes psi dihedral for first residue
                            curr_residue=residue
                            psi=compute_psi(structure,model,chain,prev_residue,curr_residue)
                            chi1=compute_chi1(structure,model,chain,prev_residue)
                            chi2=compute_chi2(structure,model,chain,prev_residue)
                            chi3=compute_chi3(structure,model,chain,prev_residue)       #right now works only for ARG, GLN, GLU, LYS and MET
                            chi4=compute_chi4(structure,model,chain,prev_residue)       #right now works only for ARG, LYS
        
                            fhandle_output.write(f'{pdbfilename[:-4]} {model.id} {chain.id} {prev_residue.id[1]} {prev_residue.resname} 999.00 {psi} 999.00 {chi1} {chi2} {chi3} {chi4}\n')
                            first=3
                            continue
        
                        if first==3:        #This block computes phi and psi dihedrals from second residue onward. At anytime in the block we have three residue variables assigned.
                            next_residue=residue
                            phi=compute_phi(structure,model,chain,prev_residue,curr_residue)
                            psi=compute_psi(structure,model,chain,curr_residue,next_residue)
                            omega=compute_omega(structure,model,chain,prev_residue,curr_residue)
                            chi1=compute_chi1(structure,model,chain,curr_residue)
                            chi2=compute_chi2(structure,model,chain,curr_residue)
                            chi3=compute_chi3(structure,model,chain,curr_residue)       #right now works only for ARG, GLN, GLU, LYS and MET
                            chi4=compute_chi4(structure,model,chain,curr_residue)       #right now works only for ARG, LYS
        
                            fhandle_output.write(f'{pdbfilename[0:4]} {model.id} {chain.id} {curr_residue.id[1]} {curr_residue.resname} {phi} {psi} {omega} {chi1} {chi2} {chi3} {chi4}\n')
        
     #                       if curr_residue.id[1]==int(xdfg):
     #                           xdfg_phi=phi;xdfg_psi=psi
     #                       if curr_residue.id[1]==int(dfg_asp):
     #                           dfg_asp_phi=phi;dfg_asp_psi=psi
     #                       if curr_residue.id[1]==int(dfg_phe):
     #                           dfg_phe_phi=phi;dfg_phe_psi=psi;dfg_phe_chi1=chi1
     #                           if dfg_phe_chi1<0:
     #                               dfg_phe_chi1=np.round((dfg_phe_chi1+360),2)
        
                            prev_residue=curr_residue
                            curr_residue=next_residue       #update residue variables
    
                if first==3:        #This block computes phi dihedral for the last residue
                    phi=compute_phi(structure,model,chain,prev_residue,curr_residue)
                    omega=compute_omega(structure,model,chain,prev_residue,curr_residue)
                    chi1=compute_chi1(structure,model,chain,curr_residue)
                    chi2=compute_chi2(structure,model,chain,curr_residue)
                    chi3=compute_chi3(structure,model,chain,curr_residue)       #right now works only for ARG, GLN, GLU, LYS and MET
                    chi4=compute_chi4(structure,model,chain,curr_residue)       #right now works only for ARG, LYS
    
                    fhandle_output.write(f'{pdbfilename[0:4]} {model.id} {chain.id} {curr_residue.id[1]} {curr_residue.resname} {phi} 999.00 {omega} {chi1} {chi2} {chi3} {chi4}\n')
    
    log.close()
    return

