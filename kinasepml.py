#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:18:13 2020

@author: vivekmodi
"""
import os, sys, subprocess, numpy as np, pandas as pd
from datetime import datetime


sys.path.append(os.getcwd()+'/scripts')
from create_dirs import create_dirs
from Bio.Blast.Applications import NcbipsiblastCommandline
from download_pdbaa import download_pdbaa
from create_blastdb import create_blastdb
from run_psiblast import run_psiblast
from read_psiblast import read_psiblast
from create_gene_dict import domain_dict
from create_gene_dict import gene_dict
from kinase_statistics import statistics
from download_cifs import download_cifs
from download_sifts import download_sifts
from split_chains import split_chains
from parse_sifts import parse_sifts
from renumber_by_uniprot import renumber_by_uniprot
from renumber_by_alignment import renumber_by_alignment
from pymol_script import pymol_genes, pymol_genes_script
from dihedrals import compute_dihedrals
from read_dihedrals import read_dihedrals
from read_mmcif_header import read_header
from compute_dist_atoms import compute_distance
from assign_labels import spatial_labels, dihedral_labels, dihedral_labels_no_chi
from extract_ligands import extract_ligands
from chain_break import chain_break
from gene_synonym import gene_synonym
from identify_mutation import identify_mutation
from representative_structure import representative_structure, representative_structure_ligand
from pseudokinases import identify_pseudokinases
from read_dbscan import read_dbscan
from identify_proteinName import identify_proteinName
#from ngl_script import write_ngl_script
from seqFormatHtml import format_seq_html
from seqFormatHtml import format_seq_text
from pymolChainSuperposition import pymolChainSuperposition
from chainColor import chainColor
from pymol_ligands import pymol_ligands, pymol_ligands_session_compress, pymol_ligands_scripts, pymol_ligands_scripts_compress
from zip_coordinate_files import zip_coordinate_file
from groupZip import groupZip
from PymolSuperpositionDomain import representativeSuperpositionGene
from createPymolSession import subListPymolSession

today=str(datetime.now())[0:10].strip()


  
if __name__ == '__main__':
    pwd=os.getcwd()
    print("Present working directory: "+pwd)
    
    (pdbaa_psiblast_dir,kinasecifs,kinasesifts,kinasechains,kinasechains_renumber_uniprot,kinasechains_renumber_alignment,pymol_sessions,kinasechains_dihedrals)=create_dirs(pwd)
    
    #download_pdbaa(pdbaa_psiblast_dir)    
    
    #create_blastdb('pdbaa',pdbaa_psiblast_dir)
    query=(pwd+'/AURKA.fasta');database=(pdbaa_psiblast_dir+'/pdbaa.db');in_pssm=(pdbaa_psiblast_dir+'/AurkaPsiblastIter6PSSM.asn')
    out=(pdbaa_psiblast_dir+'/AURKA.pdbaa.xml')
    #run_psiblast(query, database, out_pssm, out)
    psiblast_cline=NcbipsiblastCommandline(db=database,in_pssm=in_pssm, out=out,num_alignments=10000,num_descriptions=10000,outfmt=5)
    psiblast_cline()    
    pdb_uniprot_dict=dict();seqres_start=dict();seqres_end=dict();pdb_gene_dict=dict();pdb_ligand_dict=dict();pdb_gene_synonym_dict=dict()
    pdb_reso_dict=dict();domain_break_dict=dict();loop_break_dict=dict();domain_pseudo_dict=dict();uniprot_proteinName_dict=dict()
    pdb_xphi_dict=dict();pdb_xpsi_dict=dict();pdb_aspphi_dict=dict();pdb_asppsi_dict=dict();pdb_aspchi1_dict=dict();pdb_aspchi2_dict=dict()
    pdb_phephi_dict=dict();pdb_phepsi_dict=dict();pdb_phechi1_dict=dict();pdb_phechi2_dict=dict();pdb_glyphi_dict=dict();pdb_glypsi_dict=dict()
    pdb_dfgmutation_dict=dict();pdb_chainmutation_dict=dict();pdb_phos_dict=dict()
    
    psiblast_result=(pdbaa_psiblast_dir+'/AURKA.pdbaa.xml');human_output=(pdbaa_psiblast_dir+'/psiblast-human.csv');nonhuman_output=(pdbaa_psiblast_dir+'/psiblast-nonhuman.csv')
    psiblast_error=(pdbaa_psiblast_dir+'/psiblast-error.log'); excluded_output=(pdbaa_psiblast_dir+'/psiblast_excluded.log')
    
    (pdb_uniprot_dict,seqres_start_dict,seqres_end_dict,pdb_dict)=read_psiblast(psiblast_result, human_output, nonhuman_output, excluded_output, psiblast_error)
    #To check the final number of chains, len(pdb_uniprot_dict)+excluded_entries=entries_psiblast_human - repeated_entries
    
        
    (domain_gene_dict,domain_group_dict,domain_uniprot_dict,domain_uniacc_dict,domainstart_dict,domainend_dict,\
     domain_alkres_dict,domain_alknum_dict,domain_rreres_dict,domain_rrenum_dict,domain_lrlres_dict,domain_lrlnum_dict,\
     domain_hrdres_dict,domain_hrdnum_dict,domain_dfgres_dict,domain_dfgnum_dict,domain_aperes_dict,domain_apenum_dict,domain_dlwres_dict,domain_dlwnum_dict,\
     domain_isrres_dict,domain_isrnum_dict)=domain_dict(pwd)
    
    (pdb_gene_dict,pdb_domain_dict,pdb_group_dict,pdb_uniacc_dict,pdb_domainstart_dict,pdb_domainend_dict,pdb_dict)=gene_dict(pwd,pdb_uniprot_dict,seqres_start_dict,seqres_end_dict,pdb_dict)
    (domain_pseudo_dict,pdb_dict)=identify_pseudokinases(domain_gene_dict,pdb_dict)
    
    #df=pd.DataFrame(pdb_dict.values(),columns=['PDBid','Protein','Uniprot','ChainLen','StrBegin','StrEnd','Species','Resolution','Method','Rvalue','FreeRvalue','Sequence','Gene','Domain','Group','UniAcc','DomainBegin','DomainEnd','Pseudo'])
    
    download_cifs(kinasecifs,pdb_domain_dict)
    download_sifts(kinasesifts,pdb_domain_dict)
    split_chains(pwd,kinasecifs,kinasechains,pdb_domain_dict)
    parse_sifts(kinasesifts,pdb_domain_dict)
    renumber_by_uniprot(pwd,kinasechains,kinasesifts,kinasechains_renumber_uniprot,pdb_domain_dict)
    renumber_by_alignment(pwd,kinasechains_renumber_uniprot,kinasechains_renumber_alignment,pdb_domain_dict)
    compute_dihedrals(kinasechains_renumber_uniprot,kinasechains_dihedrals,pdb_domain_dict)
    #zip_coordinate_file(pwd,pdb_gene_dict,pdb_domain_dict)
    #groupZip(pwd,pdb_group_dict,pdb_domain_dict)
    (pdb_xphi_dict,pdb_xpsi_dict,pdb_aspphi_dict,pdb_asppsi_dict,pdb_aspchi1_dict,pdb_aspchi2_dict,pdb_phephi_dict,pdb_phepsi_dict,\
    pdb_phechi1_dict,pdb_phechi2_dict,pdb_glyphi_dict,pdb_glypsi_dict)=read_dihedrals(kinasechains_dihedrals,pdb_domain_dict,domain_dfgnum_dict)
    (pdb_reso_dict)=read_header(pwd,kinasecifs,pdb_domain_dict)
    (domain_break_dict,loop_break_dict)=chain_break(pwd,kinasechains_renumber_uniprot,pdb_domain_dict,domainstart_dict,domainend_dict,domain_dfgnum_dict,domain_apenum_dict)
    
    pdb_gene_synonym_dict=gene_synonym(pwd,pdb_gene_dict)
    (pdb_dfgmutation_dict,pdb_chainmutation_dict,pdb_phos_dict,pdb_first_resi,pdb_last_resi)=identify_mutation(kinasesifts,pdb_domain_dict,domain_dfgnum_dict)
    
    uniprot_proteinName_dict=identify_proteinName(domain_uniprot_dict)
    
    #write_ngl_script(pwd,pdb_domain_dict,domain_dfgnum_dict)
   
    #pymolChainSuperposition(pwd,pdb_gene_dict)
    #representativeSuperpositionGene(pwd,pdb_domain_dict)
    format_seq_text(pwd,pdb_gene_dict)
    format_seq_html(pwd,pdb_gene_dict)
    
    
    
    pdb_sb_dict=dict();pdb_pheglu4_dict=dict();pdb_phelys_dict=dict()
    fhandle_read_distances=open((pwd+'/Distances.tab'),'r')
    fhandle_append_distances=open((pwd+'/Distances.tab'),'a')
    print('Computing distance between residue pairs...')
    
    for pdbs in pdb_domain_dict:
        pdb_present=0
        fhandle_read_distances.seek(0)
        for lines in fhandle_read_distances:
            lines=lines.strip();lines=lines.split()
            if pdbs in lines[0]:
                pdb_sb_dict[pdbs]=float(lines[1]);pdb_pheglu4_dict[pdbs]=float(lines[2]);pdb_phelys_dict[pdbs]=float(lines[3])
                pdb_present=1
                break
                
        if pdb_present==0:
            pdbid=pdbs[0:4]
            chainid=pdbs[4]
            domain_name=pdb_domain_dict[pdbs]
            pdb_sb_dict[pdbs]=999;pdb_pheglu4_dict[pdbs]=999;pdb_phelys_dict[pdbs]=999
     
            pdb_sb_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_alknum_dict[domain_name]),'CB',int(domain_rrenum_dict[domain_name]),'CB')
                    
            if domain_dfgres_dict[domain_name]=='F' or domain_dfgres_dict[domain_name]=='R':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CZ',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CZ',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='L' or domain_dfgres_dict[domain_name]=='P' or domain_dfgres_dict[domain_name]=='N':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CG',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CG',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='M':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CE',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CE',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='S':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'OG',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'OG',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='H':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'NE2',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'NE2',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='V' or domain_dfgres_dict[domain_name]=='A':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CB',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CB',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='W':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CZ3',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'CZ3',int(domain_alknum_dict[domain_name]),'CA')
            if domain_dfgres_dict[domain_name]=='Y':
                pdb_pheglu4_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'OH',int(domain_rrenum_dict[domain_name])+4,'CA')
                pdb_phelys_dict[pdbs]=compute_distance(kinasechains_renumber_uniprot,pdbs,int(domain_dfgnum_dict[domain_name]),'OH',int(domain_alknum_dict[domain_name]),'CA')
            fhandle_append_distances.write(f'{pdbs} {pdb_sb_dict[pdbs]} {pdb_pheglu4_dict[pdbs]} {pdb_phelys_dict[pdbs]}\n')
    
    fhandle_read_distances.close()
    fhandle_append_distances.close()
    
    pdb_spatial_dict=dict();pdb_dihedral_dict=dict();pdb_dihedral_dict_no_chi1=dict();pdb_dbscan_dict=dict()
    pdb_spatial_dict=spatial_labels(pdb_pheglu4_dict,pdb_phelys_dict)  #Exclude structures with DFGmutations
    pdb_dihedral_dict=dihedral_labels(pdb_spatial_dict,pdb_xphi_dict,pdb_xpsi_dict,pdb_aspphi_dict,pdb_asppsi_dict,pdb_phephi_dict,pdb_phepsi_dict,pdb_phechi1_dict)
    
    pdb_dihedral_dict_no_chi1=dihedral_labels_no_chi(pdb_spatial_dict,pdb_xphi_dict,pdb_xpsi_dict,pdb_aspphi_dict,pdb_asppsi_dict,pdb_phephi_dict,pdb_phepsi_dict,pdb_phechi1_dict)
    pdb_dbscan_dict=read_dbscan(pdb_gene_dict)
    pdb_ligand_dict=extract_ligands(pwd,kinasechains,pdb_domain_dict)
              
    #pymol_ligands(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
    #pymol_ligands_session_compress(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
    #pymol_ligands_scripts(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
    #pymol_ligands_scripts_compress(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
   
    
    #pymol_genes(pwd,pymol_sessions,kinasechains_renumber_uniprot,pdb_domain_dict,domain_group_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
    #pymol_genes_script(pwd,pymol_sessions,kinasechains_renumber_uniprot,pdb_domain_dict,domain_group_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict)
    #cmd=('gzip -f '+pymol_sessions+'/*pse')     #Sessions have to zipped separetely to make sure that they are created before this command is executed
    #subprocess.call(cmd,shell=True)
    
    #representative_structure(pwd,loop_break_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict)
    #representative_structure_ligand(pwd,loop_break_dict,pdb_ligand_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_reso_dict)
    pdb_color_dict=chainColor(pdb_spatial_dict,pdb_dihedral_dict)
    
     #Create All Pymol sessions
    sessionDir='static/downloads/pymolSessions'
    #subListPymolSession(pwd,sessionDir,kinasechains_renumber_uniprot,domain_group_dict,pdb_group_dict,pdb_domain_dict,pdb_spatial_dict,pdb_dihedral_dict,pdb_ligand_dict,domain_dfgnum_dict,domain_apenum_dict,domain_break_dict,pdb_reso_dict,pdb_color_dict)

    fhandle_output=open((pwd+'/Output-'+today+'.tab'),'w')
    fhandle_databasefile=open((pwd+'/Database-input-'+today+'.tab'),'w')
    fhandle_output.write('1_PDB 2_Chain 3_Reso 4_Gene 5_Uniprot 6_Domain 7_Group 8_UniAcc 9_Pseudo 10_DomainStart 11_DomainEnd 12_DomainBreak 13_LoopBreak 14_X_Phi 15_X_Psi 16_Asp_Phi 17_Asp_Psi 18_Asp_Chi1 19_Asp_Chi2 20_Phe_Phi 21_Phe_Psi 22_Phi_Chi1 23_Phe_Chi2 24_Gly_Phi 25_Gly_Psi 26_SB 27_Phe_Glu 28_Phe_Lys 29_DFGMutation 30_ChainMutation 31_Dbscan 32_Spatial 33_Dihedral 34_Newlabel 35_Phosphorylation 36_Ligand 37_PheNum 38_Color 39_PDBFirstResi 40_PDBLastResi 41_ApeNum\n')
    fhandle_databasefile.write(f'1_Chain;2_Reso;3_Gene;4_Uniprot;5_Domain;6_Group;7_Uniacc;8_Spatial;9_Dihedral;10_Ligand;11_Synonym;12_Pseudo;13_DFGMutation;14_ChainMutation;15_Phosphorylation;16_ProteinName;17_PheNum;18_Color;19_DomainStart;20_DomainEnd;21_DomainBreak;22_LoopBreak;23_Pseudo;24_PDBFirstResi;25_PDBLastResi;26_APENum\n')
            
    for pdbs in sorted(pdb_domain_dict):
        #try:
            domain_name=pdb_domain_dict[pdbs]
            uniprotname=pdb_uniprot_dict[pdbs]
            fhandle_output.write(f'{pdbs[0:4]} {pdbs} {pdb_reso_dict[pdbs]} {pdb_gene_dict[pdbs]} {pdb_uniprot_dict[pdbs]} {pdb_domain_dict[pdbs]} {pdb_group_dict[pdbs]} {pdb_uniacc_dict[pdbs]} {domain_pseudo_dict[domain_name]} {domainstart_dict[domain_name]} {domainend_dict[domain_name]} {domain_break_dict[pdbs]} {loop_break_dict[pdbs]} {pdb_xphi_dict[pdbs]:.2f} {pdb_xpsi_dict[pdbs]:.2f} {pdb_aspphi_dict[pdbs]:.2f} {pdb_asppsi_dict[pdbs]:.2f} {pdb_aspchi1_dict[pdbs]:.2f} {pdb_aspchi2_dict[pdbs]:.2f} {pdb_phephi_dict[pdbs]:.2f} {pdb_phepsi_dict[pdbs]:.2f} {pdb_phechi1_dict[pdbs]:.2f} {pdb_phechi2_dict[pdbs]:.2f} {pdb_glyphi_dict[pdbs]:.2f} {pdb_glypsi_dict[pdbs]:.2f} {pdb_sb_dict[pdbs]:.2f} {pdb_pheglu4_dict[pdbs]:.2f} {pdb_phelys_dict[pdbs]:.2f} {pdb_dfgmutation_dict[pdbs]} {pdb_chainmutation_dict[pdbs]} {pdb_dbscan_dict[pdbs]} {pdb_spatial_dict[pdbs]} {pdb_dihedral_dict[pdbs]} {pdb_dihedral_dict_no_chi1[pdbs]} {pdb_phos_dict[pdbs]} {pdb_ligand_dict[pdbs]} {domain_dfgnum_dict[domain_name]} {pdb_color_dict[pdbs]} {pdb_first_resi[pdbs]} {pdb_last_resi[pdbs]} {domain_apenum_dict[domain_name]}\n')
            fhandle_databasefile.write(f'{pdbs};{pdb_reso_dict[pdbs]};{pdb_gene_dict[pdbs]};{pdb_uniprot_dict[pdbs]};{pdb_domain_dict[pdbs]};{pdb_group_dict[pdbs]};{pdb_uniacc_dict[pdbs]};{pdb_spatial_dict[pdbs]};{pdb_dihedral_dict[pdbs]};{pdb_ligand_dict[pdbs]};{pdb_gene_synonym_dict[pdbs]};{domain_pseudo_dict[domain_name]};{pdb_dfgmutation_dict[pdbs]};{pdb_chainmutation_dict[pdbs]};{pdb_phos_dict[pdbs]};{uniprot_proteinName_dict[uniprotname]};{domain_dfgnum_dict[domain_name]};{pdb_color_dict[pdbs]};{domainstart_dict[domain_name]};{domainend_dict[domain_name]};{domain_break_dict[pdbs]};{loop_break_dict[pdbs]};{domain_pseudo_dict[domain_name]};{pdb_first_resi[pdbs]};{pdb_last_resi[pdbs]};{domain_apenum_dict[domain_name]}\n')
            #fhandle_output.write(f'{pdbs} {pdb_sb_dict[pdbs]}')
            #fhandle_output.write(str(pdb_sb_dict[pdbs]))
            
        #except:
        #    print('Error in output:'+pdbs)
#Also print in output gene list help file for the website
    fhandle_output.close()
    fhandle_databasefile.close()
    statistics(pwd,pdb_domain_dict,pdb_gene_dict,pdb_group_dict,pdb_uniprot_dict,pdb_spatial_dict,pdb_dihedral_dict)
    
    cmd=('sort  -t ";" -k6,6 -k3,3 -k8,8 -k9,9 '+pwd+'/Database-input-'+today+'.tab'+' > temp')
    subprocess.call(cmd,shell=True)
    cmd=('mv temp '+pwd+'/Database-input-'+today+'.tab')
    subprocess.call(cmd,shell=True)
