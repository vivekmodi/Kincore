#!/usr/bin/env python
# coding: utf-8




import os, sys, csv
sys.path.append(os.getcwd()+'/scripts')

import subprocess, numpy as np, pandas as pd
import openpyxl
from datetime import datetime
from Bio import *
from create_dirs import create_dirs
from Bio.Blast.Applications import NcbipsiblastCommandline
from download_pdbaa import download_pdbaa
from create_blastdb import create_blastdb
#from run_psiblast import run_psiblast
from read_psiblast import read_psiblast
from create_motifs_file import create_motifs_file
from non_human_genename import non_human_genename
from uniprotseq import uniprotseq
from create_gene_dict import gene_dict
from get_residues_from_alignment import get_residues_from_alignment
from genes_from_trembl import genes_from_trembl
from pseudokinases import identify_pseudokinases
from identify_motifs_nonhuman import identify_motifs_nonhuman
from download_cifs import download_cifs
from download_sifts import download_sifts
from split_chains import split_chains
from parse_sifts import parse_sifts
from renumber_by_uniprot import renumber_by_uniprot
from renumber_by_alignment import renumber_by_alignment
from download_phases import download_phases
from run_phoenix import run_phoenix
from get_seq_from_cif import get_seq_from_cif
from dihedrals import compute_dihedrals
from read_dihedrals import read_dihedrals
from chain_break import chain_break
from gene_synonym import gene_synonym
from identify_mutation import identify_mutation
from seqFormatHtml import format_seq_html
from seqFormatHtml import format_seq_text
from compute_all_distances import compute_all_distances
from min_distance_from_ligand import min_distance_from_ligand
from extract_ligands import extract_ligands
from assign_labels import spatial_labels, dihedral_labels, dihedral_labels_no_chi
from chainColor import chain_color
from run_edia import run_edia
from read_edia import read_edia
from classify_ligands import classify_ligands
from chelix_disposition import chelix_disposition
from createPymolSession import subListPymolSession
from geneListHelp import geneListHelp





today=str(datetime.now())[0:10].strip()
fhandle_update_date=open('update-date.txt','w')
fhandle_update_date.write(f'{today}')
fhandle_update_date.close()
pwd=os.getcwd()
print("Present working directory: "+pwd)




def run_psiblast(pwd):
    psiblast_cline=NcbipsiblastCommandline(db=f'{pwd}/pdbaa_psiblast_dir/pdbaa.db',in_pssm=f'{pwd}/pdbaa_psiblast_dir/AurkaPsiblastIter6PSSM.asn', out=f'{pwd}/pdbaa_psiblast_dir/AURKA.pdbaa.xml',num_alignments=10000,num_descriptions=10000,outfmt=5)
    psiblast_cline()    





#Main
df=pd.DataFrame()
#create_dirs(pwd)
#download_pdbaa(pwd+'/pdbaa_psiblast_dir')
#create_blastdb('pdbaa',f'{pwd}/pdbaa_psiblast_dir')
#run_psiblast(pwd)

psiblast_result=(f'{pwd}/pdbaa_psiblast_dir/AURKA.pdbaa.xml');
excluded_output=(f'{pwd}/pdbaa_psiblast_dir/psiblast_excluded.log')
df=read_psiblast(psiblast_result, excluded_output,df)   #sequences from pdbaa also contain cloning tags

#create_motifs_file(pwd)
df=gene_dict(pwd,df)
#df=get_residues_from_alignment(pwd,df)
#df=non_human_genename(pwd,df)
df=uniprotseq(pwd,df)    
#df=genes_from_trembl(pwd,df)      #get the cases which are unreviewed; only after running above three functions
#df=identify_motifs_nonhuman(pwd,df)    #align with human HMM and identify the motifs
df=identify_pseudokinases(df)  

#download_cifs(f'{pwd}/kinasecifs',df)
#download_sifts(f'{pwd}/kinasesifts',df)
split_chains(pwd,df)
parse_sifts(f'{pwd}/kinasesifts',df)
renumber_by_uniprot(pwd,df)
renumber_by_alignment(pwd,df)

#download_phases(pwd,df)
#run_phoenix(pwd,df)
#run_edia(pwd,df)
df=read_edia(pwd,df)
compute_dihedrals(pwd,df)
df=read_dihedrals(pwd,df)
df=chain_break(pwd,df)
df=gene_synonym(pwd,df)        #modify to include non-human genes in the list
df=identify_mutation(pwd,df)
df=extract_ligands(pwd,df)
format_seq_html(pwd,df)
df=get_seq_from_cif(pwd,df)      #Always use after generating uniprot numbered files - this uses incorrect chain id, maybe try .pdb files


df=compute_all_distances(pwd,df)

df=chelix_disposition(pwd,df)

df=spatial_labels(pwd,df)  #think of excluding structures with DFGmutations
df=dihedral_labels(pwd,df)
df=dihedral_labels_no_chi(pwd,df)
df=chain_color(df)

df=classify_ligands(pwd,df)
subListPymolSession(pwd,df)   #This function also copies coordinate files
geneListHelp(pwd,df)
df_sorted=df.sort_values(['Specie','Group','Gene','Spatial','Dihedral','Ligand_label']).copy()
df_sorted.to_excel(f'Kinases_df-{today}.xlsx',index=False)   #Write excel and csv
df_sorted.to_csv(f'Kinases_df-{today}.csv',sep='\t',index=False)




