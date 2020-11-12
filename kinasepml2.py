#!/usr/bin/env python
# coding: utf-8

import os, sys
sys.path.append(os.getcwd()+'/scripts')

import pandas as pd
from datetime import datetime
from create_dirs import create_dirs
from download_pdbaa import download_pdbaa
from create_blastdb import create_blastdb
from run_psiblast import run_psiblast
from read_psiblast import read_psiblast
from create_motifs_file import create_motifs_file
from uniprotseq import uniprotseq
from create_gene_dict import gene_dict
from pseudokinases import identify_pseudokinases
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
from identify_mutation import identify_mutation
from seqFormatHtml import format_seq_html
from compute_all_distances import compute_all_distances
from extract_ligands import extract_ligands
from assign_labels import spatial_labels, dihedral_labels
from chainColor import chain_color
from run_edia import run_edia
from read_edia import read_edia
from classify_ligands import classify_ligands
from chelix_disposition import chelix_disposition
from createPymolSession import subListPymolSession
from geneListHelp import geneListHelp
from copy_ngl_files import copy_ngl_files
from create_datefile import create_datefile
from get_release_date import get_release_date
from identify_author_dfg import identify_author_dfg
from pdbe_annotation import create_json




def Main(pwd):
    today=str(datetime.now())[0:10].strip()
    df=pd.DataFrame()
    create_dirs(pwd)          #No need to run this because directories already exist
    download_pdbaa(pwd+'/pdbaa_psiblast_dir')
    create_blastdb('pdbaa',f'{pwd}/pdbaa_psiblast_dir')
    run_psiblast(f'{pwd}/pdbaa_psiblast_dir','pdbaa','AurkaPsiblastIter6PSSM.asn','AURKA.pdbaa.xml')

    df=read_psiblast(df,f'{pwd}/pdbaa_psiblast_dir/AURKA.pdbaa.xml', f'{pwd}/pdbaa_psiblast_dir/psiblast_excluded.log')   #sequences from pdbaa also contain cloning tags

    create_motifs_file(pwd)    # This function also prints a file Not_found_in-alignment.txt, which has the Uniprots in with conserved residues are missing.
    df=gene_dict(pwd,df)        # Also prints New_uniprots.txt
    df=uniprotseq(pwd,df)

    df=identify_pseudokinases(df)

    download_cifs(f'{pwd}/kinasecifs',df)
    get_release_date(f'{pwd}/kinasecifs',df)
    download_sifts(f'{pwd}/kinasesifts',df)
    split_chains(pwd,df)
    parse_sifts(f'{pwd}/kinasesifts',df)
    df=renumber_by_uniprot(pwd,df)      #Also returns deposition date
    renumber_by_alignment(pwd,df)
    df=identify_author_dfg(pwd,df)
    download_phases(pwd,df)
    run_phoenix(pwd,df)
    run_edia(pwd,df)
    df=read_edia(pwd,df)
    compute_dihedrals(pwd,df)
    df=read_dihedrals(pwd,df)
    df=chain_break(pwd,df)
    #df=gene_synonym(pwd,df)        #modify to include non-human genes in the list
    df=identify_mutation(pwd,df)
    df=extract_ligands(pwd,df)
    format_seq_html(pwd,df)
    df=get_seq_from_cif(pwd,df)      #Do not need it right now; Always use after generating uniprot numbered files - this uses incorrect chain id, maybe try .pdb files


    df=compute_all_distances(pwd,df)
    df=chelix_disposition(pwd,df)
    df=spatial_labels(pwd,df)
    df=dihedral_labels(df,0.45)
    df=chain_color(df)

    df=classify_ligands(pwd,df)
    geneListHelp(pwd,df)
    copy_ngl_files(pwd,df)
    create_datefile(pwd)
    df_sorted=df.sort_values(['Specie','Group','Gene','Spatial','Dihedral','Ligand_label']).copy()
    df_sorted.to_excel(f'Kinases_df-{today}.xlsx',index=False)   #Write excel and csv
    df_sorted.to_csv(f'Kinases_df-{today}.csv',sep='\t',index=False)

    df_datesorted=df.sort_values(['Date'],ascending=False).copy()
    df_datesorted[['Date','Specie','UniprotID','Gene','PDBid']].to_excel(f'Kinases_df_date_sorted.xlsx',index=False)

    subListPymolSession(pwd,df)   #This function also copies coordinate files

    create_json(pwd,f'Kinases_df-{today}.csv')





if __name__ == '__main__':
    pwd=os.getcwd()
    print("Present working directory: "+pwd)
    Main(pwd)
