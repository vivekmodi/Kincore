#!/usr/bin/env python
# coding: utf-8

import os, sys, subprocess
sys.path.append(os.getcwd()+'/scripts')
sys.path.append(os.getcwd()+'/funpdbe-validator')

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
from update_database import update_database
from pdbe_annotation import create_json
from validate import validate
from transfer_to_dunbrack3 import transfer_to_dunbrack3
from transfer_to_pdbe import transfer_to_pdbe
from create_fasta_with_labels import create_fasta_with_labels

AADICT={'GLY':'G','ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y',\
        'TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','ARG':'R','HIS':'H','LYS':'K',\
        'ASP':'D','GLU':'E','CYS':'C','PRO':'P','SEC':'U','TPO':'T','CME':'C','CSS':'C',\
        'MSE':'M','OCY':'C','PTR':'Y','SEP':'S','CAF':'C','LGY':'K','CAS':'C','CSO':'C','CSX':'C',\
        'MK8':'E','NEP':'H','NMM':'R','CSD':'C','CYO':'Y','OCS':'C','OCY':'C','SCS':'C','ALY':'A',\
        'KCX':'K','MHO':'M','T8L':'T','CY0':'C','UNK':'X','YTH':'T'}

def identify_working_direct():
    process=subprocess.Popen('uname -a',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    stdout,stderr=process.communicate()
    if 'vivek-XPS' in str(stdout):
        pwd='/home/vivekmodi/Applications/Flask/Kinases'  #location in laptop
    else:
        pwd='/home/vivek/Applications/Flask/Kincore'     #location in workhorse;required to start cronjob
    return pwd

def Main(pwd):
    log=open(f'{pwd}/kinasepml.log','w')
    today=str(datetime.now())[0:10].strip()
    df=pd.DataFrame()
    #create_dirs(pwd)          #No need to run this because directories already exist
    #download_pdbaa(pwd+'/pdbaa_psiblast_dir')
    create_blastdb('pdbaa',f'{pwd}/pdbaa_psiblast_dir')
    run_psiblast(f'{pwd}/pdbaa_psiblast_dir','pdbaa','AurkaPsiblastIter6PSSM.asn','AURKA.pdbaa.xml')

    df=read_psiblast(pwd,df,f'{pwd}/pdbaa_psiblast_dir/AURKA.pdbaa.xml', f'{pwd}/pdbaa_psiblast_dir/psiblast_excluded.log')   #sequences from pdbaa also contain cloning tags
    download_cifs(f'{pwd}/kinasecifs',df)
    get_release_date(f'{pwd}/kinasecifs',df)
    create_motifs_file(pwd)    # This function also prints a file Not_found_in-alignment.txt, which has the Uniprots in with conserved residues are missing.
    df=gene_dict(pwd,df)        # Also prints New_uniprots.txt
    df=uniprotseq(pwd,df)

    df=identify_pseudokinases(df)

    
    
    download_sifts(f'{pwd}/kinasesifts',df)
    try:
        split_chains(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        parse_sifts(f'{pwd}/kinasesifts',df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=renumber_by_uniprot(pwd,df)      #Also returns deposition date
    except Exception as e:
        log.write(f'{e}\n')
    try:
        renumber_by_alignment(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=identify_author_dfg(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    download_phases(pwd,df)
    try:
        run_phoenix(pwd,df)
        run_edia(pwd,df)
        df=read_edia(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        compute_dihedrals(pwd,df)
        df=read_dihedrals(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=chain_break(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    #df=gene_synonym(pwd,df)        #modify to include non-human genes in the list
    df=identify_mutation(pwd,df,AADICT)
    try:
        df=extract_ligands(pwd,df)
    except Exception as e:
       log.write(f'{e}\n')
    try:
        format_seq_html(pwd,df,AADICT)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=get_seq_from_cif(pwd,df)      #Do not need it right now; Always use after generating uniprot numbered files - this uses incorrect chain id, maybe try .pdb files
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=compute_all_distances(pwd,df)
        df=chelix_disposition(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        df=spatial_labels(pwd,df)
        df=dihedral_labels(df,0.45)
    except Exception as e:
        log.write(f'{e}\n')
        
    df=chain_color(df)
    try:
        df=classify_ligands(pwd,df)
    except Exception as e:
        log.write(f'{e}\n')
        
    geneListHelp(pwd,df)
    copy_ngl_files(pwd,df)
    create_datefile(pwd)
    df_sorted=df.sort_values(['Specie','Group','Gene','Spatial','Dihedral','Ligand_label']).copy()
    df_sorted=create_fasta_with_labels(pwd,df_sorted)    #Prints sorted fasta file with sequences from pdbaa
    df_sorted.to_excel(f'Kinases_df-{today}.xlsx',index=False)   #Write excel and csv
    df_sorted.to_csv(f'Kinases_df-{today}.csv',sep='\t',index=False)

    df_datesorted=df.sort_values(['Date'],ascending=False).copy()
    df_datesorted[['Date','Specie','UniprotID','Gene','PDBid']].to_excel(f'Kinases_df_date_sorted.xlsx',index=False)

    try:
        update_database(df)
    except Exception as e:
        log.write(f'{e}\n')
    try:
        subListPymolSession(pwd,df)   #This function also copies coordinate files
    except Exception as e:
        log.write(f'{e}\n')

    try:
        create_json(pwd,f'Kinases_df-{today}.csv')
    except Exception as e:
        log.write(f'{e}\n')
    try:
        validate(pwd,f'Kinases_df-{today}.csv')
    except Exception as e:
        log.write(f'{e}\n')
    #try:
    #    transfer_to_dunbrack3(pwd)
    #except Exception as e:
    #    log.write(f'{e}\n')
    #try:
    #    transfer_to_pdbe(pwd)
    #except Exception as e:
    #    log.write(f'{e}\n')
    
    log.close()
    


if __name__ == '__main__':
    pwd=identify_working_direct()
    
    print("Present working directory: "+pwd)
    Main(pwd)
