#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:30:25 2020

@author: vivekmodi
"""
from Bio import PDB
import gzip

def chain_break(pwd,df):
    kinasechains_renumber_uniprot=f'{pwd}/kinasechains_renumber_uniprot'
    print('Identifying chain breaks...')
    fhandle_read=open((pwd+'/Chain_break.tab'),'r')
    fhandle_append=open((pwd+'/Chain_break.tab'),'a')
    ignoremodified=open(f'{pwd}/List_modified_aminoacid.txt','r')
    log=open(f'{pwd}/kinasepml.log','a')
    
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        df.at[i,'DomainBreak']=int(999)
        df.at[i,'LoopBreak']=int(999)
        
        fhandle_read.seek(0)
        pdb_present=0
        for lines in fhandle_read:
            lines=lines.strip();lines=lines.split()
            if pdbs in lines[0]:
                df.at[i,'DomainBreak']=int(float(lines[1]))
                df.at[i,'LoopBreak']=int(float(lines[2]))
                pdb_present=1
                break
        
        if pdb_present==0:
            try:
                handle=gzip.open((kinasechains_renumber_uniprot+'/'+pdbs+'.cif.gz'),'rt')
            except:
                log.write(f'chain_break: File not found {pdbs}.cif.gz\n')
                continue
            
            parser=PDB.MMCIFParser(QUIET=True)
            structure=parser.get_structure("PDB",handle)
            #domain_name=df.at[i,'Domain']
            domain_start=df.at[i,'DomainBegin'];domain_end=df.at[i,'DomainEnd']
            loop_start=df.at[i,'DFGnum'];loop_end=df.at[i,'APEnum']
    
            for model in structure:
                for chain in model:
                    first=1;max_diff_loop=-1;max_diff_domain=-1;break_in_loop='NO';break_in_domain='NO'
                    for residue in chain:
                        distance=999;diff_domain=-1;diff_loop=-1
                        
                        ignoremodified.seek(0)
                        if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():      #Proceed only if the residue is amino acid or modified amino acid
                            
                            if first==1:            #Get the first residue separately
                                prev_residue=residue
                                first=2
                                continue
                            if first==2:            #Second residue onwards
                                curr_residue=residue
                                for atom1 in prev_residue:
                                    if str(atom1.id)=='C':
                                        for atom2 in curr_residue:
                                            if str(atom2.id)=='N':
                                                prev_resi_c=prev_residue['C']
                                                curr_resi_n=curr_residue['N']
                                                distance=prev_resi_c - curr_resi_n
                                                #diff=int(curr_residue.id[1])-int(prev_residue.id[1])-1
                                                #if diff>max_diff_chain:
                                                #    max_diff_chain=diff
        
                                if int(prev_residue.id[1])>=int(domain_start) and int(prev_residue.id[1])<=int(domain_end) and distance>2:       #If 'C' or 'N' atom is not resolved then default distance will be 999
                                    diff_domain=int(curr_residue.id[1])-int(prev_residue.id[1])-1
                                    if diff_domain>max_diff_domain:
                                        max_diff_domain=diff_domain
                                    #if diff>max_diff_chain:
                                        #max_diff_chain=diff
                                    break_in_domain='YES'
                                    #print(f'{pdbs} {max_diff_domain} {domain_start} {domain_end}')
                                    #break_in_chain='YES'
                                if (int(prev_residue.id[1])>=int(loop_start)-1 and int(prev_residue.id[1])<=int(loop_end)+1) and distance>2:       #If 'C' or 'N' atom is not resolved then default distance will be 999
                                    diff_loop=int(curr_residue.id[1])-int(prev_residue.id[1])-1
                                    if diff_loop>max_diff_loop:
                                        max_diff_loop=diff_loop
                                    break_in_loop='YES'
                                    #print(f'{pdbs} {max_diff_loop} {loop_start} {loop_end}')
        
                                #if distance>2:
                                #    break_in_chain='YES'
                                #if int(prev_residue.id[1])>=int(domain_start) and int(prev_residue.id[1])<=int(domain_end) and distance>2:
                                #    break_in_domain='YES'
                                   # print(uni_name,domain_start,domain_end)
                                #print(filename[-8:-4],prev_residue.id[1],curr_residue.id[1],chain.id,diff)      #Print as stdout
        
                                prev_residue=curr_residue
            df.at[i,'DomainBreak']=max_diff_domain;df.at[i,'LoopBreak']=max_diff_loop
            fhandle_append.write(f'{pdbs} {df.at[i,"DomainBreak"]} {df.at[i,"LoopBreak"]}\n')
    fhandle_append.close()
    fhandle_read.close() 
    log.close()
    return (df)