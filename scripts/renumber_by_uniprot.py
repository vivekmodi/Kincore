#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:47:57 2020

@author: vivekmodi
"""
#Usage from commandline python renumber_by_uniprot.py /kinasechains /kinasesifts /kinasechains_renumber_uniprot 1xyzA
#Usage from another script renumber_by_uniprot(kinasechains,kinasesifts,kinasechains_renumbered,pdbs)

import os, gzip, subprocess
from Bio import PDB, SeqIO

def renumber_by_uniprot(pwd,df):
    kinasechains=f'{pwd}/kinasechains'
    kinasesifts=f'{pwd}/kinasesifts'
    kinasechains_renumber_uniprot=f'{pwd}/kinasechains_renumber_uniprot'
    ignoremodified=open(f'{pwd}/List_modified_aminoacid.txt','r')
    log=open(f'{pwd}/kinasepml.log','a')
    
    print('Renumbering MMCIF files by Uniprot numbering scheme...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        
        filename=(kinasechains+'/'+pdbs[0:4].upper()+pdbs[4]+".cif.gz")
        if os.path.isfile(kinasechains_renumber_uniprot+"/"+pdbs[0:5]+".pdb.gz"):
            if os.path.isfile(kinasechains_renumber_uniprot+"/"+pdbs[0:5]+".cif.gz"):
                continue 
        
        try:
            handle=gzip.open(filename,"rt")
        except:
            log.write(f"renumber_by_uniprot: file does not exist: {filename}\n")
            continue
        try:
            siftshandle=gzip.open((kinasesifts+'/'+pdbs[0:4].lower()+".csv.gz"),"rt")
        except:
            log.write(f"renumber_by_uniprot: file does not exist: {pdbs[0:4]}.lower().csv.gz\n")
            continue
        
        #if not os.path.isfile(filename):
        #    print("Error: Function renumber_pdbs: file does not exist:"+filename+"\n")
        #    return
            
        #if not os.path.isfile((kinasesifts+'/'+pdbs[0:4].lower()+".csv.gz")):
        #    print("Error: Function renumber_pdbs: file does not exist: "+pdbs[0:4].lower()+".csv.gz"+"\n")
        #    return
        
               
    
        
        parser=PDB.MMCIFParser(QUIET=True)
        
        structure=parser.get_structure("pdbs[0:4]",handle)
        
        df.at[i,'Date']=structure.header['deposition_date']
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    ignoremodified.seek(0)
                    if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():
                        resid=list(residue.id)
                        resid[1]=resid[1]+1000      #Change the residue numbers to a larger value so that while renumbering two residues do not have the same number
                        residue.id=tuple(resid)
                       
    
        for model in structure:
            for chain in model:
                for residue in list(chain):     #list() is used as a hack because after using chain.detach below the next residue is skipped, chain jumps to second residue after the current residue
                    ignoremodified.seek(0)
                    if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():
                        #residue_in_uniprot=0
                        resid=list(residue.id)
                        siftshandle.seek(0)
                        for line in siftshandle:
                            line=line.split(",")
                            residue_with_insert_code=(str(resid[1]-1000)+resid[2])
                            residue_with_insert_code=residue_with_insert_code.strip()
                            if str(residue_with_insert_code)==str(line[2]) and str(line[4])==str(chain.id) and str(line[5])==str(-9999):      #Do not print residue not present in Uniprot in PDBfile
                                #print(resid[0],resid[1],resid[2])
                                #chain.detach_child((resid[0], (resid[1]), ' '))         #Check the note in loop above; Also the argument for detach_child should be a tuple like residue.id
                                chain.detach_child(tuple(resid))
                                continue
                            
                            if 'Insertion' in str(line[8]) and str(residue_with_insert_code)==str(line[2]) and str(line[4])==str(chain.id) and str(line[5])!=str(-9999):
                                resid[1]=int(line[5][0:-1])     
                                resid[2]=line[5][-1]       #Insertion code
                                residue.id=tuple(resid)
                                residue_in_uniprot=1
                                continue
                            
                            if 'Linker' in str(line[8]) and str(residue_with_insert_code)==str(line[2]) and str(line[4])==str(chain.id) and str(line[5])!=str(-9999):
                                resid[1]=int(line[5][0:-1])     
                                resid[2]=line[5][-1]       #Insertion code
                                residue.id=tuple(resid)
                                residue_in_uniprot=1
                                continue
                            
                            if str(residue_with_insert_code)==str(line[2]) and str(line[4])==str(chain.id):    #The input from line is always str. To compare between same datatypes resid is also converted to str
                                resid[1]=int(line[5])
                                residue.id=tuple(resid)
                                residue_in_uniprot=1
                                continue
                            #if residue_in_uniprot==0:
                            #    chain.detach_child((' ', (resid[1]), ' '))
    
        siftshandle.close()
        io=PDB.MMCIFIO()
        io.set_structure(structure)
        io.save(kinasechains_renumber_uniprot+'/'+pdbs[0:4].upper()+pdbs[4]+".cif")
        cmd=('gzip -f '+kinasechains_renumber_uniprot+'/'+pdbs[0:4].upper()+pdbs[4]+'.cif')
        subprocess.call(cmd, shell=True)
        
        io=PDB.PDBIO()
        io.set_structure(structure)
        io.save(kinasechains_renumber_uniprot+'/'+pdbs[0:4].upper()+pdbs[4]+".pdb")
        cmd=('gzip -f '+kinasechains_renumber_uniprot+'/'+pdbs[0:4].upper()+pdbs[4]+'.pdb')
        subprocess.call(cmd, shell=True)
        
    log.close()
    return df