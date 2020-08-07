#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:53:35 2020

@author: vivekmodi
"""
import gzip, os, subprocess
from Bio import PDB

def renumber_by_alignment (pwd,df):
    kinasechains_renumber_uniprot=f'{pwd}/kinasechains_renumber_uniprot'
    kinasechains_renumber_alignment=f'{pwd}/kinasechains_renumber_alignment'
    print('Renumbering MMCIF files by alignment column numbers...')
    ignoremodified=open(f'{pwd}/List_modified_aminoacid.txt','r')
    fhandle_column=open((pwd+'/All-organisms-alignment-residue-corresspondence.tab'),'r')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        
        #if 'HUMAN' not in df.at[i,'UniprotID']:    #For now skip non-human
        #    continue
        
        if os.path.isfile(kinasechains_renumber_alignment+"/"+pdbs[0:5]+".pdb.gz"):
            if os.path.isfile(kinasechains_renumber_alignment+"/"+pdbs[0:5]+".cif.gz"):
                continue
            
        fhandle_column.seek(0)
        res_aln=dict()
        for line in fhandle_column:
            line=line.strip();line=line.split()
            #residue_with_insert_code=(str(resid[1]-3000)+resid[2])
            #residue_with_insert_code=residue_with_insert_code.strip()
            if line[4]==df.at[i,'Domain'] and line[1]==df.at[i,'UniprotID']:     #match domain
                res_aln[int(line[6])]=int(line[7])          #res[uniprot]-->column_number -- this dictionary will be used below to assign new residue numbers)


        #print("Renumbering PDBs by alignment column numbers:"+pdbs)
    #for name in pdb_chainlist:
    #    if str(name[0:4]).lower()=="5ezv" or str(name[0:4]).lower()=="4wb7" or str(name[0:4]).lower()=="6byr" or str(name[0:4]).lower()=="3bea" or \
    #    str(name[0:4]).lower()=="3krl" or str(name[0:4]).lower()=="3krj" or str(name[0:4]).lower()=="3dpk":
    #        continue

     #   print(uniprot[name[0:5]],name[0:5])

        pdbfilename=(kinasechains_renumber_uniprot+'/'+pdbs[0:5]+'.cif.gz')
        if not os.path.isfile(pdbfilename):
            print("Error: Function renumber_by_alignement: file does not exist:"+pdbfilename+"\n")
            return
            
        parser=PDB.MMCIFParser(QUIET=True)
        handle=gzip.open(pdbfilename,'rt')
        structure=parser.get_structure("pdbs[0:4]",handle)

        for model in structure:
            for chain in model:
                for residue in chain:
                    ignoremodified.seek(0)
                    if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():
                        resid=list(residue.id)
                        resid[1]=resid[1]+3000      #Change the residue numbers to a larger value so that while renumbering two residues do not have the same number
                        residue.id=tuple(resid)

        for model in structure:
            for chain in model:
                for residue in list(chain):     #list() is used as a hack because after using chain.detach below the next residue is skipped, chain jumps to second residue after the current residue
                    ignoremodified.seek(0)
                    if residue.id[0]==' ' or (residue.id[0][2:]+'\n') in ignoremodified.readlines():
                        present_in_alignment=0
                        resid=list(residue.id)
                        #fhandle_column.seek(0)
                        #for line in fhandle_column:
                        #    line=line.strip();line=line.split()
                            #residue_with_insert_code=(str(resid[1]-3000)+resid[2])
                            #residue_with_insert_code=residue_with_insert_code.strip()
                        check_key=int(resid[1]-3000)
                        if  check_key in dict.keys(res_aln):

                            resid[1]=res_aln[resid[1]-3000]
                            #if line[0]==domain:
                            #    if str(resid[1]-3000)==str(line[2]):    #The input from line is always str. To compare between same datatypes resid is also converted to str
                            #        resid[1]=int(line[3])
                            residue.id=tuple(resid)
                            #print(pdbs,residue.id)
                            present_in_alignment=1
                            #continue
                        if present_in_alignment==0:
                            #chain.detach_child((resid[0], (resid[1]), ' '))
                            chain.detach_child(tuple(resid))

        #fhandle_column.close()
        io=PDB.PDBIO()
        io.set_structure(structure)
        #if not os.path.isdir((kinasechains_renumber_alignment+'/'+pdbs[0:5]+'.pdb.gz')"./uniprot/"+uniprot[name[0:5]]):
        #    cmd=("mkdir uniprot/"+uniprot[name[0:5]])
        #    subprocess.call(cmd,shell=True)
        io.save(kinasechains_renumber_alignment+"/"+pdbs[0:5]+".pdb")
        #cmd=(f'cp {kinasechains_renumber_alignment}/{pdbs[0:4]}.upper(){pdbs[4]}.pdb {pwd}/static/downloads/alignment-numbered')
        #subprocess.call(cmd,shell=True)
        cmd=('gzip -f '+kinasechains_renumber_alignment+'/'+pdbs[0:5]+'.pdb')
        subprocess.call(cmd, shell=True)

        io=PDB.MMCIFIO()
        io.set_structure(structure)
        io.save(kinasechains_renumber_alignment+"/"+pdbs[0:5]+".cif")
        #cmd=(f'cp {kinasechains_renumber_alignment}/{pdbs[0:4]}.upper(){pdbs[4]}.cif {pwd}/static/downloads/alignment-numbered')
        #subprocess.call(cmd,shell=True)
        cmd=('gzip -f '+kinasechains_renumber_alignment+'/'+pdbs[0:5]+'.cif')
        subprocess.call(cmd, shell=True)
        #else:
        #    io.save("./uniprot/"+uniprot[name[0:5]]+"/"+name[0:5]+".column.pdb")
    return
