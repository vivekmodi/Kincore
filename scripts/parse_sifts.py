#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:03:58 2020

@author: vivekmodi
"""


import os, gzip, subprocess
import xml.etree.ElementTree as ET

def parse_sifts(kinasesifts,df):
    print(f'Parsing Sifts xml files...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']

        filename=(kinasesifts+"/"+pdbs[0:4].lower()+".xml.gz")
        filecsv=(kinasesifts+"/"+pdbs[0:4].lower()+".csv.gz")

        prev_uninum_insertion=-9999
        prev_uninum_linker=-9999
        insertion_code='A'
        if not os.path.isfile(filename):
            print("Error: Function parse_sifts: file does not exist:"+filename+"\n")
            return

        #if os.path.isfile(filecsv):
        #    continue

        handle=gzip.open(filename,"rt")
        tree=ET.parse(handle)
        root=tree.getroot()

        fhandle=open((kinasesifts+"/"+pdbs[0:4].lower()+".csv"),"w")
        fhandle.write("PDBe_resnum,PDBe_resname,PDB_resnum,PDBres_name,PDB_Chain,Uniprot_resnum,Uniprot_resname,Uniprot_accession,Annotation,SecStr"+"\n")
        base='{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}'

        for entity in root:
            if entity.tag==(str(base)+'entity'):
                prev_uninum_insertion=-9999          #at the beginning of new chain set prev_uninum to default again
                prev_uninum_linker=-9999
                for segment in entity:
                    if segment.tag==(str(base)+'segment'):
                        for listResidue in segment:
                            if listResidue.tag==(str(base)+'listResidue'):
                                for residue in listResidue:
                                    pdbe_resnum=-9999;pdbe_resname='XXX';pdb_resnum=-9999;pdb_resname='XXX';pdb_chain='XXX';uni_resnum=-9999;uni_resname='XXX';uni_acc='XXX';\
                                    pdbe_property=list();pdbe_secstr='None'
                                    if residue.tag==(str(base)+'residue'):
                                        if residue.get('dbSource')=='PDBe':
                                            pdbe_resnum=residue.get('dbResNum');pdbe_resname=residue.get('dbResName')
                                            for crossRefDb in residue:

                                                if crossRefDb.get('dbSource')=="PDB":
                                                    pdb_resnum=crossRefDb.get('dbResNum');pdb_resname=crossRefDb.get('dbResName');pdb_chain=crossRefDb.get('dbChainId')


                                                if crossRefDb.get('dbSource')=="UniProt":
                                                    uni_resnum=crossRefDb.get('dbResNum');uni_resname=crossRefDb.get('dbResName');uni_acc=crossRefDb.get('dbAccessionId')
                                            for residueDetail in residue:
                                                if residueDetail.get('property')=='Annotation':
                                                    pdbe_property.append(residueDetail.text)
                                                if residueDetail.get('property')=='nameSecondaryStructure':
                                                    pdbe_secstr=residueDetail.text

                                    if 'Insertion' in pdbe_property and prev_uninum_insertion!=-9999:        #check for insertion and make sure it is not in the beginning of the chain by -9999
                                        fhandle.write(str(pdbe_resnum)+","+str(pdbe_resname)+","+str(pdb_resnum)+","+str(pdb_resname)+","+str(pdb_chain)+","+\
                                                      str(prev_uninum_insertion)+insertion_code+","+str(uni_resname)+","+str(uni_acc)+","+str(' '.join(pdbe_property))+","+str(pdbe_secstr)+"\n")
                                        insertion_code=chr(ord(insertion_code)+1)

                                    elif 'Linker' in pdbe_property and prev_uninum_linker!=-9999:        #check for insertion and make sure it is not in the beginning of the chain by -9999
                                        fhandle.write(str(pdbe_resnum)+","+str(pdbe_resname)+","+str(pdb_resnum)+","+str(pdb_resname)+","+str(pdb_chain)+","+\
                                                      str(prev_uninum_linker)+insertion_code+","+str(uni_resname)+","+str(uni_acc)+","+str(' '.join(pdbe_property))+","+str(pdbe_secstr)+"\n")
                                        insertion_code=chr(ord(insertion_code)+1)


                                    else:
                                        fhandle.write(str(pdbe_resnum)+","+str(pdbe_resname)+","+str(pdb_resnum)+","+str(pdb_resname)+","+str(pdb_chain)+","+\
                                                      str(uni_resnum)+","+str(uni_resname)+","+str(uni_acc)+","+str(' '.join(pdbe_property))+","+str(pdbe_secstr)+"\n")
                                        prev_uninum_insertion=uni_resnum
                                        prev_uninum_linker=uni_resnum
                                        insertion_code='A'
        fhandle.close()
        cmd=("gzip -f "+kinasesifts+"/"+pdbs[0:4].lower()+".csv")
        subprocess.call(cmd, shell=True)
