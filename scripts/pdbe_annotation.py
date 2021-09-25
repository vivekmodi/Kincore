#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:40:00 2020

@author: vivekmodi
"""
import json, sys, subprocess, os
from datetime import datetime
import pandas as pd
from Bio.SeqUtils import seq1

def create_json_dirs(pwd,pdb):
    dir_name=pdb[1:3].lower()
    if not os.path.isdir(f'{pwd}/JSON/{dir_name}'):
        cmd=f'mkdir {pwd}/JSON/{dir_name}'
        subprocess.call(cmd,shell=True)
    return dir_name

def create_json(pwd,filename):
    df=pd.read_csv(f'{pwd}/{filename}',sep='\t',header='infer')
    pdb_list=list();pdb_skip=list();group=dict();gene=dict();release_date=dict()
    year=str(datetime.now())[0:4];month=str(datetime.now())[5:7];day=str(datetime.now())[10:13].strip()
    today=f'{day}/{month}/{year}'
    for i in df.index:
        pdb=df.at[i,'PDBid'][0:4]

        if str(df.at[i,'Author_Aspnum']).lower() == 'nan':   #Do not include the structures where Phe is not resolved; Pandas read null as nan
            pdb_skip.append(pdb)
            continue

        pdb_list.append(pdb)
        release_date[pdb]=today #Annotation release date
        group[pdb]=df.at[i,'Group']
        gene[pdb]=df.at[i,'Gene']
    pdb_list=set(pdb_list)

    for pdb in pdb_list:
        if pdb in pdb_skip:
            continue

        pdbe=dict()
        pdbe["data_resource"]="Kincore"
        pdbe["resource_version"]="1.0.0"
        pdbe["resource_entry_url"]=f"http://dunbrack.fccc.edu/kincore/PDB/{pdb}"
        pdbe["model_coordinates_url"]=f"http://dunbrack.fccc.edu/kincore/static/downloads/coordinateFiles/{group[pdb]}_{gene[pdb]}_{pdb}_uniNum.zip"
        pdbe["release_date"]=f"{release_date[pdb]}"
        pdbe["pdb_id"]=f"{pdb}".lower()
        pdbe["chains"]=list()

        pdbe["evidence_code_ontology"]=[{"eco_term": "computational combinatorial evidence used in automatic assertion", "eco_code": "ECO_0000246"}]

        pdbe["sites"]=list()

        chain_label=dict()
        site_id=0

        for i in df.index:

            if pdb in df.at[i,"PDBid"]:
                site_id=site_id+1
                chain_label=df.at[i,'PDBid'][4:]
                spatial=df.at[i,'Spatial']
                dihedral=df.at[i,'Dihedral']
                dfgnum=str(df.at[i,'Author_Aspnum'])      #Asp num is used
                if '.' in dfgnum:
                    dfgnum=dfgnum[0:-2]      #Remove the trailing decimal if present

                aatype=(df.at[i,"Author_Aspres"]).upper()

                residues=[{"pdb_res_label": dfgnum,"aa_type": aatype,"site_data": [{"site_id_ref": site_id,"confidence_classification": "curated"}]}]
                sites={"site_id": site_id, "label": spatial, "additional_site_annotations": {"dihedral_label": dihedral}}

                pdbe["chains"].append({"chain_label": chain_label,"residues": residues})
                pdbe["sites"].append(sites)

        dir_name=create_json_dirs(pwd,pdb)
        fhandle_json=open(f'JSON/{dir_name}/'+pdb.lower()+'.json','w')
        json.dump(pdbe,fhandle_json,indent=2)
        fhandle_json.close()


if __name__=='__main__':
    filename=sys.argv[1]
    pwd='/home/vivek/Applications/Flask/Kincore'     #Location in workhorse
    create_json(pwd,filename)
