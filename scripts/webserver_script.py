#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:35:32 2019

@author: vivekmodi
"""

##**Make a list of functions which are common in webserver and kinasepml and write them in such a way that they could be used for both the scripts

import subprocess, sys, os 
from Bio import SeqIO, SearchIO, PDB
import dihedrals_for_webserver, chelix, dfg_conf, dfg_bkbone_conf
import numpy as np

def identify_state(pwd,pdbfilename):
    group=dict();xdfg=dict();dfg_asp=dict();dfg_phe=dict();xdfg_res=dict();dfg_asp_res=dict();dfg_phe_res=dict();xdfg_phi=dict();xdfg_psi=dict();dfg_asp_phi=dict();
    dfg_asp_psi=dict();dfg_phe_phi=dict();dfg_phe_psi=dict();dfg_phe_chi1=dict();chelix_conf=dict();dfg_label=dict();dfg_bkbone=dict();chain_list=list()
    if '.cif' in pdbfilename.lower():
        parser=PDB.MMCIFParser()
    if '.pdb' in pdbfilename.lower():
        parser=PDB.PDBParser()
    structure=parser.get_structure(pdbfilename, f'{pwd}/server/uploads/{pdbfilename}')
    model_num=0
    for model in structure:
        if model_num>0:      #Uses just the first model right now
            continue
        model_num=1
        for chain in model:
            if len(chain.get_list())<=25:
                continue
            chain_list.append(chain.id)
            first_res=extract_seq(pwd,pdbfilename,chain.id)
            run_hmmsearch(pwd,pdbfilename,chain.id)
            group[chain.id]=identify_group(pwd,pdbfilename,chain.id)
            if group[chain.id]=='NA':
                xdfg[chain.id]=999;dfg_asp[chain.id]=999;dfg_phe[chain.id]=999;xdfg_res[chain.id]=999;dfg_asp_res[chain.id]=999;dfg_phe_res[chain.id]=999
                xdfg_phi[chain.id]=999;xdfg_psi[chain.id]=999;dfg_asp_phi[chain.id]=999;dfg_asp_psi[chain.id]=999;dfg_phe_phi[chain.id]=999
                dfg_phe_psi[chain.id]=999;dfg_phe_chi1[chain.id]=999;chelix_conf[chain.id]=999;dfg_label[chain.id]='None';dfg_bkbone[chain.id]='None'
            else:
                (alk_lys,rre_glu,rre4,hrd_his,xdfg[chain.id],dfg_asp[chain.id],dfg_phe[chain.id],alk_lys_res,rre_glu_res,rre4_res,\
                 hrd_his_res,xdfg_res[chain.id],dfg_asp_res[chain.id],dfg_phe_res[chain.id])=identify_residues(pwd,first_res,group[chain.id],pdbfilename,chain.id)
                #print(f'{rre4}, {dfg_phe}, {rre4_res}, {dfg_phe_res}')
                dis_phe_rre4=compute_distance(pwd,pdbfilename,chain.id,rre4,dfg_phe[chain.id],rre4_res,dfg_phe_res[chain.id],setting='rre4-phe')
                dis_phe_lys=compute_distance(pwd,pdbfilename,chain.id,alk_lys,dfg_phe[chain.id],alk_lys_res,dfg_phe_res[chain.id],setting='lys-phe')
                dis_sb=compute_distance(pwd,pdbfilename,chain.id,alk_lys,rre_glu,alk_lys_res,rre_glu_res,setting='saltbridge')
                (xdfg_phi[chain.id],xdfg_psi[chain.id],dfg_asp_phi[chain.id],dfg_asp_psi[chain.id],dfg_phe_phi[chain.id],dfg_phe_psi[chain.id],\
                 dfg_phe_chi1[chain.id])=dihedrals_for_webserver.compute_dihedrals(pwd,pdbfilename,chain.id,xdfg[chain.id],dfg_asp[chain.id],dfg_phe[chain.id])  #do not use the same script which is used in kinasepml, modification required here
                chelix_conf[chain.id]=chelix.chelix_conformation(dis_sb)
                dfg_label[chain.id]=dfg_conf.spatial_label(dis_phe_lys,dis_phe_rre4)
                dfg_bkbone[chain.id]=dfg_bkbone_conf.dihedral_labels(pwd,dfg_label[chain.id],xdfg_phi[chain.id],xdfg_psi[chain.id],dfg_asp_phi[chain.id],\
                                                                     dfg_asp_psi[chain.id],dfg_phe_phi[chain.id],dfg_phe_psi[chain.id],dfg_phe_chi1[chain.id])
                #print('Residue ',group,pdbfilename,first_res,rre4,alk_lys,hrd_his,dfg_asp+1,dis_phe_rre4,dis_phe_lys)
    return (group,chain_list,xdfg,dfg_asp,dfg_phe,xdfg_res,dfg_asp_res,dfg_phe_res,xdfg_phi,xdfg_psi,dfg_asp_phi,dfg_asp_psi,dfg_phe_phi,dfg_phe_psi,dfg_phe_chi1,chelix_conf,dfg_label,dfg_bkbone)

def extract_seq(pwd,pdbfilename,chain_id):
    fhandle_outputseq=open(f'{pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta','w')
    fhandle_outputseq.write('>'+pdbfilename[0:-4]+chain_id+'\n')
    if '.cif' in pdbfilename.lower():
        for record in SeqIO.parse(f'{pwd}/server/uploads/{pdbfilename}','cif-atom'):      #Extracts sequence from ATOM records
            if record.id.split(':')[1]==chain_id:
                fhandle_outputseq.write(f'{record.seq}')
                first_res=record.annotations['start']               #Extract the residue number of firt residue as this information is lost in fasta file and HMM output
    if '.pdb' in pdbfilename.lower():
        for record in SeqIO.parse(f'{pwd}/server/uploads/{pdbfilename}','pdb-atom'):      #Extracts sequence from ATOM records
            if record.id.split(':')[1]==chain_id:
                fhandle_outputseq.write(f'{record.seq}')
                first_res=record.annotations['start']               #Extract the residue number of firt residue as this information is lost in fasta file and HMM output
    fhandle_outputseq.close()
    return first_res                      #MODIFY the code so that if there are multiple chains in the PDB then they are printed in separate fasta files

def run_hmmsearch(pwd,pdbfilename,chain_id):
    #hmmfilename=os.path.join('/home/vivekmodi/Applications/Flask/Kinases/server/uploads/'+'Human-PK.hmm')
    cmd=(f'hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_AGC.hmmer.txt {pwd}/server/AGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_CAMK.hmmer.txt {pwd}/server/CAMK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_CK1.hmmer.txt {pwd}/server/CK1.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_CMGC.hmmer.txt {pwd}/server/CMGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_NEK.hmmer.txt {pwd}/server/NEK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_RGC.hmmer.txt {pwd}/server/RGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_STE.hmmer.txt {pwd}/server/STE.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_TKL.hmmer.txt {pwd}/server/TKL.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_TYR.hmmer.txt {pwd}/server/TYR.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_HASP.hmmer.txt {pwd}/server/HASP.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_WNK.hmmer.txt {pwd}/server/WNK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_BUB.hmmer.txt {pwd}/server/BUB.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;\
         hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}{chain_id}_ULK.hmmer.txt {pwd}/server/ULK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}{chain_id}.fasta;')
    subprocess.call(cmd,shell=True)
    return

def identify_group(pwd,pdbfilename,chain_id):
        hmm_result_AGC=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_AGC.hmmer.txt',format='hmmer3-text')
        hmm_result_CAMK=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_CAMK.hmmer.txt',format='hmmer3-text')
        hmm_result_CK1=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_CK1.hmmer.txt',format='hmmer3-text')
        hmm_result_CMGC=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_CMGC.hmmer.txt',format='hmmer3-text')
        hmm_result_NEK=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_NEK.hmmer.txt',format='hmmer3-text')
        hmm_result_RGC=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_RGC.hmmer.txt',format='hmmer3-text')
        hmm_result_STE=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_STE.hmmer.txt',format='hmmer3-text')
        hmm_result_TKL=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_TKL.hmmer.txt',format='hmmer3-text')
        hmm_result_TYR=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_TYR.hmmer.txt',format='hmmer3-text')
        hmm_result_HASP=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_HASP.hmmer.txt',format='hmmer3-text')
        hmm_result_WNK=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_WNK.hmmer.txt',format='hmmer3-text')
        hmm_result_BUB=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_BUB.hmmer.txt',format='hmmer3-text')
        hmm_result_ULK=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_ULK.hmmer.txt',format='hmmer3-text')

        eval=dict()
        eval['AGC']=100;eval['CAMK']=100;eval['CK1']=100;eval['CMGC']=100;eval['NEK']=100;eval['RGC']=100;eval['STE']=100;eval['TKL']=100;eval['TYR']=100;
        eval['HASP']=100;eval['WNK']=100;eval['BUB']=100;eval['ULK']=100;
        for hits in hmm_result_AGC:
            for hsp in hits:
                eval['AGC']=hsp.evalue
        for hits in hmm_result_CAMK:
            for hsp in hits:
                eval['CAMK']=hsp.evalue
        for hits in hmm_result_CK1:
            for hsp in hits:
                eval['CK1']=hsp.evalue
        for hits in hmm_result_CMGC:
            for hsp in hits:
                eval['CMGC']=hsp.evalue
        for hits in hmm_result_NEK:
            for hsp in hits:
                eval['NEK']=hsp.evalue
        for hits in hmm_result_RGC:
            for hsp in hits:
                eval['RGC']=hsp.evalue
        for hits in hmm_result_STE:
            for hsp in hits:
                eval['STE']=hsp.evalue
        for hits in hmm_result_TKL:
            for hsp in hits:
                eval['TKL']=hsp.evalue
        for hits in hmm_result_TYR:
            for hsp in hits:
                eval['TYR']=hsp.evalue
        for hits in hmm_result_HASP:
            for hsp in hits:
                eval['HASP']=hsp.evalue
        for hits in hmm_result_WNK:
            for hsp in hits:
                eval['WNK']=hsp.evalue
        for hits in hmm_result_BUB:
            for hsp in hits:
                eval['BUB']=hsp.evalue
        for hits in hmm_result_ULK:
            for hsp in hits:
                eval['ULK']=hsp.evalue

        minEval=100;groupName='NA'
        for groups in ('AGC','CAMK','CK1','CMGC','NEK','RGC','STE','TKL','TYR','HASP','WNK','BUB','ULK'):
            if eval[groups]<minEval:
                minEval=eval[groups]
                groupName=groups

        return groupName

def identify_residues(pwd,first_res,group,pdbfilename,chain_id):
    lys=dict();glu=dict();glu4=dict();his=dict();xdfg=dict();asp=dict();phe=dict()
    lys={'AGC':30,'CAMK':30,'CK1':30,'CMGC':30,'NEK':30,'RGC':29,'STE':30,'TKL':28,'TYR':34,'HASP':28,'WNK':31,'BUB':30,'ULK':31}
    glu={'AGC':49,'CAMK':47,'CK1':44,'CMGC':45,'NEK':48,'RGC':45,'STE':47,'TKL':46,'TYR':51,'HASP':52,'WNK':49,'BUB':39,'ULK':48}
    glu4={'AGC':53,'CAMK':51,'CK1':48,'CMGC':49,'NEK':52,'RGC':49,'STE':51,'TKL':50,'TYR':55,'HASP':56,'WNK':53,'BUB':43,'ULK':52}
    xdfg={'AGC':141,'CAMK':141,'CK1':140,'CMGC':149,'NEK':144,'RGC':139,'STE':140,'TKL':141,'TYR':145,'HASP':203,'WNK':148,'BUB':154,'ULK':142}
    asp={'AGC':142,'CAMK':142,'CK1':141,'CMGC':150,'NEK':145,'RGC':140,'STE':141,'TKL':142,'TYR':146,'HASP':204,'WNK':149,'BUB':155,'ULK':143}
    phe={'AGC':143,'CAMK':143,'CK1':142,'CMGC':151,'NEK':146,'RGC':141,'STE':142,'TKL':143,'TYR':147,'HASP':205,'WNK':150,'BUB':156,'ULK':144}

    alk_lys=rre_glu=rre4=hrd_his=x_dfg=dfg_asp=dfg_phe=99999
    alk_lys_res=rre_glu_res=rre4_res=hrd_his_res=x_dfg_res=dfg_asp_res=dfg_phe_res='XXX'

    hmm_result=SearchIO.read(f'{pwd}/server/{pdbfilename[0:-4]}{chain_id}_{group}.hmmer.txt',format='hmmer3-text')
    for hits in hmm_result:     #extract hit from alignment in HMM output file
        for hsps in hits:
            col_num=0;hmm_index=hsps.query_start;hit_index=hsps.hit_start+first_res-1
            for hmm_res in hsps.aln[0]:
                col_num=col_num+1
                if hmm_res!='.':
                    hmm_index=hmm_index+1
                if hsps.aln[1][col_num-1]!='-':
                    hit_index=hit_index+1
                if hmm_index==lys[group] and hmm_res!='.':
                    alk_lys=hit_index
                    alk_lys_res=list({hsps.aln[1][col_num-1]})[0]       #The residue is extracted using list otherwise it is printed in ''
                if hmm_index==glu[group] and hmm_res!='.':
                    rre_glu=hit_index
                    rre_glu_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==glu4[group] and hmm_res!='.':
                    rre4=hit_index
                    rre4_res=list({hsps.aln[1][col_num-1]})[0]
                #if hmm_index==his[group] and hmm_res!='.':
                #    hrd_his=hit_index
                #    hrd_his_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==xdfg[group] and hmm_res!='.':
                    x_dfg=hit_index
                    x_dfg_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==asp[group] and hmm_res!='.':
                    #fhandle_motifs.write(f'{hits.id} {col_num} {hmm_res} {hmm_index} {hit_index} {hsps.aln[1][col_num-1]}\n')
                    dfg_asp=hit_index
                    dfg_asp_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==phe[group] and hmm_res!='.':
                    dfg_phe=hit_index
                    dfg_phe_res=list({hsps.aln[1][col_num-1]})[0]
    #print(list(dfg_phe_res)[0],list(dfg_asp_res)[0])
    return (alk_lys,rre_glu,rre4,hrd_his,x_dfg,dfg_asp,dfg_phe,alk_lys_res,rre_glu_res,rre4_res,hrd_his_res,x_dfg_res,dfg_asp_res,dfg_phe_res)

def compute_distance(pwd,pdbfilename,chain_id,res1,res2,res1_type,res2_type,setting):
    if setting=='rre4-phe':
        phe_atom_type=identify_phe_atom(res2_type)
        dis_phe_rre4=distance_atoms(pwd,pdbfilename,chain_id,res1,res2,'CA',phe_atom_type)      #Change atom names for other residue types
        return dis_phe_rre4

    if setting=='lys-phe':
        phe_atom_type=identify_phe_atom(res2_type)
        dis_phe_lys=distance_atoms(pwd,pdbfilename,chain_id,res1,res2,'CA',phe_atom_type)      #Change atom names for other residue types
        return dis_phe_lys

    if setting=='saltbridge':
        dis_sb=distance_atoms(pwd,pdbfilename,chain_id,res1,res2,'CB','CB')      #Change atom names for other residue types
        return dis_sb

def identify_phe_atom(res2_type):
    if res2_type=='F':
        return 'CZ'
    if res2_type=='L' or res2_type=='P':
        return 'CG'
    if res2_type=='M':
        return 'CE'
    if res2_type=='S':
        return 'OG'
    if res2_type=='V':
        return 'CB'
    if res2_type=='W':
        return 'CZ3'
    if res2_type=='Y':
        return 'OH'
    if res2_type=='A':
        return 'CB'

def distance_atoms(pwd,pdbfilename,chain_id,res1,res2,atm1,atm2):
    if '.pdb' in pdbfilename.lower():
        parser=PDB.PDBParser()
    if '.cif' in pdbfilename.lower():
        parser=PDB.MMCIFParser()
    structure=parser.get_structure('PDB',(pwd+'/server/uploads/'+pdbfilename))
    atom_present=0; res1=int(res1); res2=int(res2)
    print(res1,res2,atm1,atm2)
    for model in structure:
        for chain in model:
            if chain.id==chain_id:
                for residue in chain:
                    if int(residue.id[1])==int(res1) and residue.get_id()[0]==' ':
                        if residue.has_id(atm1):
                            #print(res1,residue.id,chain[res1])
                            residue1=chain[res1]
                            atom_present=atom_present+1
                    if int(residue.id[1])==int(res2) and residue.get_id()[0]==' ':
                        if residue.has_id(atm2):
                            residue2=chain[res2]
                            atom_present=atom_present+1

    if atom_present==2:
        distance=np.round((residue1[atm1]-residue2[atm2]),2)
        #print(distance)
        #print(f'{distance:.2f}')
        return distance
    else:
        return 999

#if __name__ == '__main__':
#    pdbfilename=sys.argv[1]
#    (group,xdfg,dfg_asp,dfg_phe,xdfg_res,dfg_asp_res,dfg_phe_res,xdfg_phi,xdfg_psi,dfg_asp_phi,dfg_asp_psi,dfg_phe_phi,dfg_phe_psi,dfg_phe_chi1,chelix_conf,dfg_label,dfg_bkbone)=identify_state(pdbfilename)
#    print('Conf ',group,pdbfilename,chelix_conf,dfg_label,dfg_bkbone)
