#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:38:48 2020

@author: vivekmodi
"""
import subprocess, os
from Bio import SeqIO, SearchIO

def run_hmmsearch(pwd,uniname):
    #hmmfilename=os.path.join('/home/vivekmodi/Applications/Flask/Kinases/server/uploads/'+'Human-PK.hmm')
    cmd=(f'hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_AGC.hmmer.txt {pwd}/HMM_nonhuman/AGC.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_CAMK.hmmer.txt {pwd}/HMM_nonhuman/CAMK.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_CK1.hmmer.txt {pwd}/HMM_nonhuman/CK1.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_CMGC.hmmer.txt {pwd}/HMM_nonhuman/CMGC.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_NEK.hmmer.txt {pwd}/HMM_nonhuman/NEK.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_RGC.hmmer.txt {pwd}/HMM_nonhuman/RGC.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_STE.hmmer.txt {pwd}/HMM_nonhuman/STE.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_TKL.hmmer.txt {pwd}/HMM_nonhuman/TKL.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;\
         hmmsearch -o {pwd}/HMM_nonhuman/{uniname}_TYR.hmmer.txt {pwd}/HMM_nonhuman/TYR.hmm {pwd}/HMM_nonhuman/{uniname}.fasta;')
    subprocess.call(cmd,shell=True)
    return

def identify_group(pwd,uniname):
        hmm_result_AGC=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_AGC.hmmer.txt',format='hmmer3-text')
        hmm_result_CAMK=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_CAMK.hmmer.txt',format='hmmer3-text')
        hmm_result_CK1=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_CK1.hmmer.txt',format='hmmer3-text')
        hmm_result_CMGC=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_CMGC.hmmer.txt',format='hmmer3-text')
        hmm_result_NEK=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_NEK.hmmer.txt',format='hmmer3-text')
        hmm_result_RGC=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_RGC.hmmer.txt',format='hmmer3-text')
        hmm_result_STE=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_STE.hmmer.txt',format='hmmer3-text')
        hmm_result_TKL=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_TKL.hmmer.txt',format='hmmer3-text')
        hmm_result_TYR=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_TYR.hmmer.txt',format='hmmer3-text')

        eval=dict()
        eval['AGC']=100;eval['CAMK']=100;eval['CK1']=100;eval['CMGC']=100;eval['NEK']=100;eval['RGC']=100;eval['STE']=100;eval['TKL']=100;eval['TYR']=100;
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

        minEval=100;groupName='None'
        for groups in ('AGC','CAMK','CK1','CMGC','NEK','RGC','STE','TKL','TYR'):
            if eval[groups]<minEval:
                minEval=eval[groups]
                groupName=groups

        return groupName

def identify_residues(pwd,df,i):
    group=df.at[i,'Group']
    uniname=df.at[i,'UniprotID']
    gene=df.at[i,'Gene']
    
    lys=dict();glu=dict();glu4=dict();his_asp=dict();xdfg=dict();asp=dict();phe=dict();ape=dict();hinge1=dict();gtk=dict()
    lys['AGC']=30;lys['CAMK']=30;lys['CK1']=30;lys['CMGC']=30;lys['NEK']=30;lys['RGC']=29;lys['STE']=30;lys['TKL']=28;lys['TYR']=34;
    glu['AGC']=49;glu['CAMK']=47;glu['CK1']=44;glu['CMGC']=45;glu['NEK']=48;glu['RGC']=45;glu['STE']=47;glu['TKL']=46;glu['TYR']=51;
    glu4['AGC']=53;glu4['CAMK']=51;glu4['CK1']=48;glu4['CMGC']=49;glu4['NEK']=52;glu4['RGC']=49;glu4['STE']=51;glu4['TKL']=50;glu4['TYR']=55;
    his_asp['AGC']=124;his_asp['CAMK']=122;his_asp['CK1']=120;his_asp['CMGC']=130;his_asp['NEK']=127;his_asp['RGC']=122;his_asp['STE']=123;his_asp['TKL']=124;his_asp['TYR']=128;
    xdfg['AGC']=141;xdfg['CAMK']=141;xdfg['CK1']=140;xdfg['CMGC']=149;xdfg['NEK']=144;xdfg['RGC']=139;xdfg['STE']=140;xdfg['TKL']=141;xdfg['TYR']=145;
    asp['AGC']=142;asp['CAMK']=142;asp['CK1']=141;asp['CMGC']=150;asp['NEK']=145;asp['RGC']=140;asp['STE']=141;asp['TKL']=142;asp['TYR']=146;
    phe['AGC']=143;phe['CAMK']=143;phe['CK1']=142;phe['CMGC']=151;phe['NEK']=146;phe['RGC']=141;phe['STE']=142;phe['TKL']=143;phe['TYR']=147;
    ape['AGC']=170;ape['CAMK']=169;ape['CK1']=175;ape['CMGC']=175;ape['NEK']=172;ape['RGC']=168;ape['STE']=168;ape['TKL']=173;ape['TYR']=175;
    gtk['AGC']=78;gtk['CAMK']=76;gtk['CK1']=74;gtk['CMGC']=83;gtk['NEK']=77;gtk['RGC']=74;gtk['STE']=76;gtk['TKL']=75;gtk['TYR']=80;
    hinge1['AGC']=79;hinge1['CAMK']=77;hinge1['CK1']=75;hinge1['CMGC']=84;hinge1['NEK']=78;hinge1['RGC']=75;hinge1['STE']=77;hinge1['TKL']=76;hinge1['TYR']=81;
    
    alk_lys=rre_glu=rre4=hrd_asp=x_dfg=dfg_asp=dfg_phe=ape_glu=gtk_num=hinge1_num=99999
    alk_lys_res=rre_glu_res=rre4_res=hrd_asp_res=x_dfg_res=dfg_asp_res=dfg_phe_res=ape_res=gtk_res='XXX'

    hmm_result=SearchIO.read(f'{pwd}/HMM_nonhuman/{uniname}_{group}.hmmer.txt',format='hmmer3-text')
    for hits in hmm_result:     #extract hit from alignment in HMM output file
        first_in_align=0; DomainBegin=0; DomainEnd=0
        for domain_num,hsps in enumerate(hits):     #enumerate used if more than one domains are present in the sequence
            col_num=0;hmm_index=hsps.query_start;hit_index=hsps.hit_start
            #if first_in_align==0:                #the first residue of the alignment is assumed to be the beginning of domain
            DomainBegin=hsps.hit_start
            #    first_in_align=1
                         #the last residue of the alignment is assumed to be the end of domain
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
                if hmm_index==his_asp[group] and hmm_res!='.':
                    hrd_asp=hit_index
                    hrd_asp_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==xdfg[group] and hmm_res!='.':
                    x_dfg=hit_index
                    x_dfg_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==asp[group] and hmm_res!='.':
                    dfg_asp=hit_index
                    dfg_asp_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==phe[group] and hmm_res!='.':
                    dfg_phe=hit_index
                    dfg_phe_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==ape[group] and hmm_res!='.':
                    ape_glu=hit_index
                    ape_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==gtk[group] and hmm_res!='.':
                    gtk_num=hit_index
                    gtk_res=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==hinge1[group] and hmm_res!='.':
                    hinge1_num=hit_index
                DomainEnd=hit_index
            
            if hits.domain_obs_num==1:
                df.at[i,'ALKres']=alk_lys_res
                df.at[i,'ALKnum']=alk_lys
                df.at[i,'RREres']=rre_glu_res
                df.at[i,'RREnum']=rre_glu
                df.at[i,'HRDres']=hrd_asp_res
                df.at[i,'HRDnum']=hrd_asp
                df.at[i,'DFG_Aspres']=dfg_asp_res      #this column is not present for human structures
                df.at[i,'DFGres']=dfg_phe_res
                df.at[i,'DFGnum']=dfg_phe
                df.at[i,'APEres']=ape_res
                df.at[i,'APEnum']=ape_glu
                df.at[i,'GTKnum']=gtk_num
                df.at[i,'GTKres']=gtk_res
                df.at[i,'Hinge1']=hinge1_num                   
                df.at[i,'DomainBegin']=DomainBegin
                df.at[i,'DomainEnd']=DomainEnd
            
            if hits.domain_obs_num==2:     #if Uniprot sequence has two domains then check structure boundary condition
                res1=df.at[i,'StrBegin']
                res2=df.at[i,'StrEnd']
                if ((res1+res2)/2>=int(DomainBegin) and (res1+res2)/2<=int(DomainEnd)):
                    df.at[i,'Domain']=gene+str('_')+str(domain_num+1)
                    df.at[i,'ALKres']=alk_lys_res
                    df.at[i,'ALKnum']=alk_lys
                    df.at[i,'RREres']=rre_glu_res
                    df.at[i,'RREnum']=rre_glu
                    df.at[i,'HRDres']=hrd_asp_res
                    df.at[i,'HRDnum']=hrd_asp
                    df.at[i,'DFG_Aspres']=dfg_asp_res      #this column is not present for human structures
                    df.at[i,'DFGres']=dfg_phe_res
                    df.at[i,'DFGnum']=dfg_phe
                    df.at[i,'APEres']=ape_res
                    df.at[i,'APEnum']=ape_glu
                    df.at[i,'GTKnum']=gtk_num
                    df.at[i,'GTKres']=gtk_res
                    df.at[i,'Hinge1']=hinge1_num                   
                    df.at[i,'DomainBegin']=DomainBegin
                    df.at[i,'DomainEnd']=DomainEnd
                    
    #print(list(dfg_phe_res)[0],list(dfg_asp_res)[0])
    #return (alk_lys,rre_glu,rre4,hrd_asp,x_dfg,dfg_asp,dfg_phe,alk_lys_res,rre_glu_res,rre4_res,hrd_asp_res,x_dfg_res,dfg_asp_res,dfg_phe_res,ape_glu,ape_res,gtk_num,gtk_res,hinge1_num,DomainBegin,DomainEnd)
    return df

def identify_motifs_nonhuman(pwd,df):
    print('Identifying motifs for non-human proteins...')
    #fhandle_mapping=open(f'{pwd}/SwissProtIDGeneProteinMapping.csv','r')
    for i in df.index:
        uniname=df.at[i,'UniprotID']
        if 'HUMAN' not in uniname:
            #pdbfilename=df.at[i,'PDBid']
            #first_res=int(df.at[i,'First_obs_res'])
            if not os.path.isfile(f'{pwd}/HMM_nonhuman/{uniname}.fasta'):
                #fhandle_mapping.seek(0)
                #for lines in fhandle_mapping:
                #    lines=lines.strip();lines=lines.split('\t')
                uniseq=str(df.at[i,'UniSeq'])
                    #if lines[1]==uniname:
                fhandle_fasta=open(f'{pwd}/HMM_nonhuman/{uniname}.fasta','w')   #this should be uniprot sequence from SwissProt map which will also find domain boundary
                fhandle_fasta.write(f'>{uniname}\n')
                        #print_seq=df.loc[i,'StrSeq']
                fhandle_fasta.write(f"{uniseq}")
                fhandle_fasta.close()
                #break
            #df[df.PDBid==pdbfilename]['Sequence'].to_csv(f'{pwd}/HMM_nonhuman/{pdbfilename}.fasta',index=False,header='>')
            if not os.path.isfile(f'{pwd}/HMM_nonhuman/{uniname}_AGC.hmmer.txt'):
                run_hmmsearch(pwd, uniname)
            
            groupName=identify_group(pwd, uniname)
            df.at[i,'Group']=groupName
            
            if groupName!='NA':     #Sequence matches no group e.g. Q3ZC95_BOVIN
                df=identify_residues(pwd,df,i)
                # (alk_lys,rre_glu,rre4,hrd_asp,x_dfg,dfg_asp,dfg_phe,alk_lys_res,rre_glu_res,rre4_res,hrd_asp_res,x_dfg_res,dfg_asp_res,dfg_phe_res,ape_glu,ape_res,gtk_num,gtk_res,hinge1_num,DomainBegin,DomainEnd)=identify_residues(pwd, groupName, uniname)
                # df.at[i,'ALKres']=alk_lys_res
                # df.at[i,'ALKnum']=alk_lys
                # df.at[i,'RREres']=rre_glu_res
                # df.at[i,'RREnum']=rre_glu
                # df.at[i,'HRDres']=hrd_asp_res
                # df.at[i,'HRDnum']=hrd_asp
                # df.at[i,'DFG_Aspres']=dfg_asp_res      #this column is not present for human structures
                # df.at[i,'DFGres']=dfg_phe_res
                # df.at[i,'DFGnum']=dfg_phe
                # df.at[i,'APEres']=ape_res
                # df.at[i,'APEnum']=ape_glu
                # df.at[i,'GTKnum']=gtk_num
                # df.at[i,'GTKres']=gtk_res
                # df.at[i,'Hinge1']=hinge1_num
                # #df.at[i,'Hinge2']=hinge1_num+1
                # #df.at[i,'Hinge3']=hinge1_num+2
                
                # df.at[i,'Group']=groupName     #assign to closest group
                # df.at[i,'DomainBegin']=DomainBegin
                # df.at[i,'DomainEnd']=DomainEnd
    return df