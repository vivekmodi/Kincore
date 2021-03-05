#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:55:45 2020

@author: vivekmodi
"""
from Bio import SearchIO
from datetime import datetime

def read_psiblast(pwd,df,psiblast_result, excluded_output):
    print("Reading psiblast output...")
    fhandle_psiblast=SearchIO.read(psiblast_result,"blast-xml")
    fhandle_excluded=open(excluded_output,"w")
    fhandle_newuniprots=open(f'{pwd}/New_uniprots.txt','w')
  
    index=0  
    hit_accession_list=list()
    specie_list=['HUMAN','MOUSE','RATTUS','SCROFA','BOVIN','XENLA','DROME','MACACA','SHEEP','DANIO','RABIT','CHICK']
   
    for hits in fhandle_psiblast:        #hits object does not have any evalue, only hsps have evalues
        #if index>10:
        #    continue
      
        for hsp in hits.hsps:       
            if hsp.evalue<5.0 and hsp.aln_span>125:
                
                if hits.accession in hit_accession_list:    #to make sure that the sequences which are split into two hsps are not repeated in the output
                    continue
                hit_accession_list.append(hits.accession)
               
                if any(specie in hsp.hit_description.upper() for specie in specie_list):     #any in an inbuilt function
                    
                    if hsp.hit_description.find('|')!=-1:       #Skip fusion proteins; returns position of substring, -1 means not found
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))     #modified to remove the random string in the beginning
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue
                    uniprot_name=str(hsp.hit_description.split("<")[1].split(">")[0].split("(")[0])
                    
                    if uniprot_name in ('NA','D3DSX2_HUMAN','RIOK1_HUMAN','RIOK2_HUMAN','RIOK3_HUMAN','PAN3_DROME','SG196_MOUSE','SG196_DANRE','E0W1I1_PEDHC','D3ZKP6_RAT'):
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue                           #Remove the case which do not have uniprot assigned, D3DSX2 is not the correct uniprot but comes up in some PDB annotation
                        
                          
                    if hsp.hit_id[0:4] in ('6T28','6T29','3LZBE','3LZBF','3LZBG','3LZBH','5CNOX','6PYHA','6PYHD','6TLJS','6Z1T','6Z1Q','6Z83','6Z84','6YUL','6YUM'):     # Description of these structures is in the kinasepml notes file
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue
                   
                    if hsp.hit_description.split()[6]=='yes':    #Skip CA-only chains
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue
                    if hsp.hit_description.split()[2]=='EM':     #Skip EM structures
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue
                    
                    chain_length=hsp.hit_description.split()[1]
                    method=hsp.hit_description.split()[2]
                    reso=hsp.hit_description.split()[3]
                    rvalue=hsp.hit_description.split()[4]
                    free_rvalue=hsp.hit_description.split()[5]
                    protein=' '.join(hsp.hit_description.split('<')[0].split()[7:])
                    uniprot_name=str(hsp.hit_description.split('<')[1].split('>')[0].split('(')[0])
                    #sequence=''.join(str(hsp.hit.seq).split('-'))       #this prints incomplete sequence, only corressponding to the alignment
                    specie=str(hsp.hit_description.split('[')[1].split(']')[0])
                    
                    if rvalue=='NA':
                        rvalue=999;
                    if free_rvalue=='NA':
                        free_rvalue=999;
                    if reso=='NA'    :
                        reso=999
                        
                    try:    
                        residue_range=(hsp.hit_description.split("<")[1].split(">")[0].split("(")[1]).split(")")[0]     #assign structure residue numbers to dictionary
                        res1=int(residue_range.split("-")[0])       #These residue numbers in pdbaa come from sifts database
                        res2=int(residue_range.split("-")[1])
                    except IndexError:
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:]+','+description+'\n')
                        continue      #skip the structures which do not have residue information in pdbaa, e.g. 6NSLA
        
                        
                    if uniprot_name=='PAK7_HUMAN':          #These conditions are created to replace old uniprot names in pdbaa
                        uniprot_name='PAK5_HUMAN'
                    if uniprot_name=='M3KL4_HUMAN':
                        uniprot_name='M3K21_HUMAN'
                    if uniprot_name=='MLTK_HUMAN':
                        uniprot_name='M3K20_HUMAN'
                    if uniprot_name=='J3KP20_HUMAN':
                        uniprot_name='NTRK1_HUMAN'
                    if uniprot_name=='B5BUJ6_HUMAN':
                        uniprot_name='MKNK1_HUMAN'
                    if uniprot_name=='Q506Q0_HUMAN':
                        uniprot_name='JAK2_HUMAN'
                    if uniprot_name=='RK_BOVIN':
                        uniprot_name='GRK1_BOVIN'
                    if uniprot_name=='Q6DRK7_DANRE':
                        uniprot_name='Q4KMI8_DANRE'
                    if uniprot_name=='Q7TMU1_MOUSE':
                        uniprot_name='BTK_MOUSE'
                    if uniprot_name=='D3ZMK9_RAT':
                        uniprot_name='PRAG1_RAT'
                    if uniprot_name=='O76755_DROME':
                        uniprot_name='A1Z6I7_DROME'                 #Sifts is using sequence from this uniprot
                        
                    
                       
                    df.at[index,'PDBid']=hsp.hit_id
                    df.at[index,'Protein']=protein
                    df.at[index,'UniprotID']=uniprot_name
                    df.at[index,'ChainLen']=int(chain_length)
                    df.at[index,'StrBegin']=int(res1)
                    df.at[index,'StrEnd']=int(res2)
                    df.at[index,'Specie']=specie
                    df.at[index,'Resolution']=round(float(reso),2)
                    df.at[index,'Method']=method
                    df.at[index,'Rvalue']=round(float(rvalue),2)
                    df.at[index,'FreeRvalue']=round(float(free_rvalue),2)
                   
                    index=index+1
                    
                else:
                    fhandle_newuniprots.write(f"{hsp.hit_id} {hsp.hit_description}\n")
                
                   
    fhandle_excluded.close()
    fhandle_newuniprots.close()
    
    today=str(datetime.now())[0:10].strip()
    df.to_excel(f'pdbaa_psiblast_dir/df_psiblast_{today}.xlsx',index=False)   #Write excel and csv
    df.to_csv(f'pdbaa_psiblast_dir/df_psiblast_{today}.csv',sep='\t',index=False)
    
    return df