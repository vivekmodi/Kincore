#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:55:45 2020

@author: vivekmodi
"""
from Bio import SearchIO
#from collections import defaultdict

def read_psiblast(psiblast_result, excluded_output,df):
    print("Reading psiblast output...")
    fhandle_psiblast=SearchIO.read(psiblast_result,"blast-xml")
    #fhandle_human=open(human_output,"w")
    #fhandle_nonhuman=open(nonhuman_output,"w")
    fhandle_excluded=open(excluded_output,"w")
   # fhandle_noresidue=open(psiblast_error,"w")
    index=0  
    hit_accession_list=list()
    #pdb_dict=defaultdict(list)
    for hits in fhandle_psiblast:        #hits object does not have any evalue, only hsp have evalues
        #if index==500:
        #    break

        for hsp in hits.hsps:       
            if hsp.evalue<5.0 and hsp.aln_span>125:
                
                if hits.accession in hit_accession_list:    #to make sure that the sequences which are split into two hsps are not repeated in the output
                    continue
                hit_accession_list.append(hits.accession)
               
                #if '3DJ7A'  not in hsp.hit_id:
                #   continue
                #if 'HUMAN' not in hsp.hit_description:
                #    continue
                if 'HUMAN' in hsp.hit_description.upper() or 'MOUSE' in hsp.hit_description.upper() or 'RATTUS' in hsp.hit_description.upper() or \
                    'SCROFA' in hsp.hit_description.upper() or 'BOVIN' in hsp.hit_description.upper() or 'XENLA' in hsp.hit_description.upper() or\
                    'DROME' in hsp.hit_description.upper() or 'MACACA' in hsp.hit_description.upper() or 'SHEEP' in hsp.hit_description.upper() or\
                    'DANIO' in hsp.hit_description.upper():
                    
                    if hsp.hit_description.find('|')!=-1:       #Skip fusion proteins; returns position of substring, -1 means not found
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))     #modified to remove the random string in the beginning
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue
                    uniprot_name=str(hsp.hit_description.split("<")[1].split(">")[0].split("(")[0])
                    
                    if uniprot_name in ('NA','D3DSX2_HUMAN','RIOK1_HUMAN','RIOK2_HUMAN','RIOK3_HUMAN','PAN3_DROME','SG196_MOUSE','SG196_DANRE'):
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue                           #Remove the case which do not have uniprot assigned, D3DSX2 in not the correct uniprot but comes up in some PDB annotation
                        
         #           if uniprot_name.find('HUMAN')==-1:              #-1 means string not found
         #               description=(" ".join(hsp.hit_description.split(' ')[1:]))     #modified to remove the random string in the beginning
         #               fhandle_nonhuman.write(hsp.hit_id+","+description+","+uniprot_name+","+str(hsp.evalue)+","+str(hsp.bitscore)+"\n")
         #               continue
                    
                    if hsp.hit_id[0:4]=='6T28' or hsp.hit_id[0:4]=='6T29':      #Chain id in these structures is AAA
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue
                    
                    if hsp.hit_id=='3LZBE' or hsp.hit_id=='3LZBF' or hsp.hit_id=='3LZBG' or hsp.hit_id=='3LZBH' or hsp.hit_id=='5CNOX' or hsp.hit_id=='6PYHA' or hsp.hit_id=='6PYHD':      #Entire domain in not resolved but sequence is present in SEQRES
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue
                    if hsp.hit_id=='5KHUQ' or hsp.hit_id=='6PXVA' or hsp.hit_id=='6PXVC' or hsp.hit_id=='6PXWA' or hsp.hit_id=='6PXWB' or hsp.hit_id[0:4]=='6EQI':                     #CA chain only - 5KHUQ; 6PXV/W are cryo EM full length insulin recpetor
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue
                    if hsp.hit_id=='6TLJS':          #domain is not properly folded - BUB1B_HUMAN
                        description=(" ".join(hsp.hit_description.split(' ')[1:]))
                        fhandle_excluded.write(hsp.hit_id[0:5]+','+description+'\n')
                        continue
                    
                    
                    chain_length=hsp.hit_description.split()[1]
                    method=hsp.hit_description.split()[2]
                    reso=hsp.hit_description.split()[3]
                    rvalue=hsp.hit_description.split()[4]
                    free_rvalue=hsp.hit_description.split()[5]
                    protein=' '.join(hsp.hit_description.split('<')[0].split()[7:])
                    uniprot_name=str(hsp.hit_description.split('<')[1].split('>')[0].split('(')[0])
                    sequence=''.join(str(hsp.hit.seq).split('-'))       #this prints incomplete sequence, only corressponding to the alignment
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
                        res1=0;res2=0
                        continue      #skip the structures which do not have residue information in pdbaa
          #              fhandle_noresidue.write('No residue number '+hsp.hit_id[0:5]+hsp.hit_description+"\n")
                        
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
                    #df.at[index,'Construct']=sequence    #this prints incomplete sequence, only corressponding to the alignment;
                    index=index+1
                    #description=(" ".join(hsp.hit_description.split(' ')[1:]))     #modified to remove the random string in the beginning
                    #fhandle_human.write(hsp.hit_id+","+description+","+uniprot_name+","+str(hsp.evalue)+","+str(hsp.bitscore)+"\n")
                    #count=count+1
                   

    #print("Writing hits in psiblast-human.csv...")
    #print("Total number of human chains in psiblast_human: "+str(count))
    #fhandle_human.close()
    #fhandle_nonhuman.close()
    fhandle_excluded.close()
    #fhandle_noresidue.close()   
    #print(len(pdb_uniprot_dict))
    return df
    
# if __name__ == '__main__':
#     psiblast_result=sys.argv[1];human_output=sys.argv[2];nonhuman_output=sys.argv[3];fusion_output=sys.argv[4]
#     (pdb_uniprot_dict,seqres_start,seqres_end)=read_psiblast(psiblast_result, human_output, nonhuman_output, fusion_output)