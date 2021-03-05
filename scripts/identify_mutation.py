#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:15:12 2020

@author: vivekmodi
"""
import gzip

#aadict={'G':'GLY','A':'ALA','V':'VAL','I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',\
#        'W':'TRP','S':'SER','T':'THR','N':'ASN','Q':'GLN','R':'ARG','H':'HIS','K':'LYS',\
#        'D':'ASP','E':'GLU','C':'CYS','P':'PRO','U':'SEC'}

def identify_mutation(pwd,df):      
    aadict={'GLY':'G','ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y',\
        'TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','ARG':'R','HIS':'H','LYS':'K',\
        'ASP':'D','GLU':'E','CYS':'C','PRO':'P','SEC':'U','TPO':'T','CME':'C','CSS':'C',\
        'MSE':'M','OCY':'C','PTR':'Y','SEP':'S','CAF':'C','LGY':'K','CAS':'C','CSO':'C','CSX':'C',\
        'MK8':'E','NEP':'H','NMM':'R','CSD':'C','CYO':'Y','OCS':'C','OCY':'C','SCS':'C','ALY':'A',\
        'KCX':'K','MHO':'M','T8L':'T','CY0':'C','UNK':'X'}
    print('Identifying mutations in X-D-F residues...')
    kinasesifts=f'{pwd}/kinasesifts'
    log=open(f'{pwd}/kinasepml.log','a')
    
    fhandle_modified_aa=open(f'{pwd}/List_modified_aminoacid.txt','r')
        
    for i in df.index:
        pdbs=df.at[i,'PDBid']
        df.at[i,'Insertion']='None'
        df.at[i,'Linker']='None'
        df.at[i,'X_mut']='None'
        df.at[i,'Asp_mut']='None'
        df.at[i,'Phe_mut']='None'
        df.at[i,'Gly_mut']='None'
        df.at[i,'Chain_mut']='None'
        df.at[i,'Chain_phos']='None'
        df.at[i,'Modified_aa']='None'     
        df.at[i,'First_obs_res']=0
        df.at[i,'Last_obs_res']=0
        df.at[i,'XDFGresolved']='Yes'
        
        
        
        try:
            handle=gzip.open((kinasesifts+'/'+pdbs[0:4].lower()+'.csv.gz'),'rt')        
        except:
            log.write(f'identify_mutation: File not found {pdbs[0:4]}.csv.gz\n')
            continue
        
        chain=str(pdbs[4:]) 
        dfg_phe=int(df.at[i,'DFGnum'])
        dfg_asp=dfg_phe-1;xdfg=dfg_phe-2;dfg_gly=dfg_phe+1
        first_resi=0;mut_string_chain=list();phos_string_chain=list();modified_aa_chain=list()
        inserted_string=list();linker_string=list()
        
        for lines in handle:
            lines=lines.strip();lines=lines.split(',')
            
            
            if lines[0]=='PDBe_resnum':     #skip the first line
                continue
            if str(lines[4])==chain:
                
                if 'Insertion' in lines[8]:
                    inserted_resi_num=int(lines[5][0:-1])
                    inserted_resi_name=aadict[lines[3]]
                    inserted_string.append(inserted_resi_name+str(inserted_resi_num))
                    df.at[i,'Insertion']=','.join(inserted_string)
                    continue

                if 'Linker' in lines[8]:
                    linker_resi_num=int(lines[5][0:-1])
                    linker_resi_name=aadict[lines[3]]
                    linker_string.append(linker_resi_name+str(linker_resi_num))
                    df.at[i,'Linker']=','.join(linker_string)
                    continue

                if int(lines[5])==int(xdfg) and 'Engineered mutation' in lines[8]:
                    pdbResname=aadict[lines[3]]
                    mut_string=f'{lines[6]}{lines[5]}{pdbResname}'
                    df.at[i,'X_mut']=mut_string

                if int(lines[5])==int(dfg_asp) and 'Engineered mutation' in lines[8]:
                    pdbResname=aadict[lines[3]]
                    mut_string=f'{lines[6]}{lines[5]}{pdbResname}'
                    df.at[i,'Asp_mut']=mut_string

                if int(lines[5])==int(dfg_phe) and 'Engineered mutation' in lines[8]:
                    pdbResname=aadict[lines[3]]
                    mut_string=f'{lines[6]}{lines[5]}{pdbResname}'
                    df.at[i,'Phe_mut']=mut_string

                if int(lines[5])==int(dfg_gly) and 'Engineered mutation' in lines[8]:
                    pdbResname=aadict[lines[3]]
                    mut_string=f'{lines[6]}{lines[5]}{pdbResname}'
                    df.at[i,'Gly_mut']=mut_string

                if int(lines[5])!=-9999 and 'Engineered mutation' in lines[8]:
                    pdbResname=aadict[lines[3]]
                    mut_string_chain.append(f'{lines[6]}{lines[5]}{pdbResname}')
                    df.at[i,'Chain_mut']=','.join(mut_string_chain)

                if int(lines[5])!=-9999 and (lines[3]=='TPO' or lines[3]=='SEP' or lines[3]=='PTR'):
                    pdbResname=aadict[lines[3]]
                    phos_string_chain.append(f'{pdbResname}{lines[5]}')
                    df.at[i,'Chain_phos']=','.join(phos_string_chain)

                fhandle_modified_aa.seek(0)
                if int(lines[5])!=-9999 and (lines[3]+'\n') in fhandle_modified_aa.readlines():
                    pdbResname=aadict[lines[3]]
                    modified_aa_chain.append(f'{pdbResname}{lines[5]}')
                    df.at[i,'Modified_aa']=','.join(modified_aa_chain)

                if int(lines[5])!=-9999 and 'Not_Observed' not in str(lines[8]):
                    if first_resi==0:
                        df.at[i,'First_obs_res']=lines[5]
                        first_resi=1
                    df.at[i,'Last_obs_res']=lines[5]

                if (int(lines[5])==int(xdfg) or int(lines[5])==int(dfg_asp) or int(lines[5])==int(dfg_phe) or int(lines[5])==int(dfg_gly)) and 'Not_Observed' in lines[8]:
                    df.at[i,'XDFGresolved']='No'

                #if 'Engineered mutation' in lines[8]:
                #    pdbResname=aadict[lines[3]]
                #    mut_string_chain.append(f'{lines[6]}{lines[5]}{pdbResname}')
                #    df.at[i,'Chain_mut']=','.join(mut_string_chain)

        handle.close()
    fhandle_modified_aa.close()
    
    log.close()    
    return (df)