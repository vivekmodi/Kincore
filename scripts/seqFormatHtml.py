#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:45:06 2020

@author: vivekmodi
"""
import gzip,sys,os,re
import pandas as pd


def format_seq_text(pwd,df):
    aadict={'GLY':'G','ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y',\
        'TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','ARG':'R','HIS':'H','LYS':'K',\
        'ASP':'D','GLU':'E','CYS':'C','PRO':'P','SEC':'U','TPO':'T','CME':'C','CSS':'C',\
        'MSE':'M','OCY':'C','PTR':'Y','SEP':'S','CAF':'C','LGY':'K','CAS':'C','CSO':'C','CSX':'C',\
        'MK8':'E','NEP':'H','NMM':'R','CSD':'C','CYO':'Y','OCS':'C','OCY':'C','SCS':'C','ALY':'A',\
        'KCX':'K','MHO':'M','T8L':'T','CY0':'C','UNK':'X'}

    

    print('Fomatting sequences for html...')
    for i in df.index:
        uni_pdb_dict=dict()
        pdbs=df.at[i,'PDBid']
        fhandle_sifts=(pwd+'/kinasesifts/'+pdbs[0:4].lower()+'.csv.gz')
        
        
     
        outputfile=(pwd+'/formattedSeq/'+pdbs[0:4]+'.seq')
        if os.path.isfile(outputfile):
            continue
        fhandle_output=open(outputfile,'w')
        fhandle_sifts.readline()
        firstResi='True';
        for lines in fhandle_sifts:
            lines=lines.strip();lines=lines.split(',');
            pdbResname=lines[3];pdbChain=lines[4];uniResnum=int(lines[5]);uniResname=lines[6];uniAcc=lines[7];anno=lines[8]    #variable from sifts csv file
            
            if uniResnum==-9999:
                continue
            
            if pdbs[4]!=pdbChain:
                continue
            
            if firstResi=='True':
                lastResidue=uniResnum
                firstResi='False'
                
            if uniResnum-lastResidue>1:            #if there is a gap in sequence numbering then assign '-'
                diff=uniResnum-lastResidue
                for i in range(diff):
                    lastResidue=lastResidue+1
                    uni_pdb_dict[lastResidue]='-'
                
            if uniResnum!=-9999:
                if 'Not_Observed' in anno:
                    uni_pdb_dict[uniResnum]=aadict[pdbResname].lower()              #assign lower case letters
                    lastResidue=uniResnum
                    
                else:
                    uni_pdb_dict[uniResnum]=aadict[pdbResname].upper()              #assign upper case letters
                    lastResidue=uniResnum
                    
            
        
        printNum=0;count=0
        for uniResnum in uni_pdb_dict:
            if uniResnum>printNum and lastResidue>=printNum+10:
                for i in range(9,51,10):
                    if uniResnum+i<=lastResidue:
                        printNum=uniResnum+i
                        fhandle_output.write(f'{printNum:>10}  ')
                fhandle_output.write('\n')
                        
            count=count+1
            if count%50==1 and (lastResidue-uniResnum)<10:
                fhandle_output.write(f'\n')
            if count%10==0 and count%50!=0:
                fhandle_output.write(f'{uni_pdb_dict[uniResnum]}')
                fhandle_output.write(f'  ')
            elif count%50==0:
                fhandle_output.write(f'{uni_pdb_dict[uniResnum]}')
                fhandle_output.write(f'\n')
            #elif count%50==1 and (lastResidue-printNum)<10:
            #    fhandle_output.write(f'\n')
            else:
                fhandle_output.write(f'{uni_pdb_dict[uniResnum]}')
            
        fhandle_output.close()
        fhandle_sifts.close()
  

def format_seq_html(pwd,df):
    aadict={'GLY':'G','ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y',\
        'TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','ARG':'R','HIS':'H','LYS':'K',\
        'ASP':'D','GLU':'E','CYS':'C','PRO':'P','SEC':'U','TPO':'T','CME':'C','CSS':'C',\
        'MSE':'M','OCY':'C','PTR':'Y','SEP':'S','CAF':'C','LGY':'K','CAS':'C','CSO':'C','CSX':'C',\
        'MK8':'E','NEP':'H','NMM':'R','CSD':'C','CYO':'Y','OCS':'C','OCY':'C','SCS':'C','ALY':'A',\
        'KCX':'K','MHO':'M','T8L':'T','CY0':'C','UNK':'X'}
    
    log=open(f'{pwd}/kinasepml.log','a')
    print('Fomatting sequences for html...')
    for i in df.index:
        pdbs=df.at[i,'PDBid']
       
        uni_pdb_dict=dict();
        siftsFile=(pwd+'/kinasesifts/'+pdbs[0:4].lower()+'.csv.gz')
        try:
            fhandle_sifts=gzip.open(siftsFile,'rt')
        except:
            log.write(f'seqFormatHtml: File not found {siftsFile}\n')
            continue
        
      
        outputfile=(pwd+'/formattedSeq/'+pdbs[0:4]+'.html')
        #if os.path.isfile(outputfile):
        #    continue
        fhandle_output=open(outputfile,'w')
        fhandle_output.write(f'<pre class="very-dark-text">Representative sequence - Chain {pdbs[4:]}<br>')
        fhandle_sifts.readline()
        firstResi='True'
        insertion_location=0;insertion_size=0
        for lines in fhandle_sifts:
            lines=lines.strip();lines=lines.split(',');
            if 'Insertion' in lines[8] or 'Linker' in lines[8]:
                pdbResname=lines[3];pdbChain=lines[4];uniResnum=int(lines[5][0:-1]);insertion_code=lines[5][-1];uniResname=lines[6];uniAcc=lines[7];anno=lines[8]
            else:
                pdbResname=lines[3];pdbChain=lines[4];uniResnum=int(lines[5]);uniResname=lines[6];uniAcc=lines[7];anno=lines[8];insertion_code=''
            
            if pdbs[4]!=pdbChain:
                continue
            if uniResnum==-9999:
                continue
            if firstResi=='True':
                lastResidue=uniResnum              #last residue is initialized here
                firstResi='False'
                
            if uniResnum-lastResidue>1:            #if there is a gap in sequence numbering then assign '-'
                diff=uniResnum-lastResidue
                #print(f'diff is {diff}')
                for i in range(diff):
                    lastResidue=lastResidue+1
                    uni_pdb_dict[str(lastResidue)]=str('-')
                    #lastResidue=lastResidue+1
                    #print(f'lastResidue is {i},{lastResidue}')
                    
                    #print(lastResidue,uni_pdb_dict[str(lastResidue)])
                #lastResidue=lastResidue+1
                
                    
            if 'Not_Observed' in anno:
                #uni_pdb_dict[uniResnum]=aadict[pdbResname].lower()
                #lastResidue=uniResnum
                if 'mutation' in anno:
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].lower()+'</span>'
                elif 'Insertion' in anno or 'Linker' in anno:
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].lower()+'</span>'
                    insertion_size=insertion_size+1
                    if insertion_location==0:
                        insertion_location=uniResnum
        #        elif 'Linker' in anno:
        #            uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].lower()+'</span>'
                elif pdbResname=='TPO' or pdbResname=='SEP' or pdbResname=='PTR':              #'Phosphorylation' in anno:
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#5cb85c	;">'+aadict[pdbResname].lower()+'</span>'
                else:
                    uni_pdb_dict[str(uniResnum)+insertion_code]=aadict[pdbResname].lower()
                lastResidue=uniResnum
                
            else:
                #uni_pdb_dict[uniResnum]=aadict[pdbResname].upper()
                #lastResidue=uniResnum
                if 'mutation' in anno:
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].upper()+'</span>'
                elif 'Insertion' in anno or 'Linker' in anno:
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].upper()+'</span>'
                    insertion_size=insertion_size+1
                    if insertion_location==0:
                        insertion_location=uniResnum
                #elif 'Linker' in anno:
                #    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#d9534f;">'+aadict[pdbResname].upper()+'</span>'
                elif pdbResname=='TPO' or pdbResname=='SEP' or pdbResname=='PTR':
                    uni_pdb_dict[str(uniResnum)+insertion_code]='<span style="color: white;background:#5cb85c	;">'+aadict[pdbResname].upper()+'</span>'
                else:
                    uni_pdb_dict[str(uniResnum)+insertion_code]=aadict[pdbResname].upper()
                lastResidue=uniResnum
        
        printNum=0;count=0;count_insert=0
        for uniResnum in uni_pdb_dict.keys():
            #print(uniResnum,uni_pdb_dict[str(uniResnum)])
            
            if re.search('[A-Za-z]',str(uniResnum)):           #check is insertion code is present
                insertion_code=str(uniResnum[-1])
                uniResnum=int(uniResnum[0:-1])
                count_insert=count_insert+1
                #print(uniResnum,count_insert,printNum,lastResidue)
                if int(uniResnum+count_insert)>printNum and lastResidue>=printNum+10:
                    for i in range(9,51,10):
                        if uniResnum+i+count_insert<=lastResidue:
                            printNum=uniResnum+i+count_insert
                            if uniResnum+i>=insertion_location:
                                fhandle_output.write(f'{printNum-insertion_size:>10}  ')
                            else:
                                fhandle_output.write(f'{printNum-count_insert:>10}  ')
                            #print(count_insert,printNum)
                    fhandle_output.write('<br>')
                            
                count=count+1
                if count%50==1 and (lastResidue-uniResnum+count_insert)<10:
                    fhandle_output.write(f'<br>')
                if count%10==0 and count%50!=0:
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
                    fhandle_output.write(f'  ')
                elif count%50==0:
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
                    fhandle_output.write(f'<br>')
                else:
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
                    
                
            else:
                insertion_code=''
                uniResnum=int(uniResnum)
                
                if int(uniResnum+count_insert)>printNum and lastResidue>=printNum+10:
                    for i in range(9,51,10):
                        if uniResnum+i<=lastResidue:
                            printNum=uniResnum+i+count_insert
                            if uniResnum+i>=insertion_location:
                                fhandle_output.write(f'{printNum-insertion_size:>10}  ')
                            else:
                                fhandle_output.write(f'{printNum-count_insert:>10}  ')
                    fhandle_output.write('<br>')
                            
                count=count+1
                if count%50==1 and (lastResidue-uniResnum+count_insert)<10:
                    fhandle_output.write(f'<br>')
                if count%10==0 and count%50!=0:
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
                    fhandle_output.write(f'  ')
                elif count%50==0:
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
                    fhandle_output.write(f'<br>')
                else:
                    #print(uniResnum,uni_pdb_dict[str(uniResnum)])
                    fhandle_output.write(f'{uni_pdb_dict[str(uniResnum)+insertion_code]}')
            
        fhandle_output.write(f'</pre>')
        fhandle_output.close()
        fhandle_sifts.close()
    log.close()

if __name__=='__main__':
    pwd='/home/vivekmodi/Applications/Flask/Kinases'
    df=pd.read_csv(sys.argv[1],sep='\t',header='infer')
    format_seq_html(pwd,df)