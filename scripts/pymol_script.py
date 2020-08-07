#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:21:17 2020

@author: vivekmodi
"""
import subprocess,os

def pymol_genes(pwd,pymol_sessions,kinasechains_renumber_uniprot,pdb_domain_dict,domain_group_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Writing Pymol scripts for each domain...')
    domain_list=set(pdb_domain_dict.values())           #set returns unique elements from the list
    
    fhandle_pymol_bash=open((pymol_sessions+'/run_pymol.sh'),'w')
    for domains in domain_list:                 #Run a loop over each domain
        group=domain_group_dict[domains]
        
        if os.path.isfile(pymol_sessions+'/'+group+'_'+domains+'.pse.zip'):
            continue
        
        fhandle=open((pymol_sessions+'/'+group+'_'+domains+'.pml'),'w')
        fhandle.write("bg_color white\n")
        object_list=list()
        for pdbs in pdb_domain_dict:
            
            if pdb_domain_dict[pdbs]==domains:
                dfg_phe=int(domain_dfgnum_dict[domains]);dfg_asp=dfg_phe-1;xdfg=dfg_asp-1
                filename=(pdbs+'.cif.gz')
                obj_name=(pdb_spatial_dict[pdbs]+'-'+pdb_dihedral_dict[pdbs]+'-'+pdbs)
                object_list.append(obj_name)
                #obj_name=str(pdbs)
                fhandle.write("load "+kinasechains_renumber_uniprot+"/"+filename+"\n")
                fhandle.write("set_name "+pdbs+", "+obj_name+"\n")
                fhandle.write("hide lines, "+obj_name+"\n")
                fhandle.write("show cartoon, "+obj_name+"\n")
                #fhandle.write("select res "+str(xdfg)+" and "+obj_name+" and not name n+c+o"+"\n")
                #fhandle.write("show sticks, sele\n")
                fhandle.write("select res "+str(dfg_asp)+" and "+obj_name+" and not name n+c+o"+"\n")
                fhandle.write("show sticks, sele\n")
                fhandle.write("select res "+str(dfg_phe)+" and "+obj_name+" and not name n+c+o"+"\n")
                fhandle.write("show sticks, sele\n")
                fhandle.write("spectrum count, rainbow,"+obj_name+"\n")
                #last_obj=obj_name
            #if uniprot[pdbs]==values:
            #    if pdbs=='6C9DA' or pdbs=='6C9DB':  #Theseus does not run on these two structures of MARK1, maybe because they have long Cter extensions
            #        continue
            #    tpdb="t_"+str(pdbs)+".pdb"
        object_list.sort()
        fhandle.write("remove hydrogens\nremove solvent\n")
        fhandle.write("hide spheres\nhide dots\nselect HETATM\nhide (sele)\n")
        fhandle.write("alignto "+object_list[0]+"\ncenter\n")
        fhandle.write("order *,yes\n")
        fhandle.write("save "+pymol_sessions+'/'+group+'_'+domains+".pse\n")
        fhandle.close
        #print(domains)
        fhandle_pymol_bash.write(f'pymol -c {pymol_sessions}/{group}_{domains}.pml\n')
        #cmd=('pymol -c '+pymol_sessions+'/'+group+'_'+domains+'.pml')
        #subprocess.call(cmd,shell=True)
    fhandle_pymol_bash.close()
    print('Creating Pymol sessions for each domain...')
    cmd=('bash '+pymol_sessions+'/run_pymol.sh')
    process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    process.communicate()
    process.wait()
    
    
    print('Compressing Pymol sessions for each domain...')
    domain_list=set(pdb_domain_dict.values())           #set returns unique elements from the list
    os.chdir(f'{pymol_sessions}')
    for domains in domain_list:                 #Run a loop over each domain
        group=domain_group_dict[domains]
        
        if os.path.isfile(pymol_sessions+'/'+group+'_'+domains+'.pse.zip'):
            continue
              
        cmd=(f'zip -r {group}_{domains}.pse.zip {group}_{domains}.pse')
        process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        process.wait()
        cmd=(f'rm -irf {group}_{domains}.pse')
        subprocess.call(cmd,shell=True)
    os.chdir(f'{pwd}')
    
    
def pymol_genes_script(pwd,pymol_sessions,kinasechains_renumber_uniprot,pdb_domain_dict,domain_group_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Writing Pymol scripts to share...')
    domain_list=set(pdb_domain_dict.values())           #set returns unique elements from the list
    
    for domains in domain_list:                 #Run a loop over each domain
        group=domain_group_dict[domains]
        
        if os.path.isfile(pwd+'/static/downloads/pymol-uniprot-scripts/'+group+'_'+domains+'.zip'):
            continue
        if not os.path.isdir(pwd+'/static/downloads/pymol-uniprot-scripts/'+group+'_'+domains):
            cmd=(f"mkdir {pwd}/static/downloads/pymol-uniprot-scripts/{group}_{domains}")
            subprocess.call(cmd,shell=True)
        
        fhandle=open((pwd+'/static/downloads/pymol-uniprot-scripts/'+group+'_'+domains+'/'+group+'_'+domains+'.pml'),'w')
        fhandle.write("bg_color white\n")
        object_list=list()
        for pdbs in pdb_domain_dict:
            
            if pdb_domain_dict[pdbs]==domains:
                cmd=(f'cp {pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz {pwd}/static/downloads/pymol-uniprot-scripts/{group}_{domains}')
                subprocess.call(cmd,shell=True)
                cmd=(f'gunzip -f {pwd}/static/downloads/pymol-uniprot-scripts/{group}_{domains}/{pdbs}.cif.gz')
                subprocess.call(cmd,shell=True)
                dfg_phe=int(domain_dfgnum_dict[domains]);dfg_asp=dfg_phe-1;xdfg=dfg_asp-1
                filename=(pdbs+'.cif')
                obj_name=(pdb_spatial_dict[pdbs]+'-'+pdb_dihedral_dict[pdbs]+'-'+pdbs)
                object_list.append(obj_name)
                #obj_name=str(pdbs)
                fhandle.write("load "+filename+"\n")
                fhandle.write("set_name "+pdbs+", "+obj_name+"\n")
                fhandle.write("hide lines, "+obj_name+"\n")
                fhandle.write("show cartoon, "+obj_name+"\n")
                #fhandle.write("select res "+str(xdfg)+" and "+obj_name+" and not name n+c+o"+"\n")
                #fhandle.write("show sticks, sele\n")
                fhandle.write("select res "+str(dfg_asp)+" and "+obj_name+" and not name n+c+o"+"\n")
                fhandle.write("show sticks, sele\n")
                fhandle.write("select res "+str(dfg_phe)+" and "+obj_name+" and not name n+c+o"+"\n")
                fhandle.write("show sticks, sele\n")
                fhandle.write("spectrum count, rainbow,"+obj_name+"\n")
                #last_obj=obj_name
            #if uniprot[pdbs]==values:
            #    if pdbs=='6C9DA' or pdbs=='6C9DB':  #Theseus does not run on these two structures of MARK1, maybe because they have long Cter extensions
            #        continue
            #    tpdb="t_"+str(pdbs)+".pdb"
        object_list.sort()
        fhandle.write("remove hydrogens\nremove solvent\n")
        fhandle.write("hide spheres\nhide dots\nselect HETATM\nhide (sele)\n")
        fhandle.write("alignto "+object_list[0]+"\ncenter\n")
        fhandle.write("order *,yes\n")
        #fhandle.write("save "+pymol_sessions+'/'+group+'_'+domains+".pse\n")
        fhandle.close
        #print(domains)
        #fhandle_pymol_bash.write(f'pymol -c {pymol_sessions}/{group}_{domains}.pml\n')
        #cmd=('pymol -c '+pymol_sessions+'/'+group+'_'+domains+'.pml')
        #subprocess.call(cmd,shell=True)
    #fhandle_pymol_bash.close()
    #cmd=('bash '+pymol_sessions+'/run_pymol.sh')
    #process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    #process.wait()
    
    
    print('Compressing Pymol share directories for each domain...')
    domain_list=set(pdb_domain_dict.values())           #set returns unique elements from the list
    os.chdir(f'{pwd}/static/downloads/pymol-uniprot-scripts/')
    for domains in domain_list:                 #Run a loop over each domain
        group=domain_group_dict[domains]
        
        if os.path.isfile(f'{pwd}/static/downloads/pymol-uniprot-scripts/{group}_{domains}.zip'):
            continue
        
        
        cmd=(f'zip -r {group}_{domains}.zip {group}_{domains}')
        process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        process.wait()
        #cmd=(f'rm -irf {group}_{domains}')
        #subprocess.call(cmd,shell=True)
    os.chdir(f'{pwd}')