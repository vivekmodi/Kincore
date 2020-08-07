#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:37:23 2020

@author: vivekmodi
"""

import subprocess,os

def pymol_ligands(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Writing Pymol scripts and creating sessions for each ligand...')
    ligandList=list()
    for values in pdb_ligand_dict.values():
        values=values.split(',')
        for items in values:
            ligandList.append(items)
    #ligandList=set(pdb_ligand_dict.values())
    ligandList=set(ligandList)
    fhandle_pymol_bash=open((pwd+'/static/downloads/pymol-ligands'+'/run_pymol.sh'),'w')
        
    for ligands in ligandList:
        if ligands=="No_ligand":
            continue
        if os.path.isfile(pwd+'/static/downloads/pymol-ligands/'+ligands+'.pse.zip'):
            continue
        
        fhandle=open((pwd+'/static/downloads/pymol-ligands/'+ligands+'.pml'),'w')
        fhandle.write("bg_color white\n")
        object_list=list()
        for pdbs in pdb_ligand_dict:
            domainName=pdb_domain_dict[pdbs]
            
            if ligands in pdb_ligand_dict[pdbs]:
                dfg_phe=int(domain_dfgnum_dict[domainName]);dfg_asp=dfg_phe-1;xdfg=dfg_asp-1
                filename=(pdbs+'.cif.gz')
                obj_name=(pdb_spatial_dict[pdbs]+'-'+pdb_dihedral_dict[pdbs]+'-'+pdb_gene_dict[pdbs]+'-'+pdbs)
                object_list.append(obj_name)
                #obj_name=str(pdbs)
                fhandle.write("load "+pwd+"/kinasechains_renumber_uniprot/"+filename+"\n")          #Each conformation should be different color
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
        fhandle.write("hide spheres\nhide dots\n")
        fhandle.write("alignto "+object_list[0]+"\ncenter\n")
        fhandle.write("order *,yes\n")
        fhandle.write("save "+pwd+'/static/downloads/pymol-ligands/'+ligands+".pse\n")
        fhandle.close
        #print(domains)
        fhandle_pymol_bash.write(f'pymol -c {pwd}/static/downloads/pymol-ligands/{ligands}.pml\n')
        #cmd=('pymol -c '+pymol_sessions+'/'+group+'_'+domains+'.pml')
        #subprocess.call(cmd,shell=True)
    fhandle_pymol_bash.close()
    cmd=('bash '+pwd+'/static/downloads/pymol-ligands/run_pymol.sh')
    process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    process.communicate()
    process.wait()
    print(process.returncode)
    
def pymol_ligands_session_compress(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Compressing Pymol ligand sessions...')
    ligandList=list()
    for values in pdb_ligand_dict.values():
        values=values.split(',')
        for items in values:
            ligandList.append(items)
    #ligandList=set(pdb_ligand_dict.values())
    ligandList=set(ligandList)
    #fhandle_pymol_bash=open((pwd+'/static/downloads/pymol-ligands'+'/run_pymol.sh'),'w')
        
    for ligands in ligandList:
        if ligands=="No_ligand":
            continue
        if os.path.isfile(pwd+'/static/downloads/pymol-ligands/'+ligands+'.pse.zip'):
            continue
        if os.path.isfile(pwd+'/static/downloads/pymol-ligands/'+ligands+'.pse'):
            continue
        else:
            os.chdir(f'{pwd}/static/downloads/pymol-ligands')
            cmd=(f'zip -r {ligands}.pse.zip {ligands}.pse')
            process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            process.communicate()
            process.wait()
            cmd=(f'rm -irf {ligands}.pse')
            subprocess.call(cmd,shell=True)
            os.chdir(f'{pwd}')
        
    
def pymol_ligands_scripts(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Writing Pymol scripts for each ligand...')
    ligandList=list()
    for values in pdb_ligand_dict.values():
        values=values.split(',')
        for items in values:
            ligandList.append(items)
    #ligandList=set(pdb_ligand_dict.values())
    ligandList=set(ligandList)
    #fhandle_pymol_bash=open((pwd+'/static/downloads/pymol-ligands-scripts'+'/run_pymol.sh'),'w')
        
    for ligands in ligandList:
        if ligands=="No_ligand":
            continue
        ligandDir=(pwd+'/static/downloads/pymol-ligands-scripts/'+ligands)
        if not os.path.isdir(ligandDir):
            cmd=("mkdir "+ligandDir)
            subprocess.call(cmd,shell=True)
        
        if os.path.isfile(f'{ligandDir}/{ligands}.pml'):
            continue
        fhandle=open((ligandDir+'/'+ligands+'.pml'),'w')
        fhandle.write("bg_color white\n")
        object_list=list()
        for pdbs in pdb_ligand_dict:
            domainName=pdb_domain_dict[pdbs]
            
            if ligands in pdb_ligand_dict[pdbs]:
                cmd=(f'cp {pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz {ligandDir}')
                subprocess.call(cmd,shell=True)
                cmd=(f'gunzip -f {ligandDir}/{pdbs}.cif.gz')
                subprocess.call(cmd,shell=True)
                dfg_phe=int(domain_dfgnum_dict[domainName]);dfg_asp=dfg_phe-1;xdfg=dfg_asp-1
                filename=(pdbs+'.cif')
                obj_name=(pdb_spatial_dict[pdbs]+'-'+pdb_dihedral_dict[pdbs]+'-'+pdb_gene_dict[pdbs]+'-'+pdbs)
                object_list.append(obj_name)
                #obj_name=str(pdbs)
                fhandle.write("load "+filename+"\n")          #Each conformation should be different color
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
        fhandle.write("hide spheres\nhide dots\n")
        fhandle.write("alignto "+object_list[0]+"\ncenter\n")
        fhandle.write("order *,yes\n")
        #fhandle.write("save "+pwd+'/static/downloads/pymol-ligands/'+ligands+".pse\n")
        fhandle.close
        
def pymol_ligands_scripts_compress(pwd,pdb_ligand_dict,pdb_domain_dict,pdb_gene_dict,domain_dfgnum_dict,pdb_spatial_dict,pdb_dihedral_dict):
    print('Compressing Pymol ligand directories...')
    ligandList=list()
    for values in pdb_ligand_dict.values():
        values=values.split(',')
        for items in values:
            ligandList.append(items)
    
    ligandList=set(ligandList)
        
    for ligands in ligandList:
        if ligands=="No_ligand":
            continue
        
        os.chdir(f'{pwd}/static/downloads/pymol-ligands-scripts')
        if os.path.isdir(f'{ligands}.zip'):
            continue
        cmd=(f'zip -r {ligands}.zip {ligands}')
        process=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        process.communicate()
        process.wait()
        cmd=(f'rm -irf {ligands}')
        subprocess.call(cmd,shell=True)
        os.chdir(f'{pwd}')
    