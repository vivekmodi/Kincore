#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:52:35 2020

@author: vivekmodi
"""

import subprocess,os

def subListPymolSession(pwd,df):
    print('Creating Pymol sessions for Groups, Labels...')
    for organism in ('Human','All','Nonhuman'):
        for groups in ('AGC','CAMK','CK1','CMGC','NEK','STE','TKL','TYR','OTHER','.*'):
            group_out=groups
            for ligand_label in ('TypeI','TypeI½','TypeII','TypeIII','Allosteric','No_ligand','.*'):
                lig_out=ligand_label
                if ligand_label=='TypeI':
                    dontmatch='TypeI½'
                else:
                    dontmatch='X'

                for spatial in ('DFGin','DFGinter','DFGout','None','.*'):
                    if spatial=='DFGin':
                        for dihedral in ('BLAminus','BLAplus','ABAminus','BLBminus','BLBplus','BLBtrans','None','.*'):
                            if organism=='Human':
                                subListPymol=df[(df.Specie=='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Human_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName,subListPymol.count().PDBid)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='All':
                                subListPymol=df[(df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'All_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='Nonhuman':
                                subListPymol=df[(df.Specie!='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Nonhuman_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                    if spatial=='DFGinter':
                        for dihedral in ('BABtrans','None','.*'):
                            if organism=='Human':
                                subListPymol=df[(df.Specie=='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Human_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='All':
                                subListPymol=df[(df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'All_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='Nonhuman':
                                subListPymol=df[(df.Specie!='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Nonhuman_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                    if spatial=='DFGout':
                        for dihedral in ('BBAminus','None','.*'):
                            if organism=='Human':
                                subListPymol=df[(df.Specie=='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Human_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='All':
                                subListPymol=df[(df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'All_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='Nonhuman':
                                subListPymol=df[(df.Specie!='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))  & \
                                           (df.Spatial==spatial) & (df.Dihedral.str.match(dihedral))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if dihedral=='.*':
                                    dihedral='Alldihedrals'
                                outputName=f'Nonhuman_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                    if spatial=='.*':    #condition for Allspatials_Alldihedrals
                            if organism=='Human':
                                subListPymol=df[(df.Specie=='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if spatial=='.*':
                                    spatial='Allspatials'

                                outputName=f'Human_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='All':
                                subListPymol=df[(df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if spatial=='.*':
                                    spatial='Allspatials'

                                outputName=f'All_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)

                            if organism=='Nonhuman':
                                subListPymol=df[(df.Specie!='Homo sapiens') & (df.Group.str.match(groups)) & (df.Ligand_label.str.match(ligand_label) & \
                                            (df.Ligand_label.str.contains(dontmatch)==False))].sort_values('Resolution').sort_values('DomainBreak')
                                if group_out=='.*':   #changed to get the correct outputname
                                    group_out='Allgroups'
                                if lig_out=='.*':
                                    lig_out='All'
                                if spatial=='.*':
                                    spatial='Allspatials'

                                outputName=f'Nonhuman_{group_out}_{spatial}_{dihedral}_{lig_out}'
                                #print(outputName)

                                if subListPymol.count().PDBid==0:
                                    continue
                                if subListPymol.count().PDBid<=1000:        #too many structures in these cases
                                    createPymolSession(pwd,subListPymol,outputName)
                                pymolSessionScript(pwd,subListPymol,outputName)
                                createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
                                copyCoordinateFiles(pwd,subListPymol,outputName)




    #Pymol sessions for unique queries - PDB
    print('Creating Pymol sessions for PDBs...')
    pdbList=list()
    for i in df.index:
        pdb=df.at[i,'PDBid']
        pdbList.append(pdb[0:4])
    pdbList=set(pdbList)

    for pdbName in pdbList:
        for i in df.index:
            pdb=df.at[i,'PDBid']
            if pdbName in pdb:
                subListPymol=df[df.PDBid.str.match(pdbName)]
                group=subListPymol.at[i,'Group'];gene=subListPymol.at[i,'Gene']
                break

        outputName=f'{group}_{gene}_{pdbName}'
        #if not os.path.isfile(f'{pwd}/static/downloads/pymolSessions/{outputName}.pse.zip'):
        createPymolSession(pwd,subListPymol,outputName)
        pymolSessionScript(pwd,subListPymol,outputName)
        copyCoordinateFiles(pwd,subListPymol,outputName)

    #Pymol sessions for unique queries - Gene
    print('Creating Pymol sessions for genes...')
    for domain in sorted(set(df.Domain)):

        subListPymol=df[(df.Specie=='Homo sapiens') & (df.Domain==domain)]
        if len(subListPymol)>0:
            groupname=subListPymol.head(1).get('Group').values[0]
            outputName=f'Human_{groupname}_{domain}'
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            copyCoordinateFiles(pwd,subListPymol,outputName)


        subListPymol=df[(df.Domain==domain)]
        if len(subListPymol)>0:
            groupname=subListPymol.head(1).get('Group').values[0]
            outputName=f'All_{groupname}_{domain}'
            #print(outputName)
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            copyCoordinateFiles(pwd,subListPymol,outputName)

        subListPymol=df[(df.Specie!='Homo sapiens') & (df.Domain==domain)]
        if len(subListPymol)>0:
            groupname=subListPymol.head(1).get('Group').values[0]
            outputName=f'Nonhuman_{groupname}_{domain}'
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            createPymolSessionRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            pymolScriptRepresentative(pwd,subListPymol,f'Repr_{outputName}')
            copyCoordinateFiles(pwd,subListPymol,outputName)


    #Pymol sessions for unique queries - Ligand
    print('Creating Pymol sessions for ligands...')
    ligand_list=list()
    for i in df.index:
        if ',' in df.at[i,'Ligand']:
            lig_items=df.at[i,'Ligand'].split(',')
            for lig in lig_items:
                ligandname=lig.split(':')[0]
                if ligandname=='A' or ligandname=='B':   #1BKXA has a ligand named 'A'
                    continue
                ligand_list.append(ligandname)

        else:
            ligandname=df.at[i,'Ligand'].split(':')[0]
            if ligandname=='A' or ligandname=='B':
                    continue
            ligand_list.append(ligandname)

    ligand_list=set(ligand_list)

    for ligands in ligand_list:
        subListPymol=df[(df.Specie=='Homo sapiens') & (df.Ligand.str.match(ligands))]
        if len(subListPymol)>0:
            outputName=f'Human_{ligands}'
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            copyCoordinateFiles(pwd,subListPymol,outputName)

        subListPymol=df[(df.Ligand.str.match(ligands))]
        if len(subListPymol)>0:
            outputName=f'All_{ligands}'
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            copyCoordinateFiles(pwd,subListPymol,outputName)

        subListPymol=df[(df.Specie!='Homo sapiens') & (df.Ligand.str.match(ligands))]
        if len(subListPymol)>0:
            outputName=f'Nonhuman_{ligands}'
            createPymolSession(pwd,subListPymol,outputName)
            pymolSessionScript(pwd,subListPymol,outputName)
            copyCoordinateFiles(pwd,subListPymol,outputName)



def createPymolSession(pwd,subListPymol,outputName):
    #print('createPymolSession')
    kinasechainsRenumberUniprot=f'{pwd}/kinasechains_renumber_uniprot'
    #if os.path.isfile(f'{pwd}/static/downloads/pymolSessions/{outputName}.pse.zip'):
    #    return
    if len(subListPymol)==0:
        return

    fhandle=open(f'{pwd}/static/downloads/pymolSessions/{outputName}.pml','w')
    fhandle.write("bg_color black\n")
    fhandle.write("space rgb\n")

    objectList=list()

    for i in subListPymol.index:
        pdbs=subListPymol.at[i,'PDBid']
        group=subListPymol.at[i,'Group']
        domain=subListPymol.at[i,'Domain']
        spatial=subListPymol.at[i,'Spatial']
        dihedral=subListPymol.at[i,'Dihedral']
        ligand=subListPymol.at[i,'Ligand']
        dfgPhe=int(subListPymol.at[i,'DFGnum']);dfgAsp=dfgPhe-1;xdfg=dfgAsp-1
        ape=int(subListPymol.at[i,'APEnum'])
        hrd=int(subListPymol.at[i,'HRDnum'])
        color=subListPymol.at[i,'Color']
        filename=(pdbs+'.cif.gz')
        objName=f'{group}_{domain}_{spatial}_{dihedral}_{pdbs}'
        objectList.append(objName)
        fhandle.write(f"load {kinasechainsRenumberUniprot}/{filename}\n")
        fhandle.write("remove hydrogens\nremove solvent\n")
        fhandle.write("hide spheres\nhide dots\n")
        fhandle.write(f"set_name {pdbs}, {objName}\n")
        fhandle.write(f"hide lines, {objName}\n")
        fhandle.write(f"show cartoon, {objName}\n")
        #fhandle.write(f"set_color mycolor-{dihedral}, {color}\n")
        #fhandle.write(f"set_color light-grey, (230,230,230)\n")
        #fhandle.write(f"color mycolor-{dihedral}, {objName}\n")
        #fhandle.write(f"color light-grey, {objName}\n")
        fhandle.write(f"select res {xdfg}-{dfgPhe} and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select res {hrd} and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname PTR and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname SEP and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname TPO and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"spectrum count, rainbow,{objName}\n")
        #fhandle.write(f"color mycolor-{dihedral}, (name C*) and {objName}\n")
        #fhandle.write(f"color mycolor-{dihedral}, res {dfgAsp}-{ape} and {objName}\n")

        if ',' in ligand:
            for lig_items in ligand.split(','):
                ligname=lig_items.split(':')[0]
                fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                fhandle.write(f"color white, sele\n")
        else:
            ligname=ligand.split(':')[0]
            fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
            fhandle.write(f"color white, sele\n")


        fhandle.write(f"select res {xdfg}-{ape} and {objName}\n")
        fhandle.write(f"color magenta, sele\n")
        fhandle.write(f"color nitrogen, elem N\n")
        fhandle.write(f"color oxygen, elem O\n")


    #fhandle.write("hide cartoon\n")
    fhandle.write(f"alignto {objectList[0]}\ncenter\n")
    fhandle.write(f"set orthoscopic, on\n")
    fhandle.write("order *,yes\n")
    fhandle.write(f"save {pwd}/static/downloads/pymolSessions/{outputName}.pse\n")
    fhandle.close()
    cmd=(f'pymol -c {pwd}/static/downloads/pymolSessions/{outputName}.pml')
    process=subprocess.Popen(cmd,shell=True)
    process.communicate()
    process.wait()
    os.chdir(f'{pwd}/static/downloads/pymolSessions')
    cmd=(f'zip -r {outputName}.pse.zip {outputName}.pse')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()
    cmd=(f'rm -irf {outputName}.pse')
    subprocess.call(cmd,shell=True)
    os.chdir(f'{pwd}')



def pymolSessionScript(pwd,subListPymol,outputName):
    #print('pymolSessionScript')
    kinasechainsRenumberUniprot=f'{pwd}/kinasechains_renumber_uniprot'
    #if os.path.isfile(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}.zip'):
    #    return
    if len(subListPymol)==0:
        return

    if not os.path.isdir(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}'):
        os.mkdir(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}')


    for i in subListPymol.index:   #copy cif files

        pdbs=subListPymol.at[i,'PDBid']
        cmd=(f'cp {kinasechainsRenumberUniprot}/{pdbs}.cif.gz {pwd}/static/downloads/pymolSessionScripts/{outputName}')
        subprocess.call(cmd,shell=True)
        cmd=(f'gunzip -f {pwd}/static/downloads/pymolSessionScripts/{outputName}/{pdbs}.cif.gz')
        subprocess.call(cmd,shell=True)

    fhandle=open((f'{pwd}/static/downloads/pymolSessionScripts/{outputName}/{outputName}.pml'),'w')
    fhandle.write("bg_color black\n")
    fhandle.write("space rgb\n")

    objectList=list()

    for i in subListPymol.index:
        pdbs=subListPymol.at[i,'PDBid']
        group=subListPymol.at[i,'Group']
        domain=subListPymol.at[i,'Domain']
        spatial=subListPymol.at[i,'Spatial']
        dihedral=subListPymol.at[i,'Dihedral']
        ligand=subListPymol.at[i,'Ligand']
        dfgPhe=int(subListPymol.at[i,'DFGnum']);dfgAsp=dfgPhe-1;xdfg=dfgAsp-1
        ape=int(subListPymol.at[i,'APEnum'])
        hrd=int(subListPymol.at[i,'HRDnum'])
        color=subListPymol.at[i,'Color']
        filename=(pdbs+'.cif.gz')
        objName=f'{group}_{domain}_{spatial}_{dihedral}_{pdbs}'
        objectList.append(objName)
        fhandle.write(f"load {kinasechainsRenumberUniprot}/{filename}\n")
        fhandle.write("remove hydrogens\nremove solvent\n")
        fhandle.write("hide spheres\nhide dots\n")
        fhandle.write(f"set_name {pdbs}, {objName}\n")
        fhandle.write(f"hide lines, {objName}\n")
        fhandle.write(f"show cartoon, {objName}\n")
        #fhandle.write(f"set_color mycolor-{dihedral}, {color}\n")
        #fhandle.write(f"set_color light-grey, (230,230,230)\n")
        #fhandle.write(f"color mycolor-{dihedral}, {objName}\n")
        #fhandle.write(f"color light-grey, {objName}\n")
        fhandle.write(f"select res {xdfg}-{dfgPhe} and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select res {hrd} and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname PTR and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname SEP and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"select resname TPO and {objName}\n")
        fhandle.write("show sticks, sele\n")
        fhandle.write(f"spectrum count, rainbow,{objName}\n")
        #fhandle.write(f"color mycolor-{dihedral}, (name C*) and {objName}\n")
        #fhandle.write(f"color mycolor-{dihedral}, res {dfgAsp}-{ape} and {objName}\n")

        if ',' in ligand:
            for lig_items in ligand.split(','):
                ligname=lig_items.split(':')[0]
                fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                fhandle.write(f"color white, sele\n")
        else:
            ligname=ligand.split(':')[0]
            fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
            fhandle.write(f"color white, sele\n")

        #fhandle.write(f"color nitrogen, elem N\n")
        #fhandle.write(f"color oxygen, elem O\n")

        fhandle.write(f"select res {xdfg}-{ape} and {objName}\n")
        fhandle.write(f"color magenta, sele\n")
        fhandle.write(f"color nitrogen, elem N\n")
        fhandle.write(f"color oxygen, elem O\n")


    #fhandle.write("hide cartoon\n")
    fhandle.write(f"alignto {objectList[0]}\ncenter\n")
    fhandle.write(f"set orthoscopic, on\n")
    fhandle.write("order *,yes\n")
    fhandle.write(f"save {pwd}/static/downloads/pymolSessions/{outputName}/{outputName}.pse\n")
    fhandle.close()

    #Zip directories
    os.chdir(f'{pwd}/static/downloads/pymolSessionScripts/')
    cmd=(f'zip -r {outputName}.zip {outputName}')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    cmd=(f'rm -irf {outputName}')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{pwd}')

def createPymolSessionRepresentative(pwd,subListPymol,outputName):
    #print('createPymolSessionRepresentative')
    kinasechainsRenumberUniprot=f'{pwd}/kinasechains_renumber_uniprot'
    #if os.path.isfile(f'{pwd}/static/downloads/pymolSessions/{outputName}.pse.zip'):
    #    return
    if len(subListPymol)==0:
        return
    count=dict()
    for i in subListPymol.index:    #dictionary to count only 1 structure per domain
        uniprotid=subListPymol.at[i,'UniprotID']
        domain=subListPymol.at[i,'Domain']
        spatial=subListPymol.at[i,'Spatial']
        dihedral=subListPymol.at[i,'Dihedral']

        count[uniprotid,domain,spatial,dihedral]=0

    fhandle=open(f'{pwd}/static/downloads/pymolSessions/{outputName}.pml','w')
    fhandle.write("bg_color black\n")
    fhandle.write("space rgb\n")

    objectList=list()

    for i in subListPymol.index:
            pdbs=subListPymol.at[i,'PDBid']
            group=subListPymol.at[i,'Group']
            domain=subListPymol.at[i,'Domain']
            uniprotid=subListPymol.at[i,'UniprotID']
            spatial=subListPymol.at[i,'Spatial']
            dihedral=subListPymol.at[i,'Dihedral']
            ligand=subListPymol.at[i,'Ligand']
            dfgPhe=int(subListPymol.at[i,'DFGnum']);dfgAsp=dfgPhe-1;xdfg=dfgAsp-1
            ape=int(subListPymol.at[i,'APEnum'])
            hrd=int(subListPymol.at[i,'HRDnum'])
            color=subListPymol.at[i,'Color']

            if count[uniprotid,domain,spatial,dihedral]>=1:  #count only one structure per domain
                continue
            count[uniprotid,domain,spatial,dihedral]+=1
            filename=(pdbs+'.cif.gz')
            objName=f'{group}_{domain}_{spatial}_{dihedral}_{pdbs}'
            objectList.append(objName)
            fhandle.write(f"load {kinasechainsRenumberUniprot}/{filename}\n")
            fhandle.write("remove hydrogens\nremove solvent\n")
            fhandle.write("hide spheres\nhide dots\n")
            fhandle.write(f"set_name {pdbs}, {objName}\n")
            fhandle.write(f"hide lines, {objName}\n")
            fhandle.write(f"show cartoon, {objName}\n")
            #fhandle.write(f"set_color mycolor-{dihedral}, {color}\n")
            #fhandle.write(f"set_color light-grey, (230,230,230)\n")
            #fhandle.write(f"color mycolor-{dihedral}, {objName}\n")
            #fhandle.write(f"color light-grey, {objName}\n")
            fhandle.write(f"select res {xdfg}-{dfgPhe} and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select res {hrd} and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname PTR and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname SEP and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname TPO and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"spectrum count, rainbow,{objName}\n")
            #fhandle.write(f"color mycolor-{dihedral}, (name C*) and {objName}\n")
            #fhandle.write(f"color mycolor-{dihedral}, res {dfgAsp}-{ape} and {objName}\n")

            if ',' in ligand:
                for lig_items in ligand.split(','):
                    ligname=lig_items.split(':')[0]
                    fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                    fhandle.write(f"color white, sele\n")
            else:
                ligname=ligand.split(':')[0]
                fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                fhandle.write(f"color white, sele\n")

            fhandle.write(f"select res {xdfg}-{ape} and {objName}\n")
            fhandle.write(f"color magenta, sele\n")
            fhandle.write(f"color nitrogen, elem N\n")
            fhandle.write(f"color oxygen, elem O\n")


    #fhandle.write("hide cartoon\n")
    fhandle.write(f"alignto {objectList[0]}\ncenter\n")
    fhandle.write(f"set orthoscopic, on\n")
    fhandle.write("order *,yes\n")
    fhandle.write(f"save {pwd}/static/downloads/pymolSessions/{outputName}.pse\n")
    fhandle.close()
    cmd=(f'pymol -c {pwd}/static/downloads/pymolSessions/{outputName}.pml')
    process=subprocess.Popen(cmd,shell=True)
    process.communicate()
    process.wait()
    os.chdir(f'{pwd}/static/downloads/pymolSessions')
    cmd=(f'zip -r {outputName}.pse.zip {outputName}.pse')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()
    cmd=(f'rm -irf {outputName}.pse')
    subprocess.call(cmd,shell=True)
    os.chdir(f'{pwd}')

def pymolScriptRepresentative(pwd,subListPymol,outputName):
    #print('pymolScriptRepresentative')
    kinasechainsRenumberUniprot=f'{pwd}/kinasechains_renumber_uniprot'
    #if os.path.isfile(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}.zip'):
    #    return
    if len(subListPymol)==0:
        return

    if not os.path.isdir(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}'):
        os.mkdir(f'{pwd}/static/downloads/pymolSessionScripts/{outputName}')

    count=dict()
    for i in subListPymol.index:    #dictionary to count only 1 structure per domain
        uniprotid=subListPymol.at[i,'UniprotID']
        domain=subListPymol.at[i,'Domain']
        spatial=subListPymol.at[i,'Spatial']
        dihedral=subListPymol.at[i,'Dihedral']
        count[uniprotid,domain,spatial,dihedral]=0

    fhandle=open((f'{pwd}/static/downloads/pymolSessionScripts/{outputName}/{outputName}.pml'),'w')
    fhandle.write("bg_color black\n")
    fhandle.write("space rgb\n")

    objectList=list()


    for i in subListPymol.index:
            pdbs=subListPymol.at[i,'PDBid']
            group=subListPymol.at[i,'Group']
            domain=subListPymol.at[i,'Domain']
            uniprotid=subListPymol.at[i,'UniprotID']
            spatial=subListPymol.at[i,'Spatial']
            dihedral=subListPymol.at[i,'Dihedral']
            ligand=subListPymol.at[i,'Ligand']
            dfgPhe=int(subListPymol.at[i,'DFGnum']);dfgAsp=dfgPhe-1;xdfg=dfgAsp-1
            ape=int(subListPymol.at[i,'APEnum'])
            hrd=int(subListPymol.at[i,'HRDnum'])
            color=subListPymol.at[i,'Color']
            if count[uniprotid,domain,spatial,dihedral]>=1:  #count only one structure per domain
                continue
            count[uniprotid,domain,spatial,dihedral]+=1

            cmd=(f'cp {kinasechainsRenumberUniprot}/{pdbs}.cif.gz {pwd}/static/downloads/pymolSessionScripts/{outputName}')  #copy cif files
            process=subprocess.Popen(cmd,shell=True)
            process.wait()
            cmd=(f'gunzip -f {pwd}/static/downloads/pymolSessionScripts/{outputName}/{pdbs}.cif.gz')
            process=subprocess.Popen(cmd,shell=True)
            process.wait()

            filename=(pdbs+'.cif.gz')
            objName=f'{group}_{domain}_{spatial}_{dihedral}_{pdbs}'
            objectList.append(objName)
            fhandle.write(f"load {kinasechainsRenumberUniprot}/{filename}\n")
            fhandle.write("remove hydrogens\nremove solvent\n")
            fhandle.write("hide spheres\nhide dots\n")
            fhandle.write(f"set_name {pdbs}, {objName}\n")
            fhandle.write(f"hide lines, {objName}\n")
            fhandle.write(f"show cartoon, {objName}\n")
            #fhandle.write(f"set_color mycolor-{dihedral}, {color}\n")
            #fhandle.write(f"set_color light-grey, (230,230,230)\n")
            #fhandle.write(f"color mycolor-{dihedral}, {objName}\n")
            #fhandle.write(f"color light-grey, {objName}\n")
            fhandle.write(f"select res {xdfg}-{dfgPhe} and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select res {hrd} and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname PTR and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname SEP and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"select resname TPO and {objName}\n")
            fhandle.write("show sticks, sele\n")
            fhandle.write(f"spectrum count, rainbow,{objName}\n")
            #fhandle.write(f"color mycolor-{dihedral}, (name C*) and {objName}\n")
            #fhandle.write(f"color mycolor-{dihedral}, res {dfgAsp}-{ape} and {objName}\n")

            if ',' in ligand:
                for lig_items in ligand.split(','):
                    ligname=lig_items.split(':')[0]
                    fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                    fhandle.write(f"color white, sele\n")
            else:
                ligname=ligand.split(':')[0]
                fhandle.write(f"select resname {ligname}\nshow sticks, sele\n")
                fhandle.write(f"color white, sele\n")


            fhandle.write(f"select res {xdfg}-{ape} and {objName}\n")
            fhandle.write(f"color magenta, sele\n")
            fhandle.write(f"color nitrogen, elem N\n")
            fhandle.write(f"color oxygen, elem O\n")


    #fhandle.write("hide cartoon\n")
    fhandle.write(f"alignto {objectList[0]}\ncenter\n")
    fhandle.write(f"set orthoscopic, on\n")
    fhandle.write("order *,yes\n")
    fhandle.write(f"save {pwd}/static/downloads/pymolSessions/{outputName}/{outputName}.pse\n")
    fhandle.close()

    #Zip directories
    os.chdir(f'{pwd}/static/downloads/pymolSessionScripts/')
    cmd=(f'zip -r {outputName}.zip {outputName}')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    cmd=(f'rm -irf {outputName}')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{pwd}')

def copyCoordinateFiles(pwd,subList,outputName):
    #print('copyCoordinateFiles')
    coordDir=f'{pwd}/static/downloads/coordinateFiles'

    if len(subList)==0:
        return
    #if os.path.isfile(f'{coordDir}/{outputName}_uniNum.zip'):
    #    return
    #if os.path.isfile(f'{coordDir}/{outputName}_pdbNum.zip'):
    #    return
    #if os.path.isfile(f'{coordDir}/{outputName}_alignNum.zip'):
    #    return
    if not os.path.isdir(f'{coordDir}/{outputName}_uniNum'):
        os.mkdir(f'{coordDir}/{outputName}_uniNum')
    if not os.path.isdir(f'{coordDir}/{outputName}_pdbNum'):
        os.mkdir(f'{coordDir}/{outputName}_pdbNum')
    if not os.path.isdir(f'{coordDir}/{outputName}_alignNum'):
        os.mkdir(f'{coordDir}/{outputName}_alignNum')

    for i in subList.index:
        pdbs=subList.at[i,'PDBid'];group=subList.at[i,'Group']
        domain=subList.at[i,'Domain'];spatial=subList.at[i,'Spatial'];dihedral=subList.at[i,'Dihedral']

        filename=f'{group}_{domain}_{spatial}_{dihedral}_{pdbs}'
        cmd=(f'cp {pwd}/kinasechains_renumber_uniprot/{pdbs}.cif.gz {coordDir}/{outputName}_uniNum/{filename}.cif.gz')
        subprocess.call(cmd,shell=True)
        cmd=(f'cp {pwd}/kinasechains/{pdbs}.cif.gz {coordDir}/{outputName}_pdbNum/{filename}.cif.gz')
        subprocess.call(cmd,shell=True)
        cmd=(f'cp {pwd}/kinasechains_renumber_alignment/{pdbs}.cif.gz {coordDir}/{outputName}_alignNum/{filename}.cif.gz')
        subprocess.call(cmd,shell=True)


    os.chdir(f'{coordDir}/{outputName}_uniNum')
    cmd=(f'gunzip -f *')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{coordDir}/{outputName}_pdbNum')
    cmd=(f'gunzip -f *')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{coordDir}/{outputName}_alignNum')
    cmd=(f'gunzip -f *')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{coordDir}')
    cmd=(f'zip -r {outputName}_uniNum.zip {outputName}_uniNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    cmd=(f'rm -irf {outputName}_uniNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    #os.chdir(f'{pdb}')
    cmd=(f'zip -r {outputName}_pdbNum.zip {outputName}_pdbNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    cmd=(f'rm -irf {outputName}_pdbNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    #os.chdir(f'{alignment}')
    cmd=(f'zip -r {outputName}_alignNum.zip {outputName}_alignNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    cmd=(f'rm -irf {outputName}_alignNum')
    process=subprocess.Popen(cmd,shell=True)
    process.wait()

    os.chdir(f'{pwd}')
