#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 16:22:10 2020

@author: vivekmodi
"""

fhandle=open('Database-input-2020-01-30.tab','r')
pdblist=list()
for lines in fhandle:
    lines=lines.strip();lines=lines.split()
    if lines[3]=='INSR_HUMAN':
        pdblist.append(lines[0])
fhandle.close()

pymol=open('pymol_images_insr.pml','w')

for pdbs  in pdblist:
    filename=('/home/vivekmodi/Applications/Flask/Kinases/kinasechains_renumber_uniprot/'+pdbs+'.cif.gz')
    outputfile=('/home/vivekmodi/Applications/Flask/Kinases/images/'+pdbs+'.png')
    pymol.write(f'load {filename}\n')
    
    if pdbs!='1GAGA':
        pymol.write(f'load /home/vivekmodi/Applications/Flask/Kinases/kinasechains_renumber_uniprot/1GAGA.cif.gz\n')
        pymol.write('alignto 1GAGA\ndelete 1GAGA\n')
    pymol.write("hide lines\n")
    pymol.write("hide nonbonded\n")
    pymol.write("show cartoon\n")
                #fhandle.write("select res "+str(xdfg)+" and "+obj_name+" and not name n+c+o"+"\n")
                #fhandle.write("show sticks, sele\n")
    
    #pymol.write("select res 1074+1178+1156+1164 and name C+O+N\n")
    #pymol.write("show sticks, sele\n")
    #pymol.write("orient sele\n")
    #pymol.write("spectrum count, rainbow\n")
    pymol.write("set cartoon_transparency, 0.7\n")
    pymol.write("set cartoon_color, skyblue\n")
    pymol.write("select res 1176+1177+1178+1179\nshow sticks, sele\n")
    pymol.write("color skyblue, (name C*)\n")
    pymol.write("remove hydrogens\nremove solvent\n")
    pymol.write("hide spheres\nhide dots\nselect HETATM\nhide (sele)\n")
    pymol.write('set label_position,(1,1,1)\nlabel resi 1178 and name CB, "%s%s" %(resi, resn)\nset label_size, 15\n')
    pymol.write("set ray_shadows,off\n")
    pymol.write("bg_color white\n")
    pymol.write("set_view (\
    -0.948819816,   -0.312447548,    0.046000522,\
     0.082422569,   -0.385588646,   -0.918980539,\
     0.304871321,   -0.868156612,    0.391607076,\
     0.000000000,    0.000000000,  -69.517074585,\
   -30.599418640,   36.248752594,   12.066499710,\
   -77.908790588,  216.942947388,  -20.000000000 )\n")
    pymol.write("set ambient = 0.3\n")
    pymol.write("set ray_opaque_background, 1\n")
    pymol.write(f"png {outputfile},width=900, height=900, dpi=300, ray=1\n")
    pymol.write("reinitialize\n")
pymol.close()
