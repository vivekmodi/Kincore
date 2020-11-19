#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:40:32 2020

@author: vivekmodi
"""
import subprocess

def transfer_to_dunbrack3(pwd):
    cmd=(f'rsync -avr --progress {pwd}/formattedSeq/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/formattedSeq/;\
        rsync -avr --progress {pwd}/kinasechains/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/kinasechains/;\
        rsync -avr --progress {pwd}/kinasechains_dihedrals/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_dihedrals/;\
        rsync -avr --progress {pwd}/kinasechains_renumber_alignment/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_renumber_alignment/;\
        rsync -avr --progress {pwd}/kinasechains_renumber_uniprot/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_renumber_uniprot/;\
        rsync -avr --progress {pwd}/kinasecifs/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/kinasecifs/;\
        rsync -avr --progress {pwd}/static/kinasechainsNGL/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/kinasechainsNGL;\
        scp -r {pwd}/static/downloads/coordinateFiles/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/coordinateFiles/;\
        scp -r {pwd}/static/downloads/pymolSessions/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessions/;\
        scp -r {pwd}/static/downloads/pymolSessionScripts/* vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessionScripts/;')
    
    subprocess.call(cmd,shell=True)
    
    
if __name__ == '__main__':
    pwd='/home/vivek/Applications/Flask/Kincore'    #location in workhorse
    transfer_to_dunbrack3(pwd)       