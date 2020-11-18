#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:40:32 2020

@author: vivekmodi
"""
import subprocess

def transfer_to_dunbrack3(pwd):
    cmd=(f'rsync -avr --progress {pwd}/formattedSeq/* vivek@dunbrack3.fccc.edu:/var/www/html/site/formattedSeq/;\
        rsync -avr --progress {pwd}/kinasechains/* vivek@dunbrack3.fccc.edu:/var/www/html/site/kinasechains/;\
        rsync -avr --progress {pwd}/kinasechains_dihedrals/* vivek@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_dihedrals/;\
        rsync -avr --progress {pwd}/kinasechains_renumber_alignment/* vivek@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_renumber_alignment/;\
        rsync -avr --progress {pwd}/kinasechains_renumber_uniprot/* vivek@dunbrack3.fccc.edu:/var/www/html/site/kinasechains_renumber_uniprot/;\
        rsync -avr --progress {pwd}/kinasecifs/* vivek@dunbrack3.fccc.edu:/var/www/html/site/kinasecifs/;\
        rsync -avr --progress {pwd}/static/kinasechainsNGL/* vivek@dunbrack3.fccc.edu:/var/www/html/site/static/kinasechainsNGL;\
        scp -r {pwd}/static/downloads/coordinateFiles/* vivek@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/coordinateFiles/;\
        scp -r {pwd}/static/downloads/pymolSessions/* vivek@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessions/;\
        scp -r {pwd}/static/downloads/pymolSessionScripts/* vivek@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessionScripts/;')
    
    subprocess.call(cmd,shell=True)