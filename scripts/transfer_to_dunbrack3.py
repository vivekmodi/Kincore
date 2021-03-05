#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:40:32 2020

@author: vivekmodi
"""
import subprocess

def transfer_to_dunbrack3(pwd):
    cmd=(f'rsync -avr --progress {pwd}/formattedSeq/ vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/formattedSeq/;\
        rsync -avr --progress {pwd}/static/kinasechainsNGL/ vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/kinasechainsNGL;\
        rsync -avr --progress {pwd}/static/downloads/coordinateFiles/ vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/coordinateFiles/;\
        rsync -avr --progress {pwd}/static/downloads/pymolSessions/ vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessions/;\
        rsync -avr --progress {pwd}/static/downloads/pymolSessionScripts/ vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/pymolSessionScripts/;\
        rsync -avr --progress {pwd}/update-date.txt vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/;\
        rsync -avr --progress {pwd}/data_new.sqlite vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/;\
        rsync -avr --progress {pwd}/static/downloads/fasta_with_labels/PK_labels_PDB.fasta vivekmodi@dunbrack3.fccc.edu:/var/www/html/site/static/downloads/fasta_with_labels/PK_labels_PDB.fasta')
    
    subprocess.call(cmd,shell=True)
    
    
if __name__ == '__main__':
    pwd='/home/vivek/Applications/Flask/Kincore'    #location in workhorse
    transfer_to_dunbrack3(pwd)       