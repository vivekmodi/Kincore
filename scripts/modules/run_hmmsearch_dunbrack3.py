#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:24:37 2020

@author: vivekmodi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 08:29:49 2020

@author: vivekmodi
"""

#Rename to run_hmmsearch.py after moving to server

import subprocess

def run_hmmsearch(pwd,pdbfilename,index,conf_df):
    model_id=conf_df.at[index,'Model_id']
    chain_id=conf_df.at[index,'Chain_id']
    
    cmd=(f'/var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_AGC.hmmer.txt {pwd}/server/AGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_CAMK.hmmer.txt {pwd}/server/CAMK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_CK1.hmmer.txt {pwd}/server/CK1.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_CMGC.hmmer.txt {pwd}/server/CMGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_NEK.hmmer.txt {pwd}/server/NEK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_RGC.hmmer.txt {pwd}/server/RGC.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_STE.hmmer.txt {pwd}/server/STE.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_TKL.hmmer.txt {pwd}/server/TKL.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_TYR.hmmer.txt {pwd}/server/TYR.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_HASP.hmmer.txt {pwd}/server/HASP.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_WNK.hmmer.txt {pwd}/server/WNK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_BUB.hmmer.txt {pwd}/server/BUB.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;\
         /var/www/html/site/binary/hmmer-installed/bin/hmmsearch -o {pwd}/server/{pdbfilename[0:-4]}_{chain_id}_ULK.hmmer.txt {pwd}/server/ULK.hmm {pwd}/server/uploads/{pdbfilename[0:-4]}_{chain_id}.fasta;')
    
    process=subprocess.Popen(cmd,shell=True)
    process.communicate()
    process.wait()
    #subprocess.call(cmd,shell=True)
    return
