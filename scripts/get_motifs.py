#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:33:16 2020

@author: vivekmodi
"""

def get_motifs(pwd):
    
    gene_domain='';gene=dict();group=dict();uniprot=dict();uniacc=dict();domainstart=dict();domainend=dict();alk_res=dict();alk_num=dict();rre_res=dict();rre_num=dict();
    lrl_res=dict();lrl_num=dict();hrd_res=dict();hrd_num=dict();dfg_res=dict();dfg_num=dict();ape_res=dict();ape_num=dict();dlw_res=dict();dlw_num=dict();isr_res=dict();isr_num=dict();
    
    fhandle_motifs=open((pwd+"/motifs.csv"),"r")
    fhandle_motifs.readline()       #Skip the header
    for lines in fhandle_motifs:
        lines=lines.strip("\n");lines=lines.split(",")
        gene_domain=lines[1]
        gene[gene_domain]=lines[0]
        group[gene_domain]=lines[2]
        uniprot[gene_domain]=lines[3]
        uniacc[gene_domain]=lines[4]
        domainstart[gene_domain]=lines[6]
        domainend[gene_domain]=lines[7]
        alk_res[gene_domain]=lines[8]
        alk_num[gene_domain]=lines[9]
        rre_res[gene_domain]=lines[10]
        rre_num[gene_domain]=lines[11]
        lrl_res[gene_domain]=lines[12]
        lrl_num[gene_domain]=lines[13]
        hrd_res[gene_domain]=lines[14]
        hrd_num[gene_domain]=lines[15]
        dfg_res[gene_domain]=lines[16]
        dfg_num[gene_domain]=lines[17]
        ape_res[gene_domain]=lines[18]
        ape_num[gene_domain]=lines[19]
        dlw_res[gene_domain]=lines[20]
        dlw_num[gene_domain]=lines[21]
        isr_res[gene_domain]=lines[22]
        isr_num[gene_domain]=lines[23]
        #if uni_name=='TITIN_HUMAN':
        #    alk_lys=alk_lys-32154;rre_glu=rre_glu-32154;lrl_l=lrl_l-32154;hrd_asp=hrd_asp-32154;ape_glu=ape_glu-32154;dlw_trp=dlw_trp-32154;isr_arg=isr_arg-32154
    return(gene_domain,gene,group,uniprot,uniacc,domainstart,domainend,alk_res,alk_num,rre_res,rre_num,lrl_res,lrl_num,hrd_res,hrd_num,dfg_res,dfg_num,ape_res,ape_num,\
           dlw_res,dlw_num,isr_res,isr_num)