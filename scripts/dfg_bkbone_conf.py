#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math,sys

def dihedral_labels(pwd,dfg_label,xdfg_phi,xdfg_psi,dfg_asp_phi,dfg_asp_psi,dfg_phe_phi,dfg_phe_psi,dfg_phe_chi1):
    

    xdfg_phi=float(xdfg_phi);xdfg_psi=float(xdfg_psi);dfg_asp_phi=float(dfg_asp_phi);dfg_asp_psi=float(dfg_asp_psi);
    dfg_phe_phi=float(dfg_phe_phi);dfg_phe_psi=float(dfg_phe_psi);dfg_phe_chi1=float(dfg_phe_chi1)
    
    if dfg_label=='None':
        return 'None'
    if xdfg_phi==999 or xdfg_psi==999 or dfg_asp_phi==999 or dfg_asp_psi==999 or dfg_phe_phi==999 or  dfg_phe_psi==999 or dfg_phe_chi1==999:
        return 'None'

    fhandle_dfgin=open((pwd+'/server/Cluster-centroids-DFGin-filtered-2.25'),'r')
    fhandle_dfginter=open((pwd+'/server/Cluster-centroids-DFGinter-filtered-2.25'),'r')
    fhandle_dfgout=open((pwd+'/server/Cluster-centroids-DFGout-filtered-2.25'),'r')

    if dfg_label=='DFGin':
        min_cluster=0
        for lines in fhandle_dfgin:
            lines=lines.strip("\n");lines=lines.split(" ")

            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(lines[2]))))+(1-math.cos(math.radians(xdfg_psi-float(lines[3]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))

            if cosine_dis<0.3:
                min_cluster=int(lines[1])
        fhandle_dfgin.close()
        if min_cluster==0:
            return 'None'
        if min_cluster==1:
            return 'BLAminus'
        if min_cluster==2:
            return 'BLAplus'
        if min_cluster==3:
            return 'ABAminus'
        if min_cluster==4:
            return 'BLBminus'
        if min_cluster==5:
            return 'BLBplus'
        if min_cluster==6:
            return 'BLBtrans'

    if dfg_label=='DFGinter':
        min_cluster=0
        for lines in fhandle_dfginter:
            lines=lines.strip("\n");lines=lines.split(" ")
            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(lines[2]))))+(1-math.cos(math.radians(xdfg_psi-float(lines[3]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))
            fhandle_dfginter.close()
            if cosine_dis<0.3:
                return 'BABtrans'
            else:
                return 'None'

    if dfg_label=='DFGout':
        min_cluster=0
        for lines in fhandle_dfgout:
            lines=lines.strip("\n");lines=lines.split(" ")
            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(lines[2]))))+(1-math.cos(math.radians(xdfg_psi-float(lines[3]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(lines[4]))))+(1-math.cos(math.radians(dfg_asp_psi-float(lines[5]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(lines[6]))))+(1-math.cos(math.radians(dfg_phe_psi-float(lines[7]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(lines[8])))))
            fhandle_dfgout.close()
            if cosine_dis<0.3:
                return 'BBAminus'
            else:
                return 'None'
    fhandle_dfgin.close()
    fhandle_dfgout.close()
    fhandle_dfginter.close()

if __name__ == '__main__':
    dfg_label=sys.argv[1];xdfg_phi=sys.argv[2];xdfg_psi=sys.argv[3];dfg_asp_phi=sys.argv[4];dfg_asp_psi=sys.argv[5];dfg_phe_phi=sys.argv[6];dfg_phe_psi=sys.argv[7];dfg_phe_chi1=sys.argv[8]
    dfg_bkbone=dihedral_labels(dfg_label,xdfg_phi,xdfg_psi,dfg_asp_phi,dfg_asp_psi,dfg_phe_phi,dfg_phe_psi,dfg_phe_chi1)
    print(dfg_bkbone)
