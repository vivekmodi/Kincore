#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def spatial_label(dis_phe_lys,dis_phe_rre4):
    if dis_phe_rre4<=11 and dis_phe_lys>=11 and dis_phe_rre4!=999 and dis_phe_lys!=999:
        return 'DFGin'
    elif dis_phe_rre4>11 and dis_phe_lys<=14 and dis_phe_rre4!=999 and dis_phe_lys!=999:
        return 'DFGout'
    elif dis_phe_rre4<=11 and dis_phe_lys<=11 and dis_phe_rre4!=999 and dis_phe_lys!=999:
        return 'DFGinter'
    else:
        return 'None'

if __name__ == '__main__':
    dis_phe_lys=int(sys.argv[1]); dis_phe_rre4=int(sys.argv[2])
    dfg_label=spatial_label(dis_phe_lys,dis_phe_rre4)
    print(dfg_label)
