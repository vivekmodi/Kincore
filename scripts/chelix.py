#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
#homedir='/home/vivekmodi/Applications/Flask/Kinases/'

def chelix_conformation(dis_sb):
    if dis_sb<=10:
        return 'Chelix-in'
    else:
        return 'Chelix-out'

if __name__ == '__main__':
    dis_sb=int(sys.argv[1])
    chelix_conf=chelix_conformation(dis_sb)
    print(chelix_conf)
