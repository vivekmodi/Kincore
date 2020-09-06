#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:37:51 2020

@author: vivekmodi
"""
from datetime import datetime

def create_datefile(pwd):
    today=str(datetime.now())[0:10].strip()
    fhandle_update_date=open(f'{pwd}/update-date.txt','w')
    fhandle_update_date.write(f'{today}')
    fhandle_update_date.close()