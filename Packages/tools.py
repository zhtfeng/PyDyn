#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:08:03 2022

@author: zhitao
"""

import numpy as np
import re

def extractfloat(line,atNum):
    
    '''
    A clumsy function that extracts float numbers from a string
    '''
    l = []
    
    for i in line.split(' '):
        try:
            l.append(np.float(i))
        except ValueError:
            pass
    return l[atNum]

def is_number_check(num):
    
    pattern = re.compile(r'^[-+]?[-0-9]\d*\.\d*|[-+]?\.?[0-9]\d*$')
    result = pattern.match(num)
    if result: return True
    else: return False
    
    