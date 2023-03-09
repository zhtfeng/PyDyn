#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:19:05 2022

@author: zhitao
"""

from Class import class_interface


filename = r'/Users/zhitao/OneDrive - University of California, Davis/Research Projects/Norrish-Cope/Scan/CNbond-scan/n0/MOLCAS.log'
readmolcas = class_interface.MolcasInterface().readMCPDFTlog(filename)