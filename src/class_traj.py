#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:38:12 2022

@author: zhitao
"""
import re
import numpy as np
import ase.io
import pandas as pd
from ase.units import Bohr
from ase import Atoms
from Class import tools
# import class_interface as interface

class Traj():
    
    '''
    trajectory class on top of reading xyz file
    
    '''
    
    def __init__(self,xyzFilename):
        
        '''
        The xyz file provided will be converted into a list of Atom class self.xyz
        
        '''
        self.xyzFile = xyzFilename
        self.xyz = ase.io.read(xyzFilename,index=':')
        self.atmNumber = len(self.xyz[0].get_chemical_symbols()) # look at the first frame of the traj and figure out the number of atoms
        self.name = None
    def getCommentLine(self):
        '''
        
        There might me important info in the comment line, such energy or state info of NAMD, we can pull out such info and store commentline in a list

        Returns
        -------
        commentLineList : TYPE
            a list of all comment lines

        '''
        frameNum = len(self.xyz)
        
        commentLineIndex = [(self.atmNumber+2)*x+1 for x in range(frameNum)]
        commentLineList = []
        
        with open(self.xyzFile) as f:
            
            fileList = f.readlines()
            
            for i in commentLineIndex:
                
                commentLineList.append(fileList[i])
                
        return commentLineList
        
    def writexyzBohr(self):
        
        '''
        convert original xyz into a list of pandas dataframes, with units converted to Bohr for molden file
        
        '''
        
        DFlist = []
        
        for frame in self.xyz:
            
            molinBohr = Atoms(frame.get_chemical_symbols(),frame.get_positions()/Bohr)
            DFlist.append(molinBohr)
            self.xyzinDFList = DFlist
        
        return self.xyzinDFList
    
    def getDistanceArr(self,a0,a1):
        
        distArr = []
        
        for frame in self.xyz:
            
            mol = Atoms(frame.get_chemical_symbols(),frame.get_positions())
            distArr.append(mol.get_distance(a0,a1))
        
        return np.array(distArr)
    
    def getAngleArr(self,a0,a1,a2):
        
        angleArr = []
        
        for frame in self.xyz:
            
            mol = Atoms(frame.get_chemical_symbols(),frame.get_positions())
            angleArr.append(mol.get_angle(a0,a1,a2))
        
        return np.array(angleArr)
        
    def getDihedralArr(self,a0,a1,a2,a3):
        
        dihArr = []
        
        for frame in self.xyz:
            
            mol = Atoms(frame.get_chemical_symbols(),frame.get_positions())
            dihArr.append(mol.get_dihedral(a0,a1,a2,a3))
            
        # print(dihArr)
        
        return np.array(dihArr)*np.pi/180
    
    def getFrameXYZ(self,atStep):
        
        return self.xyz[atStep]
    
class NAMDTraj(Traj):
    
    def getHoppingPoint(self):
        
        commentLineList = self.getCommentLine()
        atStateArr = np.zeros((len(commentLineList)))
        hopArr = np.zeros_like(atStateArr)
        
        for i,j in enumerate(commentLineList):
            
            atState = tools.extractfloat(j, 1)
            atStateArr[i] = atState
            
            if i != 0:
                
                hopArr[i] = atStateArr[i] - atStateArr[i-1]
    
        hopPt = np.where(hopArr<0)[-1]
        if len(hopPt) > 0 :
            
            return hopPt[0]
        
        else:
            
            return None
        
        
            
            
        
    

    
    
    

    
    
        
        
    
    
        
    
            
    

    
    
    