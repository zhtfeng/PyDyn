# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 17:18:19 2021

@author: Feng
"""

import os
import pandas as pd
import numpy as np

class Multiwfn():
    
    def __init__(self,multiwfnPath=None):
        
        self.path = multiwfnPath
    
    def execute(self,multinput,chkfile):
        
        command = self.path + ' ' + str(chkfile) + ' -isilent 1 < ' + str(multinput) + ' > ' + str(chkfile) + '.out'
        print(command)
        os.system(command)
    
    def extractfloat(self,line,atNum):
        
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
    
    
    def parseESP(self, outfile,smilesFile):
        
        espDict = {}
        chargeList = []
        atomSizeList = []
        odiList = []
        posSurfaceIndicator=False
        negSurfaceIndicator=False
        moleculeIndex = outfile[:-9]
        print(moleculeIndex)
        smilesDF = pd.read_excel(smilesFile)
        smilesList = smilesDF[smilesDF['index'] == moleculeIndex].smiles.tolist()
        if len(smilesList)> 0 :espDict['smiles'] = smilesList[0]
        with open(outfile,'r') as file:
            
            for line in list(file):
                
                if 'Volume:' in line and 'Charge:' not in line:espDict['Volume'] = self.extractfloat(line, 1)
                if ' Minimal value:' in line and 'eV' not in line:
                    espDict['minValue'] = self.extractfloat(line, 0)
                    espDict['maxValue'] = self.extractfloat(line, 1)
                if ' Positive surface area:' in line and not posSurfaceIndicator: 
                    posSurfaceIndicator = True
                    espDict['posSurArea'] = self.extractfloat(line, 1)
                if ' Negative surface area:' in line and not negSurfaceIndicator: 
                    negSurfaceIndicator = True
                    espDict['negSurArea'] = self.extractfloat(line, 1)
                if ' Overall average value:' in line:espDict['ovrAvgVal'] = self.extractfloat(line, 1)
                if ' Positive average value:' in line:espDict['posAvgVal'] = self.extractfloat(line, 1)
                if ' Negative average value:' in line:espDict['negAvgVal'] = self.extractfloat(line, 1)
                if ' Positive variance:' in line:espDict['posVarVal'] = self.extractfloat(line, 1)
                if ' Negative variance:' in line:espDict['negVarVal'] = self.extractfloat(line, 1)
                if ' Balance of charges (nu):' in line:espDict['BoC'] = self.extractfloat(line, 0)
                if 'Internal charge separation (Pi):' in line:espDict['ICS'] = self.extractfloat(line, 1)
                if 'Molecular polarity index (MPI):' in line:espDict['MPI'] = self.extractfloat(line, 1)
                if 'Nonpolar surface area' in line:espDict['nonPolSuf'] = self.extractfloat(line, 1)
                if 'Polar surface area' in line:espDict['PolSurf'] = self.extractfloat(line, 1)
                if 'HOMO, energy:' in line:espDict['homo'] = self.extractfloat(line, 2)
                if 'LUMO, energy:' in line:espDict['lumo'] = self.extractfloat(line, 2)
                if 'HOMO-LUMO gap:' in line:espDict['gap'] = self.extractfloat(line, 1)
                if 'Magnitude of dipole moment:' in line: espDict['dipole'] = self.extractfloat(line, 1)
                if 'Magnitude of the traceless quadrupole moment tensor:' in line:espDict['quadrupole']=self.extractfloat(line, 0)
                if 'Magnitude: |Q_3|=' in line:
                    espDict['octopole'] =self.extractfloat(line, 0)
                    
                if 'Corrected charge:' in line:
                    
                    chargeList.append(self.extractfloat(line, 0))
                
                if 'Volume:' in line and 'Charge:' in line:
                    
                    atomSizeList.append(self.extractfloat(line, 2))
                
                if 'Minimal value:' in line and 'eV,' in line:
                    
                    espDict['alieminValue'] = self.extractfloat(line, 0)
                    espDict['aliemaxValue'] = self.extractfloat(line, 1)
                    
                if 'Average value:' in line and 'eV,' in line: espDict['alieAvg'] = self.extractfloat(line, 1)
                if 'Variance'in line and 'eV^2,' in line: espDict['alieVar'] = self.extractfloat(line, 1)
                
                if 'Orbital delocalization index:' in line:
                    
                    odiList.append(self.extractfloat(line, 0))
                    
        
        espDict['odiHomo'] = odiList[0]
        espDict['odiLumo'] = odiList[1]    
        chargeArr =np.array(chargeList)
        posmax = np.argsort(chargeArr, axis=0)[-1]
        possecmax = np.argsort(chargeArr, axis=0)[-2]
        negmin = np.argsort(chargeArr, axis=0)[0]
        negsecmin = np.argsort(chargeArr, axis=0)[1]
        
        espDict['chrg1'] = chargeArr[posmax]
        espDict['chrg2'] = chargeArr[possecmax]
        espDict['chrg-1'] = chargeArr[negmin]
        espDict['chrg-2'] = chargeArr[negsecmin]
        
        espDict['size1'] = atomSizeList[posmax]
        espDict['size2'] = atomSizeList[possecmax]
        espDict['size-1'] = atomSizeList[negmin]
        espDict['size-2'] = atomSizeList[negsecmin]
        
        
                         
        return espDict
    
    
        

        
        
        
    
    
