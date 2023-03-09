#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 19:35:28 2019

@author: Zhitao Feng
"""

# !/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import class_interface as interface

RgasK = 0.00198588
avNum = 6.0221415E23
amuToKcal = 4.184E26

class Distribution:

    def __init__(self,shape,rng):
        
        self.shape = shape
        self.rng = rng

    def qmDistribution(self):  # Generate randArrD

        randomArrNormal = self.rng.standard_normal(self.shape) 
        randomArrNormal = randomArrNormal/np.ptp(randomArrNormal)

        return randomArrNormal

    def random_atom_direction_(self, atom_num):  # Generate randArrE

        rand_array_atom = 2 * np.random.randint(2, size=(atom_num, 3)) - 1

        return rand_array_atom


class InitializeMD:

    def __init__(self,software,freqFile, classical, classicalSpacing, temp, initialDis,
                 numimage, timestep,randomstate=None,direction='Only'):
        
        if software in ['Gaussian','gaussian','GAUSSIAN','G16','Gaussian16']: 
            
            linker = interface.GaussianInterface()
        
        molecule = linker.aseFromLogFile(freqFile)
        
        self.symbolArr = molecule.get_chemical_symbols()
        self.massList = molecule.get_masses()
        self.initialGeom = molecule.get_positions()
        
        # frequency from output files (1*Nfreq dimension)
        # reduced mass from output files  (1*Nfreq dimension)
        # Eigenvector of vibrations (Nfreq*Natom*3 dimension)
        # Force constant of vibrational mode  (1*Nfreq dimension)
        
        try:
            self.freq, self.redmass,self.freqmatx,self.forceConst = linker.readFreq(freqFile)
        except TypeError:
            print('Error: Initialization Stopped when reading Frequency calculation output')
            
        self.classical = classical
        self.spacing = classicalSpacing
        self.temp = temp
        self.initialDis = initialDis
        self.numimage = numimage
        self.atomNum = len(self.symbolArr)
        self.timestep = timestep
        self.randomstate = randomstate

    def generate(self,velUnit='ans/s',verbose=False):
        
        self.energyInNormalMode()
            
        if self.randomstate == None:
            
            seed = np.random.randint(2147483647)  # generate a seed so that all random number used for initialization is consistent across all functions
            self.rng = np.random.default_rng(seed=seed)
            
        elif isinstance(self.randomstate,(int,float)):
            
            self.rng = np.random.default_rng(seed = int(self.randomstate)) #user-defined seed for reproduction
                
        self.randomExcitation()
        self.modeDisplacement()
        self.modeVelocity(direction=1) 
        coordinateInitial = self.outputInitialXYZ()
        velocityInitial = self.outputInitialVelXYZ(velUnit)
        
        return coordinateInitial, velocityInitial
    
         
    def energyInNormalMode(self):
        
        '''
        Generate vibrational excitation states of the normal modes base on frequencies
        
        output: vibLevelNum - 1*Nfreq array of integers (vib level of mode)
                zpeinJ - 1*Nfreq array of vib energy in J per molecule
                zpeinKCal - 1*Nfreq array of vib energy in kCal/mol
        
        '''
        h = 6.626075E-34
        c = 29979245800
        avgNum = 6.0221415E23

        if self.classical != 1:

            zpeInJ = self.freq * .5 * h * c

        elif self.classical == 1:

            zpeInJ = self.freq * .5 * h * c * self.spacing

        zpeInKCal = zpeInJ * avgNum / 4184

        if self.temp < 10:

            zpeRat=np.zeros_like(self.freq)
            Q=np.ones_like(self.freq)
            tester=np.ones_like(self.freq)

        else:

            zpeRat = np.exp(-2 * zpeInKCal / (RgasK * self.temp)) - 1e-14
            Q = 1 / (1 - zpeRat)
            tester = 1 / Q
        
        self._Q = Q
        self._tester = tester
        self._zpeRat = zpeRat
        self.zpeInJ = zpeInJ
        self.zpeInKCal = zpeInKCal
        
    def randomExcitation(self):
            
        vibLevelNum = np.zeros_like(self.freq)
        tester = np.copy(self._tester)
        for i in range(len(self._tester)):
            j = 0
            while j <= 4000 * self._zpeRat[i] + 2:

                if self.rng.random(1) > tester[i]:
                    
                    tester[i] += (self._zpeRat[i] ** j) / self._Q[i]
                    vibLevelNum[i] += 1

                j += 1
                
        self.vibLevelNum = vibLevelNum

        return vibLevelNum
        
    def modeDisplacement(self):
        
        '''
        Generate initial displacement along all modes. First we calculate the maximum displacement available. then the displacement along one mode
        depends on the distribution. Summing up all of the displacement gives the ΔGeom, which will be used later
        
        output:
            modeEnInJ - total energy is one vib mode in mDyne Anstrom
            totalModeEnK - float, total energy in vibrations 
            shift - 1* Nfreq array, how many displacements along one vib mode
            geomChange - Natom*3 array, change in geom
            
            
        '''
        vibLevelNum = self.vibLevelNum

        if self.classical == 1:
            
            modeEnInJ = self.zpeInJ * 1e18 * (2 * vibLevelNum)
            modeEnInKCal = self.zpeInKCal * (2 * vibLevelNum)
        else:
            
            modeEnInJ = self.zpeInJ * 1e18 * (2 * vibLevelNum + 1)
            modeEnInKCal = self.zpeInKCal * (2 * vibLevelNum + 1)

        totalModeENK = np.sum(modeEnInKCal)

        for i, each_freq in enumerate(self.freq):

            if each_freq < 10: modeEnInJ[i] = 0 # very small vib modes are anharmonic, remove their zpe

        maxDisplacement = np.sqrt((2 * modeEnInJ / self.forceConst)) # maximum displacement if all kinetic energy are turned into potential energy

        if self.initialDis == 0: # no displacement at all

            shift = np.zeros_like(self.freq)

        elif self.initialDis == 1: # classical displacement

            randArr = self.rng.rand(len(self.freq))
            shift = maxDisplacement * (2 * randArr - 1)

        elif self.initialDis == 2:

            randArrQM = Distribution(len(self.freq),self.rng).qmDistribution()
            shift = maxDisplacement * randArrQM

        for i, each_freq in enumerate(self.freq):

            if each_freq < 10 or i+1 <= self.numimage:
                
                shift[i] = 0

        displaceAlongMode = np.zeros_like(self.freqmatx)
        geomChange = np.zeros((self.atomNum, 3))

        if self.classical != 2: # excluding no displacement condition

            for i in range(len(self.freq)):

                for j in range(self.atomNum):

                    for k in range(3):
                        
                        displaceAlongMode[i, j, k] = self.freqmatx[i, j, k] * shift[i]
                        geomChange[j, k] += displaceAlongMode[i, j, k]
                        
        self.geomshift = shift
        self.modeEninJ = modeEnInJ
        self.geomChange = geomChange

        return modeEnInJ,totalModeENK,shift,geomChange
    

    def modeVelocity(self,direction):
        
        shift = self.geomshift
        modeEn = self.modeEninJ
        imagFreqNum = np.copy(self.numimage)
        
        kinEn = 1e5 * ( modeEn - 0.5*self.forceConst * shift * shift) # kinetic energy is total energy - potantial energy from displacement
        vel = np.sqrt(2 * kinEn*avNum/ (self.redmass))

        if imagFreqNum > 1: imagFreqNum = 1 # accept only one direction for 

        for i in range(len(vel)):

            if i+1 > imagFreqNum and self.rng.random() < 0.5: vel[i] = -vel[i] # The real modes are in random directions

            elif i+1 == imagFreqNum: 

                if direction < 0:
                    
                    vel[i] = -vel[i] # controls positive or negative direction on imag vib mode

        velMode = np.zeros_like(self.freqmatx)
        velArrAtom = np.zeros((self.atomNum, 3))

        if self.classical != 2:

            for i in range(len(self.freq)):

                for j in range(self.atomNum):

                    for k in range(3):
                        
                        velMode[i, j, k] = self.freqmatx[i, j, k] * vel[i] * self.timestep
                        velArrAtom[j, k] += velMode[i, j, k]

        else:

            degFreeEnK = self.temp * RgasK
            degFreeEnJ = degFreeEnK / (avNum / 4184)
            cartEn = degFreeEnJ * 1E18
            kinEnCart = 100000 * cartEn

            for i in range(self.atomNum):

                for j in range(3):
                    
                    velArrAtom[i, j] = (2 * self.rng.choice(2, 1) - 1) * self.timestep * np.sqrt(
                        (2 * kinEnCart / (self.massList[i] / avNum)))
                    
        totalModeKE = np.sum(0.5*self.massList*np.sum(velArrAtom*velArrAtom,axis=1))/(self.timestep**2*amuToKcal) # in kcal/mol
        
        self.velArr = velArrAtom
        return velArrAtom,totalModeKE

    def outputInitialXYZ(self):
        
        outputCoord = self.initialGeom+self.geomChange
        rawxyzDict = {'Atom':self.symbolArr, 'X': outputCoord[:,0], 'Y': outputCoord[:,1],'Z': outputCoord[:,2]}
        dfCoord = pd.DataFrame(rawxyzDict,columns=['Atom','X', 'Y', 'Z'])
        
        return dfCoord
    
    def outputInitialVelXYZ(self,unit='ans/s'):
        
        if unit == 'ans/s':
            
            rawxyzDict = {'Atom':self.symbolArr, 'X': self.velArr[:,0], 'Y': self.velArr[:,1],'Z': self.velArr[:,2]}
            
        if unit == 'amu':
            
            amuVel = self.velArr/2.18769126364e16
            rawxyzDict = {'Atom':self.symbolArr, 'X': amuVel[:,0], 'Y': amuVel[:,1],'Z': amuVel[:,2]}
        
        dfVel = pd.DataFrame(rawxyzDict,columns=['Atom','X', 'Y', 'Z'])
        
        return dfVel
    

    
    
         
        
    # def writeFingerprint(self):
        
    #     outputFile = open('Vibexcitation', 'w')
    #     outputFile.writelines(self.vibLevelNum)
        
        


















