# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:01:28 2022

@author: 13237
"""
import ase.io
import numpy as np
from Class import class_interface
import os
import pandas as pd

class Multisurface():
    
    
    def __init__(self):
        
        pass
    
    def renormalizeVector(self,vectorXYZ):
        
        '''
        input: xyz file of g or h vector
        output: normalized vector on each atom and its magnitude
        
        '''
        
        vector = ase.io.read(vectorXYZ).get_positions()
        
        normVector = vector/np.linalg.norm(vector)
        
        return normVector
    
    def orthogonalizeVectors(self,hvector,gvector):
        
        
        hnorm = self.renormalizeVector(hvector)
        gnorm = self.renormalizeVector(gvector)
        natoms,_ = np.shape(hnorm)
        
        hnorm = np.reshape(hnorm, (3*natoms,1)) # reshape hnorm into a column vector
        gnorm = np.reshape(gnorm, (3*natoms,1)) # reshape gnorm into a column vector
        
        vectorspace = np.column_stack((hnorm,gnorm))
        
        ortho,_ = np.linalg.qr(vectorspace)
        
        hortho = ortho[:,0]
        gortho = ortho[:,1]
        
        return hortho,gortho
    
    def jctcorthoVectors(self,hvector,gvector):
        
        hnorm = ase.io.read(hvector).get_positions()
        gnorm = ase.io.read(gvector).get_positions()
        
        natoms,_ = np.shape(hnorm)
        
        hnorm = np.reshape(hnorm, (3*natoms)) # reshape hnorm into a column vector
        gnorm = np.reshape(gnorm, (3*natoms)) # reshape gnorm into a column vector

    
        tan2beta = 2*np.dot(hnorm,gnorm)/(np.dot(gnorm,gnorm)-np.dot(hnorm,hnorm))
        beta = 0.5*np.arctan(tan2beta)
        
        print(beta*180/np.pi)
        print(np.sin(beta))
        hortho = hnorm*np.cos(beta) - gnorm*np.sin(beta)
        gortho = gnorm*np.cos(beta) + hnorm*np.sin(beta)
        
        print(hortho)
        print(np.dot(hortho,gortho)/np.sqrt((np.dot(hortho,hortho)*np.dot(gortho,gortho))))
        
        
        hortho = hortho/np.sqrt(hortho,hortho)
        gortho = hortho/np.sqrt(hortho,hortho)
        
        return hortho,gortho
        
    def branchingInitVel(self,hvector,gvector,kineticE,stepNum):
        
        atomMass = ase.io.read(hvector).get_masses()
        atomSymbol = ase.io.read(hvector).get_chemical_symbols()

        atomNum = np.shape(atomMass)[0]
        massList = np.column_stack((atomMass.T,atomMass.T,atomMass.T))
        massList = np.reshape(massList,(3*atomNum,1))
        
        # print(massList)
        
        
        
        hortho,gortho = self.orthogonalizeVectors(hvector, gvector)
        
        hvectoXYZ = hortho.reshape((atomNum,3))
        gvectoXYZ = gortho.reshape((atomNum,3))
        
        hvectoXYZ = class_interface.Interface().writeXYZFile( pd.DataFrame({'Atom': atomSymbol, 'x':  hvectoXYZ[:,0], 'y':  hvectoXYZ[:,1], 'z':  hvectoXYZ[:,2]}), r'hvector_ortho.xyz')
        gvectoXYZ = class_interface.Interface().writeXYZFile( pd.DataFrame({'Atom': atomSymbol, 'x':  gvectoXYZ[:,0], 'y':  gvectoXYZ[:,1], 'z':  gvectoXYZ[:,2]}), r'gvector_ortho.xyz')
        
        for n,angle in enumerate(np.linspace(0,2*np.pi,stepNum)):
            
            savefolder= r'C:\Users\13237\OneDrive - University of California, Davis\Research Projects\Sesquiterpene\NAMD_setup\n' + str(n) + '\\' 
            # savefolder = r'/Users/zhitao/Desktop/ServerFile/NAMD/n'+ str(n) + '/'
            os.makedirs(savefolder)
            velArr = gortho*np.cos(angle) + hortho*np.sin(angle)
            velScaler = np.sqrt(2*kineticE/(np.sum(massList*velArr*velArr)))
            
            velArr_scaled = velArr*velScaler*0.037/2 # unit transformation from eV to Hartree
            velArr_scaled = velArr_scaled.reshape((atomNum,3))
            np.savetxt(savefolder+'veloc',velArr_scaled,delimiter='  ')
            
            
            
            
            # print(np.sum(0.5*massList*velArr_scaled**2))
            
        
        
        
        
        
    
        
        