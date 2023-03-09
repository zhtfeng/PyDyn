# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 10:09:10 2021

@author: Feng
"""

import numpy as np
import pandas as pd
import class_interface as interface
import class_mdlib as mdlib
import class_config as config
import os,sys

class MD:
    
    def __init__(self):
        
        pass
    
    def orcaMD(self,initcoord,initvel,filename):
        
        orcawrite = interface.ORCAInterface()
        orcawrite.writeMDRestartFile(initcoord,initvel,filename)
        orcawrite.writeXYZFile(initcoord,filename)
        
    def orcaMDReverse(self,initcoord,initvel,filename):
        
        initvel['X'] = -initvel['X']
        initvel['Y'] = -initvel['Y']
        initvel['Z'] = -initvel['Z']
        self.orcaMD(initcoord, initvel, filename)
        
    def orcaMDEnsembleFromGaussian(self,trajNum,gaussianOutput):
        
        setup = config.Configuration("config.txt",'out.txt',None)
        setup.read_config()

        initCond = mdlib.InitializeMD(software='Gaussian',freqFile=gaussianOutput, \
                          classical = setup.config['classical'],\
                          classicalSpacing = setup.config['classicalSpacing'],\
                          temp = setup.config['temperature'],\
                          initialDis = setup.config['Initialdis'],\
                          numimage = setup.config['numimage'],\
                          timestep = setup.config['timestep'],
                          randomstate=None)
        
        for i in range(trajNum):
            
            initcoord,initvel = initCond.generate(velUnit='amu') 
            filename = 'traj' + str(i)
            currentDir = os.getcwd()
            print(os.getcwd())
#            os.mkdir('filename')
#            os.chdir('filename')
#            os.mkdir('forward')
#            
#            os.chdir('forward')
#            
#            self.orcaMD(initcoord,initvel,filename)
#            
#            self.orcaMDReverse(initcoord, initvel, filename+'rev')
            
            
            
            
            
            
            
            
            
    
    

