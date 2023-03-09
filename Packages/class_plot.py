# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 18:23:25 2022

@author: zhitao feng @ uc davis
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import pandas as pd
import numpy as np
from matplotlib import cm
mpl.rc('font',family='Times New Roman')
plt.rcParams.update({'font.size': 16})

class Plot():
    
    def __init__(self):
        
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = 'Times New Roman'


class PlotSurface(Plot):
    
    def branchingSurface(self,surfaceData):
        
        fig = plt.figure(1, figsize=(20, 20),dpi=150)
        ax = Axes3D(fig)
        ax.set_xlabel('H Vector (NAC)',labelpad=20)
        ax.set_ylabel('G Vector (Gradient Difference)',labelpad=20)
        ax.set_zlabel('Electronic Energy (kcal/mol)',labelpad=20)
        
        df_surfaces = pd.read_excel(surfaceData, header=None)
        
        surface1 = df_surfaces.iloc[:,0]
        surface2 = df_surfaces.iloc[:,1]
        
        zeroPoint = np.min(surface2)
        print(zeroPoint)
        surface1 = (np.array(surface1)-zeroPoint)*627.5
        surface2 = (np.array(surface2)-zeroPoint)*627.5
        
        dim = int(np.sqrt(surface1.shape)[0])
        
        surface1 = surface1.reshape([dim,dim])
        surface2 = surface2.reshape([dim,dim])
        
        X = np.arange(-0.6,0.7,0.1)
        Y = np.arange(-0.6,0.7,0.1)
        imatx = np.zeros_like(surface1)
        jmatx = np.zeros_like(surface1)
        
        for i in range(0,len(X)):
            for j in range(0,len(Y)):
                
                imatx[i,j] = X[i]
                jmatx[i,j] = Y[j]
        
 
        surf = ax.plot_surface(imatx, jmatx, surface1,
                       linewidth=0, antialiased=False,cmap=cm.jet,alpha=0.9)
        surf = ax.plot_surface(imatx, jmatx, surface2,
                       linewidth=0, antialiased=False,cmap=cm.jet,alpha=0.75)
        
        ax.view_init(30, 10)
        
        # print(surface1)
        
class PolarPlot(Plot):
    
    
    def polarTraj(self,trajfile):
        
        df = pd.read_excel(trajfile)
        fig = plt.figure()
        ax = fig.add_subplot(projection='polar')
        thetaArr = df['angle']*2*np.pi/360
        energyArr = df['energy']*23
        
        scatter = ax.scatter(thetaArr,energyArr,s=60,c=df['color'])
        ax.set_ylim(0,30)
        

        print(thetaArr)
        
        
class TrajEnsemblePlot(Plot):
    
    
    def plotTraj1DTimePolar(self,ax,geomX):
        
        plot = ax.plot(np.pi*geomX/180,range(len(geomX)))

    
    def plotTraj3DPolar(self,ax,geomX,geomY,geomZ,color):
        
        plot = ax.plot(geomX,geomY,geomZ,c=color, linewidth=1.0)        
        
        
        