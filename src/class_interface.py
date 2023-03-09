# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 22:21:39 2021

@author: Zhitao Feng
"""
import os
import re
import numpy as np
import pandas as pd
import ase.io
import json
import shutil

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


class Interface():
                
    def writeXYZ(self, content, filename):
        
        with open(filename, 'w') as f:
            
            f.writelines(content)
        f.close()
    
    def writeXYZFile(self,coord,filename):
               
        with open(filename+'.xyz', 'a') as file:
            
            file.write(str(len(coord))+'\n')
            file.write('XYZ generated by Pydyn\n')
            file.writelines(coord.to_string(header=False, index=False)) 
            
    def writeXYZFile_orcaScan(self,coord,filename):
               
        with open(filename, 'a') as file:
            
            file.write(str(len(coord))+'\n')
            file.write('XYZ generated by Pydyn\n')
            file.writelines(coord.to_string(header=False, index=False)) 
            # file.writelines('\n')
            file.write('\n>\n')
        
    def xyztoASE(self,xyzFilename):
        
        self.mol = ase.io.read(xyzFilename)

        return self.mol
    
    def coordList_to_xyz(self,filename):
        
        with open(filename, 'w') as file:
            for i,coord in enumerate(self.coordList):
                
                file.write(str(len(coord))+'\n')
                file.write('XYZ generated by Pydyn\n')
                file.writelines(coord.to_string(header=False, index=False)) 
                file.write('\n')     
        
class ORCAInterface(Interface):
    
    def writeMDRestartFile(self,initcoord,initvel,filename):
        
        mdRestartFilename = str(filename) + '.mdrestart'
        atomNum = np.shape(initcoord)[0]
        forceDF = pd.DataFrame({'Atom': list(initcoord['Atom']), 'x': np.zeros(atomNum), 'y': np.zeros(atomNum), 'z': np.zeros(atomNum)})
               
        with open(mdRestartFilename,'w') as mdFile:
            
            mdFile.write('&AtomCount\n')
            mdFile.write(str(atomNum)+'\n')
            mdFile.write('&CurrentStep\n')
            mdFile.write('0\n')
            mdFile.write('&SimulationTime\n')
            mdFile.write('0.00\n')
            mdFile.write('&Positions\n')
            mdFile.writelines(initcoord.to_string(header=False, index=False))
            mdFile.write('\n&Velocities\n')
            mdFile.writelines(initvel.to_string(header=False, index=False))
            mdFile.write('\n&Forces\n')
            mdFile.writelines(forceDF.to_string(header=False, index=False))
            mdFile.write('\n')
            
    def writeORCAInput(self,coord,prefixName,inputName,destFolder):
        
        prefix = open(prefixName).readlines()
        
        with open(destFolder + str(inputName)+'.inp','w') as inputfile:
            
            inputfile.writelines(prefix)
            inputfile.writelines(coord.to_string(header=False, index=False))
            inputfile.writelines(['\n*\n\n']) 
            
    def readDATFileEnergy(self,datFilename):
        '''
            Parameters
        ----------
        datFilename : TYPE
            DESCRIPTION.

        Returns
        -------
        Numpy Arrary of energy

        '''
        energyList = []
        
        with open(datFilename) as file:
            
            for line in file:
                
                energyList.append(extractfloat(line, -1))
        
        return np.array(energyList)
    
    def convertAllXYZToXYZ(self,allxyzFilename,xyzFilename):
        
        with open(allxyzFilename) as file:
            
            xyzdata = []
            allxyzdata = file.readlines()
            
            for i in allxyzdata:
               
               if '>' not in i: xyzdata.append(i)
               
        with open(xyzFilename,'w') as file:
            
            file.writelines(xyzdata)
            
                     
            
            
        
            
class GaussianInterface(Interface):
    
    def readGeom(self,filename):
                    
        with open(filename, 'r') as file:
                      
            # First by reversing the whole file, we look for the last 'standard orientation'
            
            for lineNum,line in enumerate(reversed(list(file))):
                
                if ' Standard orientation:' in line:
                    
                    coordLineNum = -lineNum - 1 # it is in reverse order, we need negetive sign
                    break
            
            # Second we look for the symbolic Z-matric line, which will offer us the atomic symbols and number of atoms
        self.coord = []     
        with open(filename, 'r') as file:    
            flag = False
            for line in file:
 
                if ' Symbolic Z-matrix:' in line:
                    
                    flag = True
                
                if flag :

                    if len(line) < 3:break
                    self.coord.append(line)
                    
        self.coord = [i[:3] for i in self.coord[2:]]
        atomNum = len(self.coord)
        
        with open(filename, 'r') as file: 
            
            xyzList = list(file)[coordLineNum+5:coordLineNum+atomNum+5]
            
        xyzList = [list(map(float, j.split())) for i,j in enumerate(xyzList)]
        xyzList = np.array([i[3:] for i in xyzList]).T
        self.coord = pd.DataFrame({'Atom': self.coord, 'x': xyzList[0], 'y': xyzList[1], 'z': xyzList[2]})
        
        return self.coord
    
    def readIRCLQA(self,filename):
        
        atomList = []     
        with open(filename, 'r') as file:    
            
            flag = False
            for line in file:
        
                if ' Symbolic Z-matrix:' in line:
                    
                    flag = True
                
                if flag :
        
                    if len(line) < 3:break
                    atomList.append(line)
            
                    
        atomList = [i[:3] for i in atomList[2:]]
        atomNum = len(atomList)
        forwardcoordList = []
        reversecoordList = []
        forwardFlag = True
        with open(filename, 'r') as file:
            file = list(file)
            for lineNum,line in enumerate(file):
                
                if 'Calculation of FORWARD path complete' in line:
                    
                    forwardFlag = False

                
                if 'Cartesian Coordinates (Ang):' in line:
                    
                    xyz = list(file)[lineNum+5:lineNum+atomNum+5]
                    xyz = [list(map(float, j.split())) for i,j in enumerate(xyz)]
                    xyz = np.array([i[2:] for i in xyz]).T
                    coord = pd.DataFrame({'Atom': atomList, 'x': xyz[0], 'y': xyz[1], 'z': xyz[2]})
                    
                    if forwardFlag: 
                        forwardcoordList.append(coord)

                    else: 
                        reversecoordList.append(coord)

        reversecoordList.reverse()
        self.coordList = reversecoordList+ forwardcoordList 
                
    def readScan(self,filename):
        
        pass
        

    def readFreq(self,filename):

    # Input filename

    # Output: freq - list of frequencies of all the modes
    # reduced mass - list of reduced mass of each mode
    # freq_matx - a [num_of_mode,num_of_atom,3] 3D array that contains the cartesian info
    # of all the vibration modes
    # force_list - list of all the force constants

        line_num_list = []
        freq_list, reduced_mass_list, force_list = [], [], []
        freq_file = open(filename)
        freq_lines = freq_file.readlines()
        matx_disorder = []

        for line_num, each_line in enumerate(freq_lines):

            re.split(r'[\s]', each_line)
            #           print(each_line)
            if 'Frequencies' in each_line:
                line_num_list.append(line_num)

        atom_num = int((line_num_list[1] - line_num_list[0] - 7) / 3)
        line_num_list = line_num_list[:(len(line_num_list) - atom_num + 2)]
        freq_matx = np.zeros((3 * atom_num - 6, atom_num, 3))
        raw_matx = np.zeros((3 * atom_num, 3 * atom_num - 6))
        which_freq = 0

        for i, each_line_num in enumerate(line_num_list):

            freqline = re.findall("-?\d+\.*\d*", freq_lines[each_line_num])
            reduced_mass_line = re.findall("-?\d+\.*\d*", freq_lines[each_line_num + 1])
            force_line = re.findall("-?\d+\.*\d*", freq_lines[each_line_num + 2])
            freq_list += list(map(float, freqline))
            reduced_mass_list += list(map(float, reduced_mass_line))
            force_list += list(map(float, force_line))

            for atom_linenum, matx_lines in enumerate(range(each_line_num + 5, each_line_num + 5 + 3 * atom_num)):

                matx_disorder = list(map(float, re.findall("-?\d+\.*\d*", freq_lines[matx_lines])[3:]))

                for item in range(len(matx_disorder)):
                    try:
                        raw_matx[atom_linenum, item + which_freq] = matx_disorder[item]
                    except IndexError:
                        
                        print('An error occured when reading vibrational mode eigenvectors, Please check whether freq=hpmodes is in your input')
                        return None
                    
            which_freq += len(matx_disorder)

        row_num, col_num = np.shape(raw_matx)

        freq_list = np.array(freq_list)
        for cols in range(col_num):

            for rows in range(row_num):

                freq_matx[cols, int(rows / 3), rows % 3] = raw_matx[rows, cols]

        freq_modified = np.zeros_like(freq_list)
        for i, each_freq in enumerate(freq_list):

            if each_freq < 0:

                freq_modified[i] = 2

            else:

                freq_modified[i] = each_freq
        freq_file.close()

        return freq_modified, np.array(reduced_mass_list), freq_matx, np.array(force_list)
    
    def LogFileToXYZ(self,inputFilename,outputFilename):
        
        self.readGeom(inputFilename)
        self.writeXYZFile(self.coord,outputFilename)
        
    def aseFromLogFile(self,inputFilename):
        
        self.LogFileToXYZ(inputFilename, 'temp.xyz')
        molecule = ase.io.read('temp.xyz')

        return molecule
    
    def outputElecEnergy(self,filename):
        
        with open(filename, 'r') as file:
                      
            # First by reversing the whole file, we look for the last 'SCF Done'
            
            for lineNum,line in enumerate(reversed(list(file))):
                
                if 'SCF Done:' in line:
                    
                    scfLineNum = -lineNum - 1 # it is in reverse order, we need a negetive sign
                    break
                
                
class MolcasInterface(Interface):
    
    def __init__(self,atomNum = 0,filename=None):
    
        self.outfile = filename
        # self.iterList = self.IterNum()
        self.atomNum = atomNum
    def writeBranchingSpace_to_allXYZ(self,geomXYZ,hvectorXYZ,gvectorXYZ,stepsize,ptNumberDirection):
        
        
        initcoord = self.xyztoASE(geomXYZ).get_positions()
        hvector = self.xyztoASE(hvectorXYZ).get_positions()
        gvector = self.xyztoASE(gvectorXYZ).get_positions()
        openfile = open('neb.xyz','w')
            
                    
        originPoint = initcoord - hvector*ptNumberDirection*stepsize-gvector*ptNumberDirection*stepsize
        pt = 0
        for i in np.arange(0,2*ptNumberDirection+1):
            
            for j in np.arange(0,2*ptNumberDirection+1):
                
                currentCoord = originPoint + i*stepsize*hvector + j*stepsize*gvector
                currentXYZ = pd.DataFrame({'Atom': self.xyztoASE(geomXYZ).get_chemical_symbols() , 'x': currentCoord[:,0], 'y': currentCoord[:,1], 'z': currentCoord[:,2]})
                self.writeXYZFile_orcaScan(currentXYZ, 'neb.xyz')
                
    def IterNum(self): # output the line num of the start of each iteration, for later use when splitting
        
        if self.outfile == 'None':
            
            return None
    
        iterList = []
        
        with open(self.outfile, 'r') as file:
            
            for lineNum,line in enumerate(list(file)):
                
                if '>>> EXPORT ITER = ' in line:
                    
                    iterList.append(lineNum)
        
        self.iterList = iterList
                    
        
        return iterList
    
    def sliceOutput(self,sliceIndex):
        self.iterList = self.IterNum()
        with open(self.outfile, 'r') as file:
            
            file = list(file)
            
            fileSlice = file[self.iterList[sliceIndex]:self.iterList[sliceIndex+1]] # slice out the block of corresponding iteration
        
        
        return fileSlice
    
    
    def getPyramidJsonData(self,startIndex,xyzFile='MOLCAS.md.xyz',Hamitonian='mcpdft',saveDir=r'D:\MolcasRead'):
        self.iterList = self.IterNum()
        finalIndex = startIndex
        with open(xyzFile, 'r') as file:
            
            sliceLength=self.atomNum+2
            file = list(file)
            num_of_slice = int(len(file)/sliceLength)
            # if Hamitonian == 'mcpdft': num_of_slice -= 1
            
            print(num_of_slice)
            
            for i in range(num_of_slice):
                
                startLine = i*sliceLength
                endLine = (i+1)*sliceLength
                
                xyzCartesian = file[startLine+2:endLine]
                                             
                xyzCartesian = pd.DataFrame(xyzCartesian)
                xyzCartesian.replace(r'\n','',regex=True,inplace = True) # remove all the \n character at the endo of each line
                
                # print(xyzCartesian)
                
                with open(saveDir+'/xyz'+str(i+startIndex)+'.txt', 'w') as f:
                    
                    xyzCartesian = xyzCartesian.to_string(header=False, index=False)

                    f.write(xyzCartesian)
                    
                    
        with open(self.outfile, 'r') as file:
            
            file = list(file)
                    
        print(len(self.iterList))
        for iteration in range(len(self.iterList)-1):
        # for iteration in range(10):
            # currentSlice = self.sliceOutput(iteration)
            currentSlice = file[self.iterList[iteration]:self.iterList[iteration+1]]
            state=0
            totalGradArr = []
            totalRASSFenergyArr = []
            for lineNum,line in enumerate(currentSlice):
                
                
                if 'Molecular gradients' in line:
                    state += 1 
                    xyzGrad = np.zeros((self.atomNum,3))
                    firstline = lineNum+8
                    lastline = lineNum+self.atomNum+7
                    
                    for i,j in enumerate(currentSlice[firstline:lastline+1]):
                        
                        print(j)
                        xyzGrad[i,0] = extractfloat(j,0)
                        xyzGrad[i,1] = extractfloat(j,1)
                        xyzGrad[i,2] = extractfloat(j,2)
                    
                    totalGradArr.append(xyzGrad)
                    
                if '::    RASSCF root number  ' in line and Hamitonian =='casscf':
                    
                    totalRASSFenergyArr.append(extractfloat(line, 1))
                    
                elif '::    CMS-PDFT Root' in line and Hamitonian == 'mcpdft':
                    
                    totalRASSFenergyArr.append(extractfloat(line, 1))
                    
            # print(totalRASSFenergyArr)
                    
            totalGradArr = np.array(totalGradArr) 
            totalRASSFenergyArr = np.array(totalRASSFenergyArr)
            totalGradArr = totalGradArr.reshape((state*self.atomNum*3,1))
            np.savetxt(saveDir+'/grad'+str(startIndex+iteration)+'.txt',totalGradArr)
            np.savetxt(saveDir+'/energy'+str(startIndex+iteration)+'.txt',totalRASSFenergyArr)
            finalIndex += 1

                    
        return finalIndex
                    
    def loadXYZ_Pyramid(self,filename):
        
        with open(filename,'r') as f:
            
           coord =  f.readlines()
           
           coord= pd.DataFrame(coord)
           coord.replace(r'\n','',regex=True,inplace = True) # remove all the \n character at the end of each line
           print(coord)
           coord = coord.values.tolist()
        
        return coord
        
        
    def compileJsonDataFile(self,molNum=1,folder=r'D:\MolcasRead'):
        
        totalGrad = []
        totalRASSFenergy = []
        totalXYZ = []
        totalNAC = []
        totalSOC = []
        totalDict = {'natom':self.atomNum, 'nstate': 2, 'nnac': 0, 'nsoc': 0}
        for i in range(molNum):
            
            energyFile = folder + r'/energy'+str(i)+'.txt'
            cartesianFile = folder +r'/xyz'+str(i)+'.txt'
            gradFile = folder + r'/grad'+str(i)+'.txt'
            
            energyArr = np.loadtxt(energyFile)
            gradArr = np.loadtxt(gradFile)
            gradArr = gradArr.reshape((2,self.atomNum,3))
            nacArr = []
            socArr = []
            coordArr = np.loadtxt(cartesianFile,dtype=np.str)
            # coordList = self.loadXYZ_Pyramid(cartesianFile)
            
            
            totalRASSFenergy.append(energyArr)
            totalGrad.append(gradArr)
            totalXYZ.append(coordArr)
            totalNAC.append(nacArr)
            totalSOC.append(socArr)
            
            
        totalDict['xyz'] = totalXYZ    
        totalDict['energy'] = totalRASSFenergy
        totalDict['grad'] = totalGrad
        totalDict['nac'] = totalNAC
        totalDict['soc'] = totalSOC
        
        totalDF = pd.DataFrame.from_dict(totalDict,orient='index')
        totalDF = pd.DataFrame.to_json(totalDF[0])
        
        with open(folder+'/egJson.json','w') as f:
            
            f.write(totalDF)
            
    def readMCPDFTOptimizationlog(self,nstate):
        
        with open(self.outfile,'r') as file:
            
            # we need to the energies of the final step, so reverse the whole file
            energyArr = np.zeros([nstate])
            reversedLogfile = list(file)
            reversedLogfile.reverse()
            for lineNum,line in enumerate(reversedLogfile):
                                
                if '::    CMS-PDFT Root    ' + str(nstate) in line: 
                    
                    lastline = lineNum # figure out where the last line of energy is

                    break
            for i in range(1,nstate+1):

                energyArr[i-1] = extractfloat(reversedLogfile[lastline+i-1], 1) # extract all states
        
            
            return np.flip(energyArr)
        
    def readCASSCFOptimizationlog(self,nstate):
        
        with open(self.outfile,'r') as file:
            
            # we need to the energies of the final step, so reverse the whole file
            energyArr = np.zeros([nstate])
            reversedLogfile = list(file)
            reversedLogfile.reverse()
            for lineNum,line in enumerate(reversedLogfile):
                                
                if '::    RASSCF root number  ' + str(nstate) in line: 
                    
                    lastline = lineNum # figure out where the last line of energy is

                    break
            for i in range(1,nstate+1):

                energyArr[i-1] = extractfloat(reversedLogfile[lastline+i-1], 1) # extract all states
        
            
            return np.flip(energyArr)
        
    def checkHappyLanding(self):
        
        with open(self.outfile,'r') as file:

            reversedLogfile = list(file)
            reversedLogfile.reverse()
            
            for lineNum,line in enumerate(reversedLogfile):
                
                if ' Happy landing! ' in line:

                    return True
            
            return False
        
class MoldenInterface(Interface):
    
    
    def __init__(self,moldenFilename):
        
        '''
        
        iniitializing the molden interface class gives you a list of readlines of molden data
        
        '''
        
        with open(moldenFilename) as file:
            
            self.molden = list(file)
            
    def replaceCoord(self,newCoord):
        
        '''
        a function that replace the FR_COORD section to a new coordinate, for ML dataset sampling
        input:
            
            newCoord: ase.Atoms format of xyz
        '''
        
        # first we figure out where the coordinates are in the molden file
        
        atomNum = len(newCoord)
        print(atomNum)
        
        for lineNum,line in enumerate(self.molden):
            
            if '[FR-COORD]' in line: 
                
                start = lineNum+1
                end = lineNum + atomNum +1
                
                break
            
        newMolden = self.molden.copy()
        atmIndex = 0
        
        for i in range(start,end):
            
            print(newMolden[i])
            newMolden[i] = str(newCoord.get_chemical_symbols()[atmIndex]) + ' ' + str(newCoord.get_positions()[atmIndex,0])+ ' '+str(newCoord.get_positions()[atmIndex,1])+ ' '+str(newCoord.get_positions()[atmIndex,2])+'\n'
            atmIndex += 1
            # print(newCoord[atmIndex])
        
        return newMolden
            
class SharcInterface(Interface):
    
    def getIndexPoint(self,fileList):

        indexNumList = []
        for lineNum,line in enumerate(fileList):
            if line[:5] == 'Index':
                indexNumList.append(lineNum)
    
        return indexNumList
        
    def parseInitcond(self,inicondFilename,atCond,reverse=False):
            
        with open(inicondFilename) as f:
            
            file = f.readlines()
            
        indexNumList = self.getIndexPoint(file)
        
        num = indexNumList[atCond]
        atomNum = extractfloat(file[2], 0)
        # print(atomNum)
        CoordArr = [file[i] for i in range(num+2,num+2+int(atomNum))]

        SumFile = []

        for i in CoordArr:
            line = re.split(r'[\s]',i)
            while '' in line: line.remove('')
            SumFile.append(line)
        df = pd.DataFrame(SumFile,columns=list('ABCDEFGHI'))
        coord = df.loc[:,list('ABCDEF')]
        coord.to_csv('geom',header=None, index=None, sep=' ', mode='w')
        if not reverse:
            vel = df.loc[:,list('GHI')]
            vel.to_csv('veloc',header=None, index=None, sep=' ', mode='w')
        else:
            vel = df.loc[:,list('GHI')]
            vel = vel.to_numpy().astype(np.float)
            vel = -vel
            np.savetxt("veloc", vel, delimiter=" ")
            
        return None
            
    def initializeSharcTraj(self,numCond,initcondPath,outputDir,templateDir):
                
        for i in range(numCond):
            
            currentDir=outputDir + 'TRAJ_0'+str(i)
            
            if os.path.isdir(currentDir): 
                shutil.rmtree(currentDir)
                
            shutil.copytree(templateDir, currentDir)
            os.chdir(currentDir)
            self.parseInitcond(initcondPath,i)
        
        os.chdir('..')
        
        return None
            
            
            
            
            
            
            
            
        
            
                    
                               
                

                     
                    
                
                
            
            
    
        
        
        

                
                
                
                
                
            
            
                
                
    

