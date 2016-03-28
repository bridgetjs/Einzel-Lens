from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functions as fc
loopmin=5
loopmax=6

oldstdout = sys.stdout
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'

N=20
InitialSize=fc.sigmax
Pathname=fc.Pathname
stdxerror=[]
Tfitarray=[]
ScreenPosArray=[]
LensPosArray=[]
eTfittedarray=[]

Temp=50
LensZ=50
ScreenPos=100

ParentFolderName='Pos=%d' %LensZ

AfileName='AVals/Adata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
BfileName='BVals/Bdata(%s,Screen=%d)..txt' %(ParentFolderName,ScreenPos)
        
fc.ABSet(AfileName,BfileName)
print fc.Afile,fc.Bfile
Afit,Bfit=fc.ABFIT()

stdx=[]
stdxerror=[]
V=[]
UI=[]

for i in range(loopmin,loopmax):
    sys.stdout=oldstdout
    VBackPlate=0-250*i;
    FolderName='VBack=%d' %VBackPlate
        
    #Make and change to the corresponding directory
    Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
    if os.path.exists(Path):
        os.chdir(Path);
        print os.getcwd()
    else:
        print "Shit's fucked with ",FolderName
        break

    if not os.path.exists("Plots"):
        os.makedirs("Plots")
        fc.Fisher()
        
    stdx2=[] #Array for single voltage standard deviations
    UI2=[] #Array for single voltage initial energies
    for j in range(0,2):
        beam=fc.BeamGenerator(N,Temp,InitialSize,fc.InitialPos)#N,T,sigma,central Z position
        fc.BeamDynWriter("beamdynTraj.in","beam","Snaps",beam,90)
        
        print "Running GPT in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        fc.GPTCall("beamdynTraj.in","time","std1.txt")
        
        fc.Plotter("std1.txt",FolderName,'xz','yz','Tz')


os.chdir("../../../../VaryParameters");



for i in range(loopmin,loopmax):
    sys.stdout=oldstdout
    VBackPlate=0-250*i;
    FolderName='VBack=%d' %VBackPlate
    
    #Make and change to the corresponding directory
    Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
    if os.path.exists(Path):
        os.chdir(Path);
        print os.getcwd()
    else:
        print "Shit's fucked with ",FolderName
        break
    
    if not os.path.exists("Plots"):
        os.makedirs("Plots")
        fc.Fisher()
    
    stdx2=[] #Array for single voltage standard deviations
    UI2=[] #Array for single voltage initial energies
    for j in range(1,2):
        beam=fc.BeamGenerator(N,Temp,InitialSize,fc.InitialPos)#N,T,sigma,central Z position
        fc.BeamDynWriter("beamdynTest.in","beam","Screens",beam,ScreenPos)

        print "Running GPT in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        fc.GPTCall("beamdynTest.in","position","std1.txt")

        finxarray,E=fc.Plotter("std1.txt",FolderName,'inx','finx','RealBunch')
        stdx2.append(np.std(finxarray))#Calculate std of final positions
        UI2.append(E) #Add intial energy to array
        print finxarray
    UI.append(np.mean(UI2)) #Calculate mean energy
    stdx.append(np.mean(stdx2)) #Calculate mean standard deviation
    stdxerror.append(np.std(stdx2))#Calculate error on the mean

    print 'U=',np.mean(UI2),'keV (',np.mean(stdx2)*1e3, '+/-', np.std(stdx2)*1e3, ') mm'
    
os.chdir("../../../../VaryParameters");
plt.show()
