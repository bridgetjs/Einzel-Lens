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
loopmin=1
loopmax=21

oldstdout = sys.stdout
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'

N=2000
InitialSize=fc.sigmax
Pathname=fc.Pathname
stdxerror=[]
Tfitarray=[]
ScreenPosArray=[]
LensPosArray=[]
eTfittedarray=[]

Temp=10
LensZ=50
ScreenPos=100

ParentFolderName='Tuned2'

stdx=[]
stdxerror=[]
V=[]
UI=[]
UIe=[]

for i in range(loopmin,loopmax):
    sys.stdout=oldstdout
    VBackPlate=0-250*i;
    V.append(abs(VBackPlate)/1e3)
    FolderName='VBack=%d' %VBackPlate
        
    #Make and change to the corresponding directory
    Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
    if os.path.exists(Path):
        os.chdir(Path);
        print os.getcwd()
    else:
        print "Shit's fucked with ",Path
        break

    if not os.path.exists("Plots"):
        os.makedirs("Plots")
        fc.Fisher()
        
    stdx2=[] #Array for single voltage standard deviations
    UI2=[] #Array for single voltage initial energies
    for j in range(0,1):
        
        beam=fc.BeamGenerator(N,Temp,InitialSize,fc.InitialPos)#N,T,sigma,central Z position
        fc.BeamDynWriter("beamdynTraj.in","beam","Snaps",beam)
        
        print "Running GPT in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        fc.GPTCall("beamdynTraj.in","time","std1.txt",'std')
#        fc.Plotter("std1.txt",FolderName,'xz','yz','Tz','avgs','show')

        fc.GPTCall("beamdynMRFCTest.in","position","std1.txt",'std')
        ScreenSTD,E=fc.Plotter("std1.txt",FolderName,'RealBunchStds')
        stdx2.append(ScreenSTD)#Calculate std of final positions
        UI2.append(E) #Add intial energy to array
#        print E
#
#    print UI2
    UI.append(np.mean(UI2)) #Calculate mean energy
    UIe.append(np.std(UI2))
    os.chdir("../../../../VaryParameters");

plt.plot(V,UI,'.')
line = optimization.curve_fit(fc.line, V, UI)
print line[0][0], np.sqrt(np.diag(line[1])[0])
plt.plot(np.asarray(V),np.asarray(V)*line[0][0]+line[0][1],'k',label=r'U=(%1.4f $\pm$ %1.4f )$|V_{acc}|$' %(line[0][0],np.sqrt(np.diag(line[1])[0])))
plt.xlabel(r'$|V_{acc}|$ (kV)',fontsize=14)
plt.ylabel(r'$U_{bunch}$ (keV)',fontsize=14)
plt.legend(loc=2,fontsize=14)
print UI,UIe
plt.show()




