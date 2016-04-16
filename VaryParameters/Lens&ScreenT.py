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
A=1
Aplot=1
Eplot=0
oldstdout = sys.stdout
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
TempRange=[25]

#TempRange=np.linspace(25,150,6)
Number=[1000,2000,5000,7500,10000]
N=500
InitialSize=fc.sigmax
stdxerror=[]
Tfitarray=[]
ScreenPosArray=[]
LensPosArray=[]
eTfittedarray=[]
Narray=[]
Pathname="SF_Files/LensPosVar"
LensZRange=[50]
#fc.LensZRange
ScreenPosRange=[100]
#fc.ScreenPosRange
ScreenPos=100
LensZ=50

for Temp in TempRange:
    
    for ScreenPos in fc.ScreenPosRange:
   
        for LensZ in fc.LensZRange:
            
            stdx=[]
            stdxerror=[]
            V=[]
            UI=[]
            ParentFolderName='Pos=%d' %LensZ

            if ScreenPos<=LensZ: continue
            
            print 'Running T Simulation with Lens at %d and Screen at %d' %(LensZ, ScreenPos)
            
            AfileName='AVals/Adata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
            BfileName='BVals/Bdata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
            
            fc.ABSet(AfileName,BfileName)
            print fc.Afile,fc.Bfile
            Afit,Bfit=fc.ABFIT()
            for i in range(loopmin,loopmax):
                sys.stdout=oldstdout
                VBackPlate=0-250*i;
                V.append(-VBackPlate*1e-3)
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
                for j in range(0,10):
                    beam=fc.BeamGenerator(N,Temp,InitialSize,fc.InitialPos)#N,T,sigma,central Z position
                    fc.BeamDynWriter("beamdynMRFCTest.in","beam","Screens",beam,ScreenPos)
                
                    print "Running GPT in "+ FolderName
                #                DynamicFile,     GroupBy,    Outtxt
                    fc.GPTCall("beamdynMRFCTest.in","position","std1.txt",'std')
                    
                    ScreenSTD,E=fc.Plotter("std1.txt",FolderName,'RealBunchStds')
                    stdx2.append(ScreenSTD)#Calculate std of final positions
                    UI2.append(E) #Add intial energy to array

                UI.append(np.mean(UI2)) #Calculate mean energy
                stdx.append(np.mean(stdx2)) #Calculate mean standard deviation
                stdxerror.append(np.std(stdx2))#Calculate error on the mean

                print 'U=',np.mean(UI2),'keV (',np.mean(stdx2)*1e3, '+/-', np.std(stdx2)*1e3, ') mm'
                
                os.chdir("../../../../VaryParameters")
                    
#            plt.figure(30)
#            plt.errorbar(np.asarray(UI),np.asarray(stdx),np.asarray(stdxerror),fmt='.',markersize=0)
            Tfit=optimization.curve_fit(fc.Model2, UI, stdx,sigma=stdxerror, maxfev=10000,p0=Temp)
            Tfitted=Tfit[0][0]
            eTfitted=np.sqrt(np.diag(Tfit[1])[0])
            print 'Temperature = %2.3f +/- %2.3f' %(Tfitted,eTfitted)
#            plt.plot(UI,fc.PlotModel(UI,Tfitted,Afit,Bfit,InitialSize),label='T=(%2.1f $\pm$ %2.1f)K' %(Tfitted,eTfitted))
            Tfitarray.append(Tfitted)
            eTfittedarray.append(eTfitted)
            ScreenPosArray.append(ScreenPos)
            LensPosArray.append(LensZ)
            
            plt.figure(31)
            plt.errorbar(LensZ,Tfitted,eTfitted,fmt='.')


    plt.figure(31)
    plt.axhline(y=Temp)

print Narray,Tfitarray,eTfittedarray
S=LensPosArray,ScreenPosArray,Tfitarray,eTfittedarray
datafile=open('Screen&Lens.txt','w')
np.savetxt(datafile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
datafile.close()

#plt.figure(30)
#plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
#plt.ylabel('Rms Size of Bunch (m)')
#plt.legend()
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.eps'%InitialSize)
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.png'%InitialSize)


plt.figure(31)
plt.xlabel('Number of particles')
plt.ylabel('Fitted Temperature (K)')
plt.axis([10, 100, -25, 50])
plt.savefig('FittedTempvsN.eps')
plt.show()
#Test of git
#Test 2
