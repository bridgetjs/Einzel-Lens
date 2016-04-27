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
import timeit
int=0

loopmin=1
loopmax=21
A=1
Aplot=1
Eplot=0
oldstdout = sys.stdout
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
#TempRange=[10,10,10,20,20,20,30,30,30,40,40,40,50,50,50]
TempRange=[10,20,30]
#TempRange=np.linspace(25,150,6)
Number=[2000,3000,4000,5000]
N=2000
InitialSize=fc.sigmax
stdxerror=[]
Tfitarray=[]
ScreenPosArray=[]
LensPosArray=[]
eTfittedarray=[]
Narray=[]
SeptArray=[]
Pathname=fc.Pathname
LensZRange=[50]
#fc.LensZRange
ScreenPosRange=[100]
#fc.ScreenPosRange
ScreenPos=100
LensZ=50
Flighttubelength=[]
for Temp in TempRange:
   
#    for V2 in fc.VoltageRange:
#    for N in Number:

    for k in range(0,2):
        start_time = timeit.default_timer()

        if k==0:
            ParentFolderName='Tuned2'
            ScreenPos=100
      
        if k==1:
            ParentFolderName='HalfLength'
            ScreenPos=50


        
#        ParentFolderName='Tuned2'
#        Flighttubelength.append(ScreenPos)
#        ParentFolderName='Sept=%1.2f' %PlateSept
#        ParentFolderName='Sept=%1.2f' %PlateSept
#        ParentFolderName='HalfLength'
        stdx=[]
        stdxerror=[]
        V=[]
        UI=[]
#        ParentFolderName='Pos=%d' %LensZ



#        if ScreenPos<=LensZ: continue

        print 'Running T Simulation with Lens at %d and Screen at %d' %(LensZ, ScreenPos)
        
        AfileName='AVals/Adata(%s).txt' %(ParentFolderName)
        BfileName='BVals/Bdata(%s).txt' %(ParentFolderName)
        
        fc.ABSet(AfileName,BfileName)
        print fc.Afile,fc.Bfile
        
        Afit,Bfit=fc.ABFIT()
        
        if int==0:
            fc.ABCheck()
            int+=1
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
                sys.exit()
            
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
                
        plt.figure(36)
        plt.errorbar(np.asarray(UI),np.asarray(stdx),np.asarray(stdxerror),fmt='.',markersize=0)
        Tfit=optimization.curve_fit(fc.Model2, UI, stdx,sigma=stdxerror, maxfev=100000,p0=Temp)
        Tfitted=Tfit[0][0]
        eTfitted=np.sqrt(np.diag(Tfit[1])[0])
        print 'Temperature = %2.3f +/- %2.3f' %(Tfitted,eTfitted)
        plt.plot(UI,fc.PlotModel(UI,Tfitted,Afit,Bfit,InitialSize),label='T=(%2.1f $\pm$ %2.1f)K' %(Tfitted,eTfitted))
        Tfitarray.append(Tfitted)
        eTfittedarray.append(eTfitted)
#        ScreenPosArray.append(ScreenPos)
#        LensPosArray.append(LensZ)
#        Narray.append(N)
#        SeptArray.append(PlateSept)
        plt.figure(35)
        plt.errorbar(ScreenPos,Tfitted,eTfitted,fmt='.')
        elapsed = timeit.default_timer() - start_time
        print elapsed
            
    plt.figure(35)
    plt.axhline(y=Temp)

#print SeptArray,Tfitarray,eTfittedarray
#S=Flighttubelength,Tfitarray,eTfittedarray
#datafile=open('FlightTube.txt','w')
#np.savetxt(datafile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
#datafile.close()

plt.figure(36)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.eps'%InitialSize)
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.png'%InitialSize)

#
plt.figure(35)
plt.xlabel('Flight Tube Length (cm)')
plt.ylabel('Fitted Temperature (K)')
plt.axis([0, 100, 0, 100])
#plt.savefig('FittedTemp.eps')
plt.show()
##Test of git
##Test 2
