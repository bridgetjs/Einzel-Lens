from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import functions as fc
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
loopmin=1
loopmax=21

oldstdout = sys.stdout
ParentFolderName="Acc2"
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
TempRange=[25]
#TempRange=np.linspace(25,150,6)
Temp=25
N=1000#Plenty
AfileName='../GPT_A/Adata.txt'
BfileName='../GPT_B/Bdata.txt'
fc.ABSet(AfileName,BfileName)
Afit,Bfit=fc.ABFIT()
Number=[50,100,250,500,750,1000,1250,1500,1750,2000,2500,3500,5000,6000,7500,8500,10000]

InitialSize=fc.sigmax

for indx,N in enumerate(Number):
    print N
    V=[]
    UI=[]
    UF=[]
    
    stdxerror=[]
    Tfitarray=[]
    eTfittedarray=[]
    for j in range(0,1):
        stdx=[]
        V=[]
        UI=[]
        for i in range(loopmin,loopmax):
            sys.stdout=oldstdout
            VBackPlate=0-250*i;
            V.append(-VBackPlate*1e-3)
            FolderName='VBack=%d' %VBackPlate

            #Make and change to the corresponding directory
            Path="../" + ParentFolderName+"/"+FolderName
            if os.path.exists(Path):
                os.chdir(Path);
            else:
                print "Shit's fucked with ",FolderName
                break
            
            if not os.path.exists("Plots"):
                os.makedirs("Plots")
            fc.Fisher()
            
            stdx2=[] #Array for single voltage standard deviations
            UI2=[] #Array for single voltage initial energies
            for j in range(0,10):
                beam=fc.BeamGenerator(N,Temp,InitialSize,0.025)#N,T,sigma,central Z position
                fc.BeamDynWriter("beamdynMRFCTest.in","beam","Screens",beam)
            
                print "Running GPT in "+ FolderName
            #                DynamicFile,     GroupBy,    Outtxt
                fc.GPTCall("beamdynMRFCTest.in","position","std1.txt")
                
                finxarray,E=fc.Plotter("std1.txt",FolderName,i,'inx','finx','RealBunch')
                stdx2.append(np.std(finxarray))#Calculate std of final positions
#            UI.append(E) #Add intial energy to array

            UI.append(np.mean(E)) #Calculate mean energy
            stdx.append(np.mean(stdx2)) #Calculate mean standard deviation
            stdxerror.append(np.std(stdx2))#Calculate error on the mean

            print 'U=',np.mean(E),'keV (',np.mean(stdx2)*1e3, '+/-', np.std(stdx2)*1e3, ') mm'
        

    #        stdx.append(fc.rms(finxarray))
    #


            os.chdir("../../TSim");

        #print sigma

        sigmafile=open("stdxdata.txt",'w')
        S= UI , stdx
    #    np.savetxt(sigmafile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
        plt.figure(30)
#        plt.plot(np.asarray(UI),np.asarray(stdx),'.',markersize=0)
        plt.errorbar(np.asarray(UI),np.asarray(stdx),np.asarray(stdxerror),fmt='.',markersize=0)
        Tfit=optimization.curve_fit(fc.Model2, UI, stdx,sigma=stdxerror, maxfev=10000,p0=Temp)
#        print Tfit
        Tfitted=Tfit[0][0]
        eTfitted=np.sqrt(np.diag(Tfit[1])[0])
#        print Tfitted,eTfitted
        print 'Temperature = %2.3f +/- %2.3f' %(Tfitted,eTfitted)
        plt.plot(UI,fc.PlotModel(UI,Tfitted,Afit,Bfit,InitialSize),colourstring[indx],label='T=(%2.1f $\pm$ %2.1f)K' %(Tfitted,eTfitted))
        Tfitarray.append(Tfitted)
        eTfittedarray.append(eTfitted)
        plt.figure(31)
        plt.errorbar(N,Tfitted,eTfitted)

#    plt.figure(31)
#    print Tfitarray,[1/e**2 for e in eTfittedarray],eTfittedarray
#    TWeight,eTWeight=np.average(Tfitarray,weights=[1/e**2 for e in eTfittedarray],returned=True)
#    plt.errorbar(N,TWeight,eTWeight)

plt.figure(31)
plt.axhline(y=Temp)
plt.xlabel('Number of Particles')
plt.ylabel('Fitted Temperature (K)')
plt.savefig('FittedTempvsN.eps')
plt.axis([0, 12000, 0, 50])

plt.figure(30)
plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
plt.ylabel('Rms Size of Bunch (m)')
plt.legend()
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.eps'%InitialSize)
#plt.savefig('rmsSize&T(initalsize=%1.2e)910.png'%InitialSize)
plt.show()