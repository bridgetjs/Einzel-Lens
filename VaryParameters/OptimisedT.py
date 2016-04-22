from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
from itertools import izip
sys.path.append('../')
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functions as fc
loopmin=1
loopmax=21

oldstdout = sys.stdout

TempRange=[10,10,10,30,30,30,50,50,50]

#TempRange=np.linspace(25,150,6)
Number=[2000,3000,4000,5000]

InitialSize=fc.sigmax
stdxerror=[]
Tfitarray=[]
ScreenPosArray=[]
LensPosArray=[]
eTfittedarray=[]
sTfittedarray=[]
Narray=[]
Pathname=fc.Pathname
LensZRange=[50]
#fc.LensZRange
ScreenPosRange=[100]
#fc.ScreenPosRange
ScreenPos=100
LensZ=50

Tfit10=[]
Tfit30=[]
Tfit50=[]

eTfit10=[]
eTfit30=[]
eTfit50=[]

sTfit10=[]
sTfit30=[]
sTfit50=[]

qTfit10=[]
qTfit30=[]
qTfit50=[]


for N in Number:
    T10=[]
    wT10=[]

    T30=[]
    wT30=[]

    T50=[]
    wT50=[]
    
    for Temp in TempRange:
        
        stdx=[]
        stdxerror=[]
        V=[]
        UI=[]
        ParentFolderName='Tuned2'
        
#        if ScreenPos<=LensZ: continue

        print 'Running T Simulation with Lens at %d and Screen at %d' %(LensZ, ScreenPos)
        print 'N=%d,T=%d' %(N,Temp)
        AfileName='AVals/Adata(Tuned2).txt'
        BfileName='BVals/Bdata(Tuned2).txt'
        
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
                print "Shit's fucked with ",Path
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
                
        plt.figure(30)
        plt.errorbar(np.asarray(UI),np.asarray(stdx),np.asarray(stdxerror),fmt='.',markersize=0)
        Tfit=optimization.curve_fit(fc.Model2, UI, stdx,sigma=stdxerror, maxfev=10000,p0=Temp)
        Tfitted=Tfit[0][0]
        eTfitted=np.sqrt(np.diag(Tfit[1])[0])
        print 'Temperature = %2.3f +/- %2.3f' %(Tfitted,eTfitted)
        plt.plot(UI,fc.PlotModel(UI,Tfitted,Afit,Bfit,InitialSize),label='T=(%2.1f $\pm$ %2.1f)K' %(Tfitted,eTfitted))
#        Tfitarray.append(Tfitted)
#        eTfittedarray.append(eTfitted)
        ScreenPosArray.append(ScreenPos)
        LensPosArray.append(LensZ)
        Narray.append(N)

        if Temp==10:
            T10.append(Tfitted)
            wT10.append(1/eTfitted**2)
        if Temp==30:
            T30.append(Tfitted)
            wT30.append(1/eTfitted**2)
        if Temp==50:
            T50.append(Tfitted)
            wT50.append(1/eTfitted**2)

#        plt.figure(31)
#        plt.errorbar(N,Tfitted,eTfitted,fmt='.')
#        plt.axhline(y=Temp)


    Tfit10.append(np.average(T10,weights=wT10))
    Tfit30.append(np.average(T30,weights=wT30))
    Tfit50.append(np.average(T50,weights=wT50))

    sTfit10.append(np.sqrt(3/2)*np.sqrt(1/3 * np.sum([ (T-10)**2 for T in T10])) )
    sTfit30.append(np.sqrt(3/2)*np.sqrt(1/3 * np.sum([ (T-30)**2 for T in T30])))
    sTfit50.append(np.sqrt(3/2)*np.sqrt(1/3 * np.sum([ (T-50)**2 for T in T50])))

    eTfit10.append(1/np.sqrt(np.sum(wT10)))
    eTfit30.append(1/np.sqrt(np.sum(wT30)))
    eTfit50.append(1/np.sqrt(np.sum(wT50)))


qTfit10=[np.sqrt(a**2+b**2) for a,b in izip(eTfit10,sTfit10) ]
qTfit30=[np.sqrt(a**2+b**2) for a,b in izip(eTfit30,sTfit30) ]
qTfit50=[np.sqrt(a**2+b**2) for a,b in izip(eTfit50,sTfit50) ]



print  Number*3, Tfit10+Tfit30+Tfit50 ,eTfit10+eTfit30+eTfit50 ,sTfit10+sTfit30+sTfit50 ,qTfit10+qTfit30+qTfit50
S= Number*3, Tfit10+Tfit30+Tfit50 ,eTfit10+eTfit30+eTfit50 ,sTfit10+sTfit30+sTfit50 ,qTfit10+qTfit30+qTfit50
datafile=open('NandTAvg.txt','w')
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
plt.errorbar(Number,np.asarray(Tfit10),np.asarray(eTfit10),fmt='r.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit30),np.asarray(eTfit30),fmt='g.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit50),np.asarray(eTfit50),fmt='b.',markersize=0)
plt.axhline(y=10)
plt.axhline(y=30)
plt.axhline(y=50)
plt.axis([0, max(Number)+10, 0, 70])
plt.savefig('FittedTempvsN.eps')

plt.figure(32)
plt.xlabel('Number of particles')
plt.ylabel('Fitted Temperature (K) - standard deviation error')
plt.errorbar(Number,np.asarray(Tfit10),np.asarray(sTfit10),fmt='r.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit30),np.asarray(sTfit30),fmt='g.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit50),np.asarray(sTfit50),fmt='b.',markersize=0)
plt.axhline(y=10)
plt.axhline(y=30)
plt.axhline(y=50)
plt.axis([0, max(Number)+10, 0, 70])
plt.savefig('FittedTemp(standard deviation)vsN.eps')

plt.figure(33)
plt.xlabel('Number of particles')
plt.ylabel('Fitted Temperature (K) - quadrature errors')
plt.errorbar(Number,np.asarray(Tfit10),np.asarray(qTfit10),fmt='r.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit30),np.asarray(qTfit30),fmt='g.',markersize=0)
plt.errorbar(Number,np.asarray(Tfit50),np.asarray(qTfit50),fmt='b.',markersize=0)
plt.axhline(y=10)
plt.axhline(y=30)
plt.axhline(y=50)
plt.axis([0, max(Number)+10, 0, 70])
plt.savefig('FittedTemp(quadrature)vsN.eps')


plt.show()
#Test of git
#Test 2
