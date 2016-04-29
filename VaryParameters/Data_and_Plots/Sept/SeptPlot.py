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

def data4(filename): #function to return the data file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2];   #read data assuming data is in the exepcted format
    return x,y,z
dir='Data_and_Plots/Sept/'
Arangelist=[]
Brangelist=[]
RatioList=[]
Sept_List=[]
SeptDiffs=[]
DiffList=[]

TempRange=[10,30]

S,T,eT=np.loadtxt('Data_and_Plots/Sept/Data.txt',unpack=True)

plt.figure(1)
plt.errorbar(S,T,eT,fmt='.')
for Temp in TempRange: plt.axhline(y=Temp)
plt.axis([0, 2, 0, 40])
plt.xlabel('Plate Separation,s, (cm)')
plt.ylabel('Temperature (K)')
plt.savefig(dir+'fits.eps')
for PlateSept in fc.TotSeptRange:
    ParentFolderName='Sept=%1.2f' %PlateSept
    
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName
        
    UA,A,eA=np.loadtxt(AfileName,unpack=True)
    UB,B,eB=np.loadtxt(BfileName,unpack=True)

    fc.ABSet(AfileName,BfileName)
    Afit,Bfit=fc.ABFIT()
    #    fc.ABCheck()
    Arange=np.max(A)-np.min(A)
    #    eA=np.sqrt(  eA[A.index(max(A))]**2   + eA[A.index(min(A))]**2    )
    Brange=np.max(B)-np.min(B)
    if len(Afit.roots())==1 and len(Bfit.roots())==1:
        DiffList.append((Afit.roots() - Bfit.roots())[0])
        SeptDiffs.append(PlateSept)
    
    Arangelist.append(Arange)
    Brangelist.append(Brange)
    RatioList.append((Brange/Arange)**2)
    Sept_List.append(PlateSept)

plt.figure(2)
plt.plot(Sept_List,Brangelist,'.')
plt.xlabel('Plate Separation (cm)')
plt.ylabel('Range of B')
plt.xlim([0,2.5])
plt.savefig(dir+'BRange.eps')


plt.figure(3)
plt.plot(Sept_List,Arangelist,'.')
plt.xlabel('Plate Separation (cm)')
plt.ylabel('Range of A')
plt.xlim([0,2.5])
plt.savefig(dir+'ARange.eps')


plt.figure(4)
plt.plot(SeptDiffs,DiffList,'.')
plt.xlim([0,2.5])
plt.xlabel('Plate Separation (cm)')
plt.ylabel(r'Difference between $A_0$ $B_0$ (keV)')
plt.savefig(dir+'A0 B0 difference.eps')

