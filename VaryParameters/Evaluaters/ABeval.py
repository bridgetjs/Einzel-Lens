from __future__ import division
import sys
sys.path.append('../')
import functions as fc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

LensZRange=fc.LensZRange
ScreenPosRange=fc.ScreenPosRange

def data(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r')
    data=np.loadtxt(file)
    file.close()
    D=data[:,1]
    return D

def Alldata(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r')
    data=np.loadtxt(file)
    file.close()
    U=data[:,0]
    D=data[:,1]
    return U,D

for indx,ScreenPos in enumerate(ScreenPosRange):
    Arangelist=[]
    Brangelist=[]
    RatioList=[]
    Zlist=[]
    for indx2,LensZ in enumerate(LensZRange):
        
        if ScreenPos<=LensZ: continue
        
        ParentFolderName='Pos=%d' %LensZ
        
        AfileName='AVals/Adata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
        BfileName='BVals/Bdata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
        
        
        UA,A=Alldata(AfileName)
        UB,B=Alldata(BfileName)
        
        fc.ABSet(AfileName,BfileName)
        Afit,Bfit=fc.ABFIT()
        
        Arange=np.max(A)-np.min(A)
        Brange=np.max(B)-np.min(B)

        Arangelist.append(Arange)
        Brangelist.append(Brange)
        RatioList.append((Brange/Arange)**2)
        Zlist.append(LensZ)
    
    
        plt.figure(6)
        plt.plot(UA,Afit(UA))
        plt.plot(UA,A,'.')
        plt.xlabel('U (keV)')
        plt.ylabel('A (no units)')
    
        plt.figure(5)
        plt.plot(UB,Bfit(UB))
        plt.plot(UB,B,'.')
        plt.xlabel('U (keV)')
        plt.ylabel('B m/rad')

    plt.figure(1)
    plt.plot(Zlist,Brangelist,fc.colourstring[indx],label='Pos=%d' %ScreenPos)
    plt.figure(2)
    plt.plot(Zlist,RatioList,fc.colourstring[indx],label='Pos=%d' %ScreenPos)
    plt.figure(3)
    plt.plot(Zlist,Arangelist,fc.colourstring[indx],label='Pos=%d' %ScreenPos)


plt.figure(1)
plt.legend()
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of B')

plt.figure(2)
plt.xlabel('Lens Position(m)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')
plt.legend()

plt.figure(3)
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of A')
plt.legend()


plt.show()