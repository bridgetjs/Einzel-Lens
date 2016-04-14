from __future__ import division
import sys
sys.path.append('../../')
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
    eD=data[:,2]
    return D,eD


def Alldata(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r')
    data=np.loadtxt(file)
    file.close()
    U=data[:,0]
    D=data[:,1]
    eD=data[:,2]
    return U,D,eD


Arangelist=[]
Brangelist=[]
RatioList=[]
Sept_List=[]


AfileName='../Avals/Adata(Optimised).txt'
BfileName='../BVals/Bdata(Optimised).txt'


UA,A,eA=Alldata(AfileName)
UB,B,eB=Alldata(BfileName)

fc.ABSet(AfileName,BfileName)

Afit,Bfit=fc.ABFIT()

plt.figure(1)
plt.plot(UA,Afit(UA))
plt.plot(UA,A,'.')
plt.xlabel('U (keV)')
plt.ylabel('A (no units)')


plt.figure(2)
plt.plot(UB,Bfit(UB))
plt.plot(UB,B,'.')
plt.xlabel('U (keV)')
plt.ylabel('B m/rad')

Arange=np.max(A)-np.min(A)
#    eA=np.sqrt(  eA[A.index(max(A))]**2   + eA[A.index(min(A))]**2    )
Brange=np.max(B)-np.min(B)

print Arange,Brange

#Arangelist.append(Arange)
#Brangelist.append(Brange)
#RatioList.append((Brange/Arange)**2)
#Sept_List.append(PlateSept)

#plt.figure(1)
#plt.plot(Sept_List,Brangelist,fc.colourstring[1])
#plt.figure(2)
#plt.plot(Sept_List,RatioList,fc.colourstring[2])
#plt.figure(3)
#plt.plot(Sept_List,Arangelist,fc.colourstring[3])
#
#
#plt.figure(1)
#plt.legend()
#plt.xlabel('Plate Seperation (cm)')
#plt.ylabel('Range of B')
#plt.savefig('BSept.eps')
#
#plt.figure(2)
#plt.xlabel('Plate Seperation (cm)')
#plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')
#plt.savefig('BA_Sept.eps')
#
#plt.figure(3)
#plt.xlabel('Plate Seperation (cm)')
#plt.ylabel('Range of A')
#plt.legend()
#plt.savefig('A_Sept.eps')
plt.show()