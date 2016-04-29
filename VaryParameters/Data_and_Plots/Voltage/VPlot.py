from __future__ import division
import sys
sys.path.append('../')
import functions as fc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

dir='Data_and_Plots/Voltage/'
def data(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r')
    data=np.loadtxt(file)
    file.close()
    D=data[:,1]
    eD=data[:,2]
    return D,eD
Arangelist=[]
Brangelist=[]
RatioList=[]
ValList=[]
ValDiffs=[]
DiffList=[]

for indx,V2 in enumerate(fc.VoltageRange):
    
    ParentFolderName='Voltage=%d' %V2
    
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName
    
    A,eA=data(AfileName)
    B,eB=data(BfileName)
    
    fc.ABSet(AfileName,BfileName)
    Afit,Bfit=fc.ABFIT()
#    fc.ABCheck()
    Arange=np.max(A)-np.min(A)
#    eA=np.sqrt(  eA[A.index(max(A))]**2   + eA[A.index(min(A))]**2    )
    Brange=np.max(B)-np.min(B)
    if len(Afit.roots())==1 and len(Bfit.roots())==1:
        DiffList.append((Afit.roots() - Bfit.roots())[0])
        ValDiffs.append(V2*1e-3)

    Arangelist.append(Arange)
    Brangelist.append(Brange)
    RatioList.append((Brange/Arange)**2)
    ValList.append(V2*1e-3)

plt.figure(1)
plt.plot(ValList,Brangelist,'.')
plt.figure(3)
plt.plot(ValList,Arangelist,'.')


plt.figure(1)
plt.xlabel('Central Voltage (kV)')
plt.ylabel('Range of B')
plt.xlim([2.5,5.5])
plt.savefig(dir+'BRange.eps')


plt.figure(3)
plt.xlabel('Central Voltage (kV)')
plt.ylabel('Range of A')
plt.xlim([2.5,5.5])
plt.savefig(dir+'ARange.eps')


plt.figure(4)
plt.plot(ValDiffs,DiffList,'.')
plt.xlim([2.5,5.5])
plt.xlabel('Central Voltage (kV)')
plt.ylabel(r'Difference between $A_0$ $B_0$ (keV)')
plt.savefig(dir+'A0 B0 difference.eps')

#plt.show()