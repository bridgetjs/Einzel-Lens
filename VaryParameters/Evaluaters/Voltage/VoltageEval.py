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
Arangelist=[]
Brangelist=[]
RatioList=[]
VList=[]

for indx,V2 in enumerate(fc.VoltageRange):

    
    ParentFolderName='Voltage=%d' %V2
    
    AfileName='../AVals/Adata(Voltage=%d).txt' %(V2)
    BfileName='../BVals/Bdata(Voltage=%d).txt' %(V2)
    
    
    A,eA=data(AfileName)
    B,eB=data(BfileName)
    
    
    Arange=np.max(A)-np.min(A)
#    eA=np.sqrt(  eA[A.index(max(A))]**2   + eA[A.index(min(A))]**2    )
    Brange=np.max(B)-np.min(B)

    Arangelist.append(Arange)
    Brangelist.append(Brange)
    RatioList.append((Brange/Arange)**2)
    VList.append(V2)

print eA,eB
plt.figure(1)
plt.plot(VList,Brangelist,'.')
plt.axis([3900,5100,5,5.5])


plt.figure(2)
plt.plot(VList,RatioList,'.')
plt.axis([3900,5100,0,0.5])

plt.figure(3)
plt.plot(VList,Arangelist,'.')
plt.axis([3900,5100,9,10])

plt.figure(1)
plt.legend()
plt.xlabel('Central Plate Voltage (V)')
plt.ylabel('Range of B')
plt.savefig('BVoltage_fine.eps')

plt.figure(2)
plt.xlabel('Central Plate Voltage (V)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')
plt.savefig('BAVoltage_fine.eps')

plt.figure(3)
plt.xlabel('Central Plate Voltage (V)')
plt.ylabel('Range of A')
plt.legend()
plt.savefig('AVoltage_fine.eps')


plt.show()