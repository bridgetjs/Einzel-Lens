from __future__ import division
import sys
sys.path.append('../')
import functions as fc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
LensZRange=fc.LensZRange

def data(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r')
    data=np.loadtxt(file)
    file.close()
    D=data[:,1]
    return D

for LensZ in LensZRange:
    ParentFolderName='Pos=%d' %LensZ
    AfileName='Avals/Adata(%s).txt' %ParentFolderName
    A=data(AfileName)
    
    Arange=np.max(A)-np.min(A)

    BFolderName='Pos=%d' %LensZ
    BfileName='Bvals/Bdata(%s).txt' %BFolderName
    B=data(BfileName)
    Brange=np.max(B)-np.min(B)
#    print LensZ,Brange,Arange

    plt.figure(1)
    plt.plot(1e-2*LensZ,Brange,'k.')

    plt.figure(2)
    plt.plot(1e-2*LensZ,(Brange/Arange)**2,'k.')

    plt.figure(3)
    plt.plot(1e-2*LensZ,Arange,'k.')
plt.figure(1)
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of B')

plt.figure(2)
plt.xlabel('Lens Position(m)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')

plt.figure(3)
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of A')


plt.show()