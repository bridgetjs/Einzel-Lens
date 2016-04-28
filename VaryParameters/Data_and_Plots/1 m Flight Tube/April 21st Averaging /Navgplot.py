from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math

def data4(filename): #function to return the data file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2];   #read data assuming data is in the exepcted format
    return x,y,z,data[:,3],data[:,4]

TempRange=[10,30,50]

N,T,eT,sT,q=data4('NandTAvg.txt')

for i,No in enumerate(N):
    plt.figure(1)
    plt.errorbar(N[i],T[i],eT[i],fmt='.')

    plt.figure(2)
    plt.errorbar(N[i],T[i],sT[i],fmt='.')

    plt.figure(3)
    plt.errorbar(N[i],T[i],q[i],fmt='.')

for Temp in TempRange:
    plt.figure(1)
    plt.axhline(y=Temp)

    plt.figure(2)
    plt.axhline(y=Temp)

    plt.figure(3)
    plt.axhline(y=Temp)

plt.figure(1)
plt.axis([0, 5500, 0, 60])
plt.xlabel('N particles')
plt.ylabel('Fitted Temperature (K)')
plt.savefig('AveragedT.eps')
plt.savefig('AveragedT.png')

plt.figure(2)
plt.axis([0, 5500, 0, 60])
plt.xlabel('N particles')
plt.ylabel('Fitted Temperature (K) - standard deviation ')
plt.savefig('AveragedTstd.eps')
plt.savefig('AveragedTstd.png')

plt.figure(3)
plt.axis([0, 5500, 0, 60])
plt.xlabel('N particles')
plt.ylabel('Temperature (K)')
plt.ylabel('Fitted Temperature (K) - quadrature')
plt.savefig('AveragedTquadrature.eps')
plt.savefig('AveragedTquadrature.png')
