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
    return x,y,z

TempRange=[10,20,30,40,50]

N,T,eT=data4('N.txt')

for i,No in enumerate(N):
    plt.figure(1)
    plt.errorbar(N[i],T[i],eT[i],fmt='.')
for Temp in TempRange:
    plt.figure(1)
    plt.axhline(y=Temp)

plt.figure(1)
plt.axis([0, 5500, 0, 60])
plt.xlabel('N particles')
plt.ylabel('Temperature (K)')
plt.show()
