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

TempRange=[10,30]

S,T,eT=np.loadtxt('Data_and_Plots/Seperation/Data.txt',unpack=True)

plt.figure(1)
plt.errorbar(S,T,eT,fmt='.')
for Temp in TempRange: plt.axhline(y=Temp)

plt.figure(1)
plt.axis([0, 2, 0, 40])
plt.xlabel('Plate Separation,s, (cm)')
plt.ylabel('Temperature (K)')
plt.show()
