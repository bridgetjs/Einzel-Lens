from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functions as fc

Screen=50
colourstring='bgrcmykbgrcmybgrcmybgrcmybgrcmybgrcmy'

def data4(filename): #function to return the data file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2];   #read data assuming data is in the exepcted format
    return x,y,z


N,T,eT=data4('N.txt')
for i,No in enumerate(N):
    plt.figure(1)
    plt.errorbar(N[i],T[i],eT[i],fmt='.')

plt.figure(1)
plt.axhline(y=25)
plt.axhline(y=50)
plt.axis([0, 11000, 0, 100])
plt.show()
