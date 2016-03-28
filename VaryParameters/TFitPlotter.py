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
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'

def data4(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2]; w=data[:,3]   #read data assuming data is in the exepcted format
    return x,y,z,w

datefile='Screen&LensTvalues.txt'
Screens,Lenses,T,eT=data4(datefile)
i=0
for indx,Screen in enumerate(Screens):
    
#    if Lenses[indx]==30:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[0])
#    if Lenses[indx]==40:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[1])
    if Lenses[indx]==50:
        plt.figure(1)
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[2])
#    if Lenses[indx]==60:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[3])
#    if Lenses[indx]==70:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[4])
#    if Lenses[indx]==80:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[5])
#    if Lenses[indx]==90:
#        plt.figure(1)
#        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[6])

plt.figure(1)
plt.axhline(y=50)
plt.axis([25, 125, 25, 75])
plt.xlabel('Screen Position ')
plt.ylabel('Fitted Temperature (K)')
plt.show()