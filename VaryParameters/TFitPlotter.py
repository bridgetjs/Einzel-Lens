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

def data4(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2]; w=data[:,3]   #read data assuming data is in the exepcted format
    return x,y,z,w

datefile='Screen&Lens.txt'
#datefile='Tin=25.txt'
Lenses,Screens,T,eT=data4(datefile)
i=0

for indx,Screen in enumerate(Screens):
    if eT[indx]<4: print Lenses[indx],Screens[indx],T[indx],'+/-',eT[indx]



for indx,Screen in enumerate(Screens):
   
    if Lenses[indx]==30:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[0])
    if Lenses[indx]==40:
        plt.figure(1)
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[1])
    
    if Lenses[indx]==50:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[2])
    if Lenses[indx]==60:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[3])
    if Lenses[indx]==70:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[4])
    if Lenses[indx]==80:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[5])
    if Lenses[indx]==90:
        plt.figure(1)
        
        plt.errorbar(Screens[indx],T[indx],eT[indx],fmt='.'+colourstring[6])

plt.figure(1)
plt.plot(0,0,'.'+colourstring[0],label='Lens=30cm')
plt.plot(0,0,'.'+colourstring[1],label='Lens=40cm')
plt.plot(0,0,'.'+colourstring[2],label='Lens=50cm')
plt.plot(0,0,'.'+colourstring[3],label='Lens=60cm')
plt.plot(0,0,'.'+colourstring[4],label='Lens=70cm')
plt.plot(0,0,'.'+colourstring[5],label='Lens=80cm')
plt.plot(0,0,'.'+colourstring[6],label='Lens=90cm')
plt.axhline(y=25)
plt.axis([20, 110, 0, 100])
plt.xlabel('Screen Position (cm)')
plt.ylabel('Fitted Temperature (K)')
plt.legend(loc=2)
#plt.savefig('ScreenPosition.eps')

for indx,Lens in enumerate(Lenses):
    
    if Screens[indx]==30:
        plt.figure(2)
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[0])
    if Screens[indx]==40:
        plt.figure(2)
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[1])
    
    if Screens[indx]==50:
        plt.figure(2)
        
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[2])
    if Screens[indx]==60:
        plt.figure(2)
        
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[3])
    if Screens[indx]==70:
        plt.figure(2)

        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[4])
    if Screens[indx]==80:
        plt.figure(2)
        
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[5])
    if Screens[indx]==90:
        plt.figure(2)
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[6])
    if Screens[indx]==100:
        plt.figure(2)
        plt.errorbar(Lenses[indx],T[indx],eT[indx],fmt='.'+colourstring[7])

plt.figure(2)
plt.plot(0,0,'.'+colourstring[0],label='Screen=30cm')
plt.plot(0,0,'.'+colourstring[1],label='Screen=40cm')
plt.plot(0,0,'.'+colourstring[2],label='Screen=50cm')
plt.plot(0,0,'.'+colourstring[3],label='Screen=60cm')
plt.plot(0,0,'.'+colourstring[4],label='Screen=70cm')
plt.plot(0,0,'.'+colourstring[5],label='Screen=80cm')
plt.plot(0,0,'.'+colourstring[6],label='Screen=90cm')
plt.plot(0,0,'.'+colourstring[7],label='Screen=100cm')
plt.axhline(y=25)
plt.axis([10, 100, 0, 100])
plt.xlabel('Lens Position (cm)')
plt.ylabel('Fitted Temperature (K)')
plt.legend(loc=2)

plt.show()
#plt.savefig('ScreenPosition.eps')

