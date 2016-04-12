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
Apperture_List=[]

for indx,ApertureSize in enumerate(fc.ApertureRange):
    
    ParentFolderName='Aperture=%1.1f' %ApertureSize
    
    AfileName='../Avals/Adata(Aperture=%1.1f).txt' %(ApertureSize)
    BfileName='../BVals/Bdata(Aperture=%1.1f).txt' %(ApertureSize)
    
    
    A,eA=data(AfileName)
    B,eB=data(BfileName)
    
    
    Arange=np.max(A)-np.min(A)
#    eA=np.sqrt(  eA[A.index(max(A))]**2   + eA[A.index(min(A))]**2    )
    Brange=np.max(B)-np.min(B)

    Arangelist.append(Arange)
    Brangelist.append(Brange)
    RatioList.append((Brange/Arange)**2)
    Apperture_List.append(ApertureSize)

plt.figure(1)
plt.plot(Apperture_List,Brangelist,fc.colourstring[1])
plt.figure(2)
plt.plot(Apperture_List,RatioList,fc.colourstring[2])
plt.figure(3)
plt.plot(Apperture_List,Arangelist,fc.colourstring[3])


plt.figure(1)
plt.legend()
plt.xlabel('Aperture Size(cm)')
plt.ylabel('Range of B')
#plt.savefig('BApp.eps')

plt.figure(2)
plt.xlabel('Aperture Size(cm)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')
plt.savefig('BA_App.eps')

plt.figure(3)
plt.xlabel('Aperture Size(cm)')
plt.ylabel('Range of A')
plt.legend()
plt.savefig('A_apperture.eps')


plt.show()