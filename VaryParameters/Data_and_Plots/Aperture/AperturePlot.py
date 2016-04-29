from __future__ import division
import sys
sys.path.append('../')
import functions as fc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

dir='Data_and_Plots/Aperture/'
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
ApertureDiffs=[]
DiffList=[]

for indx,ApertureSize in enumerate(fc.ApertureRange):
    
    ParentFolderName='Aperture=%1.2f' %ApertureSize
    
    AfileName='Avals/Adata(Aperture=%1.2f).txt' %(ApertureSize)
    BfileName='BVals/Bdata(Aperture=%1.2f).txt' %(ApertureSize)
    
    
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
        ApertureDiffs.append(ApertureSize)

    Arangelist.append(Arange)
    Brangelist.append(Brange)
    RatioList.append((Brange/Arange)**2)
    Apperture_List.append(ApertureSize)

plt.figure(1)
plt.plot(Apperture_List,Brangelist,'.')
plt.figure(2)
plt.plot(Apperture_List,RatioList,'.')
plt.figure(3)
plt.plot(Apperture_List,Arangelist,'.')


plt.figure(1)
plt.xlabel('Aperture Size(cm)')
plt.ylabel('Range of B')
plt.xlim([1.5,5])
plt.savefig(dir+'BRange.eps')


plt.figure(3)
plt.xlabel('Aperture Size(cm)')
plt.ylabel('Range of A')
plt.xlim([1.5,5])
plt.savefig(dir+'ARange.eps')


plt.figure(4)
plt.plot(ApertureDiffs,DiffList,'.')
plt.xlim([1.5,5])
plt.xlabel('Aperture Size (cm)')
plt.ylabel(r'Difference between $A_0$ $B_0$ (keV)')
plt.savefig(dir+'A0 B0 difference.eps')

#plt.show()