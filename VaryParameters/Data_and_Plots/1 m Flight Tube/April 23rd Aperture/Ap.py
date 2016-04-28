import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

import os
print os.getcwd()

#Tlist=[10,20,30]
#A,T,eT=fc.data3('Data_and_Plots/April 23rd Aperture/ApertureData.txt')
#
#plt.figure(5)
#for indx,Av in enumerate(A):
#    plt.errorbar(A[indx],T[indx],eT[indx],fmt='.')
#for T in Tlist: plt.axhline(y=T)
#plt.xlabel('ApertureSize (cm)')
#plt.ylabel('Fitted Temperature (K)')
#plt.axis([1.75,5.25,0,35])


Arangelist=[]
Brangelist=[]
DiffList=[]
Aperturelist=[]
Rlist=[]

for ApertureSize in fc.ApertureRange:
    
    ParentFolderName='Aperture=%1.2f' %ApertureSize
    
    fc.GPTrun(ParentFolderName,'A',InitialZ=0.025)
    fc.GPTrun(ParentFolderName,'B',InitialZ=0.025)
    
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName
    fc.ABSet(AfileName,BfileName)

    UB,B,eB=fc.data3(BfileName)
    UA,A,eA=fc.data3(AfileName)

    Afit,Bfit=fc.ABFIT()

    if len(Afit.roots())==1:
        DiffList.append((Afit.roots() - Bfit.roots())[0])
    
    Arange=np.max(A)-np.min(A)
    Brange=np.max(B)-np.min(B)
        
    Arangelist.append(Arange)
    Brangelist.append(Brange)
    Rlist.append((Brange/Arange)**2)

plt.figure(3)
plt.plot(fc.ApertureRange,Brangelist,'.',label='Range of B')
plt.plot(fc.ApertureRange,Arangelist,'.',label='Range of A')
plt.plot(fc.ApertureRange,Rlist,'.',label=r'$\frac{R(B)}{R(A)}^2}$ ')
plt.xlim([1.75,5.25])
plt.xlabel('Aperture Radius (cm)')
plt.legend()
plt.ylabel('Value')


plt.figure(4)
plt.plot(fc.ApertureRange[0:len(DiffList)],DiffList,'.')
plt.xlim([1.75,5.25])
plt.xlabel('Plate Seperation (cm)')
plt.ylabel(r'Difference between $A_0$ $B_0$ (keV)')

plt.show()
