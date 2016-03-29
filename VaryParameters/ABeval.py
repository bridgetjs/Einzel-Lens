from __future__ import division
import sys
sys.path.append('../')
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
    return D

for indx,ScreenPos in enumerate(ScreenPosRange):
    Arangelist=[]
    Brangelist=[]
    RatioList=[]
    Zlist=[]
    for indx2,LensZ in enumerate(LensZRange):
        
        if ScreenPos<=LensZ: continue
        
        ParentFolderName='Pos=%d' %LensZ
        
        AfileName='AVals/Adata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
        BfileName='BVals/Bdata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
        
        
        A=data(AfileName)
        B=data(BfileName)
        
        
        Arange=np.max(A)-np.min(A)
        Brange=np.max(B)-np.min(B)

        Arangelist.append(Arange)
        Brangelist.append(Brange)
        RatioList.append((Brange/Arange)**2)
        Zlist.append(LensZ)
    #    print LensZ,Brange,Arange
#        if indx2==0:
#            plt.figure(1)
#            plt.plot(1e-2*LensZ,Brange,fc.colourstring[indx]+'.',label='Pos=%d' %ScreenPos)
#
#            plt.figure(2)
#            plt.plot(1e-2*LensZ,(Brange/Arange)**2,fc.colourstring[indx]+'.',label='Pos=%d' %ScreenPos)
#
#            plt.figure(3)
#            plt.plot(1e-2*LensZ,Arange,fc.colourstring[indx]+'.',label='Pos=%d' %ScreenPos)
#        else:
#            plt.figure(1)
#            plt.plot(1e-2*LensZ,Brange,fc.colourstring[indx]+'.')
#
#            plt.figure(2)
#            plt.plot(1e-2*LensZ,(Brange/Arange)**2,fc.colourstring[indx]+'.')
#
#            plt.figure(3)
#            plt.plot(1e-2*LensZ,Arange,fc.colourstring[indx]+'.')

    plt.figure(1)
    plt.plot(Zlist,Brangelist,fc.colourstring[indx],label='Pos=%d' %ScreenPos)
    plt.figure(2)
    plt.plot(Zlist,RatioList,fc.colourstring[indx],label='Pos=%d' %ScreenPos)
    plt.figure(3)
    plt.plot(Zlist,Arangelist,fc.colourstring[indx],label='Pos=%d' %ScreenPos)


plt.figure(1)
plt.legend()
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of B')

plt.figure(2)
plt.xlabel('Lens Position(m)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2$')
plt.legend()

plt.figure(3)
plt.xlabel('Lens Position(m)')
plt.ylabel('Range of A')
plt.legend()


plt.show()