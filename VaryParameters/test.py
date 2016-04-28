import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

z=0.05
T=10
ParentFolder='Pos=25'
fc.GPTrun(ParentFolder,'A',ScreenPos=50,InitialZ=z)
fc.GPTrun(ParentFolder,'B',ScreenPos=50,InitialZ=z)

AfileName='AVals/Adata(%s).txt' %ParentFolder
BfileName='BVals/Bdata(%s).txt' %ParentFolder
fc.ABSet(AfileName,BfileName)
fc.ABCheck()
Tfit,eTfit=fc.GPTrun(ParentFolder,'T',ScreenPos=50,Temp=T,N=2000,InitialZ=z)

#
#
#ParentFolder='Aperture=5.00'
#
#
#zl=[0.025,0.05,0.075]
#T=10
#fittedT=[]
#efittedT=[]
#
#for z in zl:
#    
#    fc.GPTrun(ParentFolder,'A',InitialZ=z)
#    fc.GPTrun(ParentFolder,'B',InitialZ=z)
#
#    AfileName='AVals/Adata(%s).txt' %ParentFolder
#    BfileName='BVals/Bdata(%s).txt' %ParentFolder
#
#    fc.ABSet(AfileName,BfileName)
#    Tfit,eTfit=fc.GPTrun(ParentFolder,'T',Temp=T,N=2000,InitialZ=z)
#    
#    fittedT.append(Tfit)
#    efittedT.append(eTfit)
#
#fc.Saver('InitialPos.txt',zl,fittedT,efittedT)
#
plt.figure(25)
plt.xlabel('Kinetic Energy after acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.show()
#plt.savefig('Diri.eps')
#
#plt.figure(1)
#plt.errorbar(zl,fittedT,efittedT,fmt='.')
#
#
#plt.show()

