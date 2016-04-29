import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import time
import os
import numpy as np
Tlist=[10,10,10]
fittedT=[]
efittedT=[]
Aperture=[]
V=[]
dir='Data_and_Plots/Voltage/'
if not os.path.exists(dir):
    os.makedirs(dir)

for V2 in fc.VoltageRange:
    
    ParentFolderName='Voltage=%d' %V2
    fc.GPTrun(ParentFolderName,'A',ScreenPos=50,InitialZ=0.025)
    fc.GPTrun(ParentFolderName,'B',ScreenPos=50,InitialZ=0.025)
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName

    fc.ABSet(AfileName,BfileName)
    Afit,Bfit=fc.ABFIT()
#    fc.ABCheck()

    for T in Tlist:
        Tfit,eTfit=fc.GPTrun(ParentFolderName,'T',Temp=T,N=2000,ScreenPos=50,InitialZ=0.025,BunchSize=2e-3)
        fittedT.append(Tfit)
        efittedT.append(eTfit)
        V.append(V2)

fc.Saver(dir+'Data(3x10).txt',V,fittedT,efittedT)

plt.figure(2)
plt.errorbar(np.asarray(V)*1e-3,fittedT,efittedT,fmt='.')
plt.axis([2.5,5.5,0,40])
for T in Tlist: plt.axhline(y=T)
plt.xlabel('Central Plate Voltage (kV) ')
plt.ylabel('Fitted temperature (K)')
plt.savefig(dir+'Tfit(3x10).eps')

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig(dir+'Curves.eps')
plt.show()