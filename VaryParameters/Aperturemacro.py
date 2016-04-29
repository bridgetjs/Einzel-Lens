import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import time
import os

Tlist=[10,20,30]
fittedT=[]
efittedT=[]
Aperture=[]
V=[]
dir='Data_and_Plots/Aperture/'
if not os.path.exists(dir):
    os.mkdirs(dir)

for ApertureSize in fc.ApertureRange:
    
    ParentFolderName='Aperture=%1.2f' %ApertureSize

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
        Aperture.append(ApertureSize)

fc.Saver(dir+'Data.txt',Aperture,fittedT,efittedT)

plt.figure(1)
plt.errorbar(Aperture,fittedT,efittedT,fmt='.')
plt.axis([1.5,5,0,40])
for T in Tlist: plt.axhline(y=T)
plt.xlabel('Aperture size (cm)')
plt.ylabel('Fitted temperature (K)')
plt.savefig(dir+'Tfits.eps')

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig(dir+'Curves.eps')
plt.show()