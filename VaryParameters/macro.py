import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import time

Tlist=[10,20,30]
fittedT=[]
efittedT=[]
Aperture=[]
for ApertureSize in fc.ApertureRange:
    
    ParentFolderName='Aperture=%1.2f' %ApertureSize

    fc.GPTrun(ParentFolderName,'A')
    fc.GPTrun(ParentFolderName,'B')
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName

    fc.ABSet(AfileName,BfileName)
    Afit,Bfit=fc.ABFIT()
    for T in Tlist:
        Tfit,eTfit=fc.GPTrun(ParentFolderName,'T',Temp=T,N=2000)
        fittedT.append(Tfit)
        efittedT.append(eTfit)
        Aperture.append(ApertureSize)

fc.Saver('Data and Plots/April 23rd Aperture/ApertureDataSept=1.txt',Aperture,fittedT,efittedT)

plt.figure(1)
plt.errorbar(Aperture,fittedT,efittedT,fmt='.')
for T in Tlist: plt.axhline(y=T)
plt.xlabel('ApertureSize (cm)')
plt.ylabel('Fitted Temperature (K)')

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.show()