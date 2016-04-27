import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import time

Tlist=[10,20,30]
fittedT=[]
efittedT=[]
Aperture=[]
Range=[2,2.5,3,3.5,4,4.5,5,5.5,6,6.5]
Range=[3,3.5,4,6]
V=[]
#for ApertureSize in fc.ApertureRange:
#    
#    ParentFolderName='Aperture=%1.2f' %ApertureSize
#
for V2 in fc.VoltageRange:
    ParentFolderName='Pos=%d' %LensZ
#    ParentFolderName='V=%d' V2

    fc.GPTrun(ParentFolderName,'A',ScreenPos=50,InitialZ=0.025)
    fc.GPTrun(ParentFolderName,'B',ScreenPos=50,InitialZ=0.025)
    AfileName='AVals/Adata(%s).txt' %ParentFolderName
    BfileName='BVals/Bdata(%s).txt' %ParentFolderName

    fc.ABSet(AfileName,BfileName)
    Afit,Bfit=fc.ABFIT()
#    fc.ABCheck()

    if ScreenPos<=LensZ: continue
    
    for T in Tlist:
        Tfit,eTfit=fc.GPTrun(ParentFolderName,'T',Temp=T,N=2000,ScreenPos=50,InitialZ=0.025,BunchSize=2e-3)
        fittedT.append(Tfit)
        efittedT.append(eTfit)
        Aperture.append(ApertureSize)

fc.Saver('Data_and_Plots/Position/Data.txt',V,fittedT,efittedT)

plt.figure(1)
plt.errorbar(Aperture,fittedT,efittedT,fmt='.')
plt.axis([1.5,7,0,35])
for T in Tlist: plt.axhline(y=T)
plt.xlabel('Aperture size (cm)')
plt.ylabel('Fitted temperature (K)')
plt.savefig('Data_and_Plots/Position/Tfits.eps')

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig('Data_and_Plots/Position/Curves.eps')
plt.show()