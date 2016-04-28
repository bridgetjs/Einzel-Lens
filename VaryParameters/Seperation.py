import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import time

Tlist=[10,30]
fittedT=[]
efittedT=[]
Sept=[]

V=[]
for PlateSept in fc.TotSeptRange:
    ParentFolderName='Sept=%1.2f' %PlateSept


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
        Sept.append(PlateSept)

fc.Saver('Data_and_Plots/Seperation/Data.txt',Sept,fittedT,efittedT)

plt.figure(1)
plt.errorbar(Sept,fittedT,efittedT,fmt='.')
plt.axis([1.5,5,0,40])
for T in Tlist: plt.axhline(y=T)
plt.xlabel('Aperture size (cm)')
plt.ylabel('Fitted temperature (K)')
plt.savefig('Data_and_Plots/Seperation/Tfits.eps')

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig('Data_and_Plots/Seperation/Curves.eps')
plt.show()