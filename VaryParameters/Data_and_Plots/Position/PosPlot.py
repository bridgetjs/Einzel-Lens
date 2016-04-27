import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

colourstring=fc.colourstring
Tlist=[10]
Z,Screen,fittedT,efittedT=np.loadtxt('Data_and_Plots/Position/Data.txt',unpack=True)


for i,Lens in enumerate(fc.LensZRange):
    for indx,ZL in enumerate(Z):
        if ZL==Lens:
            plt.figure(1)
            plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt=colourstring[i])
    plt.figure(1)
    plt.plot(0,0,'.'+colourstring[i],label='Lens=%dcm' %Lens)

plt.figure(1)
for T in Tlist: plt.axhline(y=T)
plt.axis([14, 51, 0, 40])
plt.xlabel('Screen Position (cm)')
plt.ylabel('Fitted Temperature (K)')
plt.legend(loc=2)


Z=[]
Screen=[]

Arangelist=[]
Brangelist=[]
DiffList=[]
Lensdiffs=[]
Screendiffs=[]
Aperturelist=[]
Rlist=[]

for LensZ in fc.LensZRange:
    
    for ScreenZ in fc.ScreenPosRange:
        
        if ScreenZ<=LensZ: continue
        
        ParentFolderName='Pos=%d' %LensZ
        
        
        fc.GPTrun(ParentFolderName,'A',ScreenPos=ScreenZ,InitialZ=0.025)
        fc.GPTrun(ParentFolderName,'B',ScreenPos=ScreenZ,InitialZ=0.025)
        AfileName='AVals/Adata(%s).txt' %(ParentFolderName)
        BfileName='BVals/Bdata(%s).txt' %(ParentFolderName)
        
        fc.ABSet(AfileName,BfileName)
        Afit,Bfit=fc.ABFIT()
        
        UA,A,eA=np.loadtxt(AfileName,unpack=True)
        UB,B,eB=np.loadtxt(BfileName,unpack=True)

        if len(Afit.roots())==1 and len(Bfit.roots())==1:
            DiffList.append((Afit.roots() - Bfit.roots())[0])
            Screendiffs.append(ScreenZ)
        Arange=np.max(A)-np.min(A)
        Brange=np.max(B)-np.min(B)
        
        Arangelist.append(Arange)
        Brangelist.append(Brange)
        Rlist.append((Brange/Arange)**2)
        Z.append(LensZ)
        Screen.append(ScreenZ)

fc.Saver('Ranges.txt',Z,Screen,Arangelist,Brangelist,Rlist,DiffList,Screendiffs)
lab=0

for i,Lens in enumerate(fc.LensZRange):
    for indx,ZL in enumerate(Z):
        if ZL==Lens:
            plt.figure(2)
            plt.plot(Screen[indx],Arangelist[indx],'.'+colourstring[i])
            plt.figure(3)
            plt.plot(Screen[indx],Brangelist[indx],'.'+colourstring[i])
            if indx<len(DiffList):
                plt.figure(4)
                plt.plot(Screendiffs[indx],DiffList[indx],'.'+colourstring[i])
            plt.figure(5)
            plt.plot(Screen[indx],Rlist[indx],'.'+colourstring[i])
    for n in range(2,6):
        plt.figure(n)
        plt.plot(0,0,'.'+colourstring[i],label='Lens=%d cm' %Lens)


plt.figure(2)
plt.xlabel('Screen Position (cm)')
plt.ylabel('Range of A (no units)')
plt.xlim([14,51])
plt.legend(loc=2)

plt.figure(3)
plt.xlabel('Screen Position (cm)')
plt.ylabel('Range of B (m/rad)')
plt.xlim([14,51])
plt.legend(loc=2)

plt.figure(4)
plt.xlabel('Screen Position (cm)')
plt.ylabel(r'Difference between $A_0$ $B_0$ (keV)')
plt.xlim([14,51])
plt.legend(loc=2)

plt.figure(5)
plt.xlabel('Screen Position (cm)')
plt.ylabel(r'$\frac{R(B)}{R(A)}^2}$ (m/rad)')
plt.xlim([14,51])
plt.legend(loc=2)


plt.show()


