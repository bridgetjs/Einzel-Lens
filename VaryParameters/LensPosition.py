import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

Tlist=[10]
fittedT=[]
efittedT=[]
Aperture=[]
Z=[]
Screen=[]

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
        #    fc.ABCheck()
        
#        for T in Tlist:
#            Tfit,eTfit=fc.GPTrun(ParentFolderName,'T',Temp=T,N=2000,ScreenPos=ScreenZ,InitialZ=0.025,BunchSize=2e-3)
#            fittedT.append(Tfit)
#            efittedT.append(eTfit)
#            Z.append(LensZ)
#            Screen.append(ScreenZ)

fc.Saver('Data_and_Plots/Position/Data.txt',Z,Screen,fittedT,efittedT)

plt.figure(25)
plt.xlabel('Kinetic Energy after 7.5 cm of acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig('Data_and_Plots/Position/Curves.eps')


for indx,Screens in enumerate(Screen):
    
    if Z[indx]==15:
        plt.figure(1)
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[0])
    if Z[indx]==20:
        plt.figure(1)
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[1])
    
    if Z[indx]==25:
        plt.figure(1)
        
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[2])
    if Z[indx]==30:
        plt.figure(1)
        
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[3])
    if Z[indx]==35:
        plt.figure(1)
        
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[4])
    if Z[indx]==40:
        plt.figure(1)
        plt.errorbar(Screen[indx],fittedT[indx],efittedT[indx],fmt='.'+colourstring[5])


plt.figure(1)
plt.plot(0,0,'.'+colourstring[0],label='Lens=15cm')
plt.plot(0,0,'.'+colourstring[1],label='Lens=20cm')
plt.plot(0,0,'.'+colourstring[2],label='Lens=25cm')
plt.plot(0,0,'.'+colourstring[3],label='Lens=30cm')
plt.plot(0,0,'.'+colourstring[4],label='Lens=35cm')
plt.plot(0,0,'.'+colourstring[5],label='Lens=40cm')

for T in Tlist: plt.axhline(y=T)
plt.axis([10, 55, 0, 20])
plt.xlabel('Screen Position (cm)')
plt.ylabel('Fitted Temperature (K)')


