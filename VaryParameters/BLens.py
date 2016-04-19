from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functions as fc

loopmin=1
loopmax=21

oldstdout = sys.stdout
Pathname=fc.Pathname
#ParentFolderName="AcceleratorPlates"
LensZRange=fc.LensZRange
ScreenPosRange=fc.ScreenPosRange
#LensZRange=[50]
#ScreenPosRange=[100]
TempRange=[25]
InitialPos=fc.InitialPos
InitialSize=fc.sigmax
LensZ=50
ScreenPos=100
B_list=[]
UI=[]
V=[]

for k in range(0,1):
    ParentFolderName='Tuned2'
#for PlateSept in fc.TotSeptRange:
#    ParentFolderName='Sept=%1.2f' %PlateSept
#for V2 in fc.VoltageRange:
#    ParentFolderName='Optimised'
#    ParentFolderName='Voltage=%d' %V2

    B_list=[]
    eB_list=[]
    UI=[]
    V=[]
    
    if ScreenPos<=LensZ: continue

    print 'Running BStudy with Lens at %d and Screen at %d' %(LensZ, ScreenPos)

    for i in range(loopmin,loopmax):
        sys.stdout=oldstdout
        VBackPlate=0-250*i;
        V.append(-VBackPlate*1e-3)
        FolderName='VBack=%d' %VBackPlate
        
        #Make and change to the corresponding directory
        Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
        if os.path.exists(Path):
            os.chdir(Path);
#                print os.getcwd()
        else:
            print "Shit's fucked with ",FolderName
            break
        if not os.path.exists("Plots"):
            os.makedirs("Plots")
    
        fc.BeamDynWriter("beamdynB.in","B","Screens",0,ScreenPos,InitialPos)
        
        fc.Fisher()
        
        print "Running GPT in "+ FolderName
        #           DynamicFile,     GroupBy,    Outtxt
        fc.GPTCall("beamdynB.in","position","stdB.txt")
        
        U_i,B,eB=fc.Plotter("stdB.txt",FolderName,i,'inx','div0','finx','B')
        B_list.append(B)
        eB_list.append(eB)
        UI.append(U_i)

        os.chdir("../../../../VaryParameters");
#    print os.getcwd()

    BfileName='BVals/Bdata(%s).txt' %ParentFolderName

    Bfile=open(BfileName,'w')
    S= UI , B_list, eB_list
    np.savetxt(Bfile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
    Bfile.close()
#    plt.figure(1)
#    plt.plot(UI,B_list,'.',label=ParentFolderName)


#plt.figure(1)
#plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
#plt.ylabel('B (m/rad)')
#plt.legend()
#plt.savefig("B_with_T_LensPos.eps")
#plt.show()
