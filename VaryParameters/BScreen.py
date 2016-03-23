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
#LensZRange=[50]
TempRange=[25]

InitialSize=fc.sigmax
B_list=[]
UI=[]
V=[]

for LensZ in LensZRange:
    ParentFolderName='Pos=%d' %LensZ
    if LensZ==75.5:
        ParentFolderName='Pos=75Working'
    B_list=[]
    UI=[]
    V=[]
    for i in range(loopmin,loopmax):
        sys.stdout=oldstdout
        VBackPlate=0-250*i;
        V.append(-VBackPlate*1e-3)
        FolderName='VBack=%d' %VBackPlate
        
        #Make and change to the corresponding directory
        Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
        if os.path.exists(Path):
            os.chdir(Path);
            print os.getcwd()
        else:
            print "Shit's fucked with ",FolderName
            break
        if not os.path.exists("Plots"):
            os.makedirs("Plots")
    
        fc.BeamDynWriter("beamdynB.in","B","Screens",0,0.95)
        
        fc.Fisher()
        
        print "Running GPT in "+ FolderName
        #           DynamicFile,     GroupBy,    Outtxt
        fc.GPTCall("beamdynB.in","position","std1.txt")
        
        U_i,B=fc.Plotter("std1.txt",FolderName,i,'inx','div0','finx','B')
        B_list.append(B)
        UI.append(U_i)

        os.chdir("../../../../VaryParameters");
#    print os.getcwd()

    BfileName='BVals/Bdata(%s).txt' %ParentFolderName
    Bfile=open(BfileName,'w')
    S= UI , B_list
    np.savetxt(Bfile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
    Bfile.close()
    plt.figure(1)
    plt.plot(UI,B_list,'.',label=ParentFolderName)


plt.figure(1)
plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
plt.ylabel('B (m/rad)')
plt.legend()
plt.savefig("B_with_T_LensPos.eps")
plt.show()
