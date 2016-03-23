from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functionsSR as fc

loopmin=1
loopmax=21

oldstdout = sys.stdout
ParentFolderName="Acc2"

TempRange=[25]

InitialSize=fc.sigmax
A_list=[]
UI=[]
V=[]
for i in range(loopmin,loopmax):
    sys.stdout=oldstdout
    VBackPlate=0-250*i;
    V.append(-VBackPlate*1e-3)
    FolderName='VBack=%d' %VBackPlate
    
    #Make and change to the corresponding directory
    Path="../" + ParentFolderName+"/"+FolderName
    if os.path.exists(Path):
        os.chdir(Path);
        print os.getcwd()
    else:
        print "Shit's fucked with ",FolderName
        break
    if not os.path.exists("Plots"):
        os.makedirs("Plots")
    
    fc.BeamDynWriter("beamdynAfcTest.in","A","Screens",0)

    fc.Fisher()
    
    print "Running GPT in "+ FolderName
        #           DynamicFile,     GroupBy,    Outtxt
    fc.GPTCall("beamdynAfcTest.in","position","std1.txt")
        
    U_i,A=fc.Plotter("std1.txt",FolderName,i,'inx','finx','A')
    A_list.append(A)
    UI.append(U_i)

    os.chdir("../../TSim");

os.chdir("../GPT_A")
Afile=open('Adata.txt','w')
S= UI , A_list
np.savetxt(Afile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
Afile.close()
plt.figure(1)
plt.plot(UI,A_list,'b.')
plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
plt.ylabel('A (no units)')
plt.savefig("A_with_T.eps")
