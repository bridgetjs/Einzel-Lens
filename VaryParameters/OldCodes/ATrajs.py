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
import functions as fc

loopmin=14
loopmax=17

oldstdout = sys.stdout
ParentFolderName="AcceleratorPlates"

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
    
    fc.BeamDynWriter("beamdynAfcTest.in","A","Touts",0)

    fc.Fisher()
    
    print "Running GPT in "+ FolderName
        #           DynamicFile,     GroupBy,    Outtxt
    fc.GPTCall("beamdynAfcTest.in","time","std1.txt","Efields")
        
    fc.Plotter("std1.txt",FolderName,i,'xz')


    os.chdir("../../TSim");

plt.show()
