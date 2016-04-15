from __future__ import print_function
from matplotlib import cm
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import os
import functions as fc
loopmin=1
loopmax=21
oldstdout = sys.stdout


####### TO VARY
ApertureSize=2
V2=4000;
V1=0;
LensZRange=fc.LensZRange
PlateSept=1; # Z seperation between the ring
FlightTubeLength= 100; # length of the Cavity (cm)

GrandFolder=fc.GrandFolder
Pathname="SF_Files/"+GrandFolder

#if not os.path.exists(Pathname):
#    os.makedirs(Pathname)
#os.chdir(Pathname)

LensZ=50
ApertureSizeRange=fc.ApertureRange

###############################################  Prepare OUTSF7 Files for GPT    ###############################################################
for PlateSept in fc.SeptRange:
    ParentFolderName='Sept=%1.1f' %PlateSept
    
    for i in range(loopmin,loopmax):
        
        VBackPlate=0-250*i;
        
        
        #Print a string with the Voltages of the Plates
        FolderName2='VBack=%d' %VBackPlate
        
        #Make and change to the corresponding directory
        infile=Pathname + "/"+ParentFolderName+"/"+FolderName2+"/"+"OUTSF7.TXT" # Infile path
        outfile=Pathname + "/"+ParentFolderName+"/"+FolderName2+"/"+"SHORTOUTSF7.TXT" # Outfile path
        numline=34 # skip the header
        print('Reading from', infile,' Processing for GPT input') #
        #Header line
        p="         R             Z                Er         Ez             absE             V\n"

        headerfinish=False

        fin=np.loadtxt(infile,skiprows=numline)
        np.savetxt(outfile, fin, fmt='%1.5e', delimiter='\t', header = '\tR\t\tZ\t\tEr\t\tEz\t\tabsE\t\tV', comments='')

sys.stdout=oldstdout
