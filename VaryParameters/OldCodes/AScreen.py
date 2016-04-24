from __future__ import division
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
Pathname="SF_Files/LensPosVar"
ParentFolderName="AcceleratorPlates"
LensZRange=fc.LensZRange
ScreenPosRange=fc.ScreenPosRange
#LensZRange=[50]
#ScreenPosRange=[100]
TempRange=[25]

InitialSize=fc.sigmax
A_list=[]
UI=[]
V=[]

for ScreenPos in ScreenPosRange:

    for LensZ in LensZRange:
        ParentFolderName='Pos=%d' %LensZ
        A_list=[]
        UI=[]
        V=[]
        
        if ScreenPos<=LensZ: continue
        
        print 'Running AStudy with Lens at %d and Screen at %d' %(LensZ, ScreenPos)
        
        for i in range(loopmin,loopmax):
            sys.stdout=oldstdout
            VBackPlate=0-250*i;
            V.append(-VBackPlate*1e-3)
            FolderName='VBack=%d' %VBackPlate
            
            #Make and change to the corresponding directory
            Path="../" +Pathname +"/"+ ParentFolderName+"/"+FolderName
          
            if os.path.exists(Path):
                os.chdir(Path);
            else:
                print "Shit's fucked with ",FolderName
                break
            if not os.path.exists("Plots"):
                os.makedirs("Plots")
            
            fc.BeamDynWriter("beamdynAfcTest.in","A","Screens",0,ScreenPos,fc.InitialPos)

            fc.Fisher()
            
            print "Running GPT in "+ FolderName
                #           DynamicFile,     GroupBy,    Outtxt
            fc.GPTCall("beamdynAfcTest.in","position","std1.txt")
                
            U_i,A,eA=fc.Plotter("std1.txt",FolderName,i,'inx','finx','A')
            A_list.append(A)
            UI.append(U_i)

            os.chdir("../../../../VaryParameters");

    #    os.chdir("AVals")
        AfileName='Avals/Adata(%s,Screen=%d).txt' %(ParentFolderName,ScreenPos)
        Afile=open(AfileName,'w')
        S= UI , A_list
        np.savetxt(Afile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
        Afile.close()
#        plt.figure(1)
#        plt.plot(UI,A_list,'.',label=ParentFolderName)
#        Afit=si.UnivariateSpline(UI,A_list)
#        plt.plot(UI,Afit(UI))
        #    os.chdir("../")

#    plt.figure(1)
#    plt.xlabel('Kinetic Energy after 2.5 cm of acceleration (keV)')
#    plt.ylabel('A (no units)')
#    plt.legend()
#    plt.savefig("A_with_T_(%s,Screen=%d).eps" %(ParentFolderName,ScreenPos))

