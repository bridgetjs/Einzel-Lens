from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
import scipy.optimize as optimization

loopmin=1
loopmax=2
A=0
B=1
Aplot=0
oldstdout = sys.stdout
ParentFolderName="AcceleratorPlates"
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
A_List=[]
V=[]
for i in range(loopmin,loopmax):
    sys.stdout=oldstdout
    VBackPlate=0-100*i;
    V.append(-VBackPlate*1e-3)
    FolderName='VBack=%d' %VBackPlate
    
    
    #Make and change to the corresponding directory
    Path="../" + ParentFolderName+"/"+FolderName
    if os.path.exists(Path):
        os.chdir(Path);
    else:
        print "Shit's fucked with ",FolderName
        break
    
    fid2=open("Bscan.mr",'w')
    sys.stdout=fid2
    print('xmom -1e-3 1e-3 1e-4')
    fid2.close()
    sys.stdout=oldstdout

    beamdyn=open("Bbeamdyn.in",'w')
    sys.stdout = beamdyn
    print('#Beam Dynamics File')
    print(' ') ;
    print ('Eo   =1e-4;                 # Energy [eV]')
    print('G    = 1-qe*Eo/(me*c*c) ;	# Corresponding Lorentz factor G []')
    print('Beta = sqrt(1-G^-2) ;		# Corrsponding Normalized velocity []')
    print('dE   = 1e-10 ;              # Energy spread [fraction of Eo]')
    print('I    = 0.02 ;               # Beam current [A]')

    print('zlen = 0.002; #bunch length')
    print('radius = 0.002; #bunch radius')

    print('nps  = 25;		        # Number of particles []')
    print('Qtot = -1.e-12 ;		# Total charge in bunch [C]')
    print('nmac = (Qtot/qe)/nps ;		# Initial nmacro')


    if A:
        print('setstartline("beam",nps,-radius,0 ,0.025,radius,0,0.025,0,0,0,me,qe,1);')
    if B:
        print('setstartpar("beam",0,0,0.025,xmom,0,0,me,qe,1);')

#    print('addxdiv("beam",0,div);')
    print('map2D_E("wcs","z",0,"fieldmap.gdf","R","Z","Er","Ez",1) ;')
    print('Zcol=0;')
    print('copperscatter("wcs","z",0, "copper",0.01*nmac,2*nmac) ;')

    print('accuracy(5) ; ')
    #print('screen("wcs","I",.9);')
    #print('screen("wcs","I",0.05);')
    print('snapshot(0,2e-7,5e-11);')
    beamdyn.close();
    sys.stdout=oldstdout

    def func(x, a, b):
        return a*x + b

    if not os.path.exists("Plots"):
        os.makedirs("Plots")

    os.system("fishfile -o fieldmap.gdf -g boundaries.in EinzelLens.am")

    file=open("boundaries.in",'r')
    outfile=open("boundaries1.in",'w')
    lines = file.readlines()
    fin = lines[:-5]
    outfile.writelines(fin)
    file.close()
    outfile.close()
    print "Running GPT"
    os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100\
              \n mr -o result.gdf Bscan.mr gpt_track Bbeamdyn.in boundaries1.in  GPTLICENSE=1539721513 \
              \n gdftrans -o traj.gdf result.gdf time x y rxy z Bx By Bz G \
              \n gdf2a -o traj.txt traj.gdf")


    input_file='traj.txt'

    fin=open(input_file,'r')

    new_part=0


    thispartx=[]
    thisparty=[]
    thispartrxy=[]
    thispartz=[]
    thispartGamma=[]
    thisparttime=[]
    thispartG=[]
    thispartT=[]
    thispartbx=[]
    thispartby=[]
    thispartbz=[]
    thispartwrongT=[]

    inxarray=[]
    finxarray=[]
    initialT=[]
    finalT=[]
    for l in fin:
        
        
        if(len(l.split())==0):
            
            plt.figure(1)
            plt.plot(thispartz,thispartx)
            
#            plt.figure(2)
#            plt.plot(thispartz,thisparty)
#            
#            plt.figure(3)
#            plt.plot(thispartz,thispartGamma)
#            
#            plt.figure(4)
#            plt.plot(thispartz,thispartG)

            
#            inx=float(thispartx[0])
#            finx=float(thispartx[1])
#            inxarray.append(inx*1e3)
#            finxarray.append(finx*1e3)
#            
#            
#            gammainitial=thispartGamma[0]
#            gammafinal=thispartGamma[1]
#            Tfinal=thispartT[1]
#            Tinitial=thispartT[0]
#            finalT.append(Tfinal)
#            initialT.append(Tinitial)
#            
#            print Tfinal,Tinitial

            thispartx=[]
            thisparty=[]
            thispartrxy=[]
            thispartz=[]
            thispartGamma=[]
            thisparttime=[]
            thispartT=[]
            thispartbx=[]
            thispartby=[]
            thispartbz=[]
            thispartwrongT=[]
            thispartG=[]
            new_part=0
            continue
        
        if(l.split()[0] == 'ID'): continue
        if(l.split()[0] == 'xrad'): continue
        if(l.split()[0] == 'x'):
            new_part=1
            continue
        
        if(new_part==1):
            
            x=l.split()[0]
            y=l.split()[1]
            rxy=l.split()[2]
            z=l.split()[3]
            bx=float(l.split()[4])
            by=float(l.split()[5])
            bz=float(l.split()[6])
            G=float(l.split()[7])
            b=math.sqrt(bx**2+by**2+bz**2)
            gamma=1/math.sqrt(1-b**2)
            
            T=511*(gamma-1)
            wrongT=511*(G-1)
            
            thispartx.append(x)
            thisparty.append(y)
            thispartrxy.append(rxy)
            thispartz.append(z)
            
            thispartGamma.append(gamma)
            thispartT.append(T)
            thispartG.append(G)
            thispartwrongT.append(wrongT)
            #thispartbx.append(bx) thispartby.append(by)  thispartbz.append(bz)


    plt.figure(1)
    plt.xlabel('z (m)')
    plt.ylabel('x (m)')
#    
#    plt.figure(2)
#    plt.xlabel('z (m)')
#    plt.ylabel('y (m) ')
#
#
#    plt.figure(3)
#    plt.xlabel('z (m)')
#    plt.ylabel('Gamma ')
#    
#    plt.figure(4)
#    plt.xlabel('z (m)')
#    plt.ylabel('Gamma -wrong')
#
#    plt.figure(i)
#    plt.plot(inxarray,finxarray,'r.')
#    plt.xlabel('initial x (mm)')
#    plt.ylabel('final x (mm) ')
#    a0=([0.0, 0.0])
#    b= optimization.curve_fit(func, inxarray, finxarray, a0)
#    A_List.append(b[0][0])
#    print "The matrix element is",b[0][0]
#    plt.plot(np.asarray(inxarray),b[0][0]*np.asarray(inxarray),colourstring[i])
#    plotname="../../GPT_A/"+FolderName
#    plt.savefig(plotname+ ".eps")
#    #plt.savefig(plotname+ ".png")
#    plt.close()
#
if Aplot:
    

    
    
    Afile=open('Adata.txt','w')
    S= V , A_List
    np.savetxt(Afile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)

    plt.figure(loopmax+1)
    plt.plot(V,A_List,'b.')
    plt.xlabel('Accelerator Plate Voltage (kV)')
    plt.ylabel('A (no units)')
    plt.savefig("A_with_V.eps")
#plt.show()
    os.chdir("../../GPT_A");

plt.show()


