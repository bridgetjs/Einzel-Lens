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
loopmin=1
loopmax=14
A=1
Aplot=1
Eplot=0
oldstdout = sys.stdout
ParentFolderName="AcceleratorPlates"
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
A_List=[]
V=[]
UI=[]
UF=[]
sigma=[]
sigma2=[]
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'

Afit,Bfit=fc.ABFIT()


e=1.602e-19
sigmax=2e-3

def Model(U,T):
    sigma=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
    return sigma

def PlotModel(Energy,T):
    Sigmaarray=[]
    for U in Energy:
        sigmaU=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
        Sigmaarray.append(sigmaU)
    return Sigmaarray


def rms(vec):
    ssum=0
    for i in vec:
        ssum+=i**2
    return math.sqrt(ssum/len(vec))

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
    Tbunch=10
    GBwidth=math.sqrt(1.381e-23*Tbunch/9.11e-31)/3e8
    beamdyn=open("beamdynTsim.in",'w')
    sys.stdout = beamdyn
    print('#Beam Dynamics File')
    print(' ') ;
    print ('Eo   =1e-4;                 # Energy [eV]')
    print('G    = 1-qe*Eo/(me*c*c) ;	# Corresponding Lorentz factor G []')
    print('Beta = sqrt(1-G^-2) ;		# Corrsponding Normalized velocity []')
    print('dE   = 1e-10 ;              # Energy spread [fraction of Eo]')
    print('I    = 0.02 ;               # Beam current [A]')
    print('Tbunch=100;')
    print('GBwidth=G*sqrt(1.381e-23*Tbunch/(me*c*c));')

    print('zlen = 0.002; #bunch length')
    print('radius = 0.002; #bunch radius')

    print('nps  = 100;		        # Number of particles []')
    print('Qtot = -1.e-12 ;		# Total charge in bunch [C]')
    print('nmac = (Qtot/qe)/nps ;		# Initial nmacro')

    print('setparticles("beam",nps,me,qe,Qtot) ;')
    print('setzdist("beam","g",0, zlen, 3, 3) ;')
    print('setrxydist("beam","g", 0, radius, 0, 3) ;')
    print('setphidist("beam","u",0,2*pi) ;')
    print('setGBzdist( "beam", "g", 0, GBwidth, 3, 3 ) ;')
    print('setGBrxydist( "beam", "g", 0, GBwidth, 0, 3) ;')
    print('setGBphidist("beam","u",0,2*pi) ;')
    print('setGBxemittance( "beam", 2e-9 ) ;')
    print('setGByemittance( "beam", 2e-9 ) ;')
    print('setGdist("beam","u", G, dE*(G-1) ) ;')

    print('map2D_E("wcs","z",0,"fieldmap.gdf","R","Z","Er","Ez",1) ;')
    print('Zcol=0;')
    print('copperscatter("wcs","z",0, "copper",0.01*nmac,2*nmac) ;')

    print('accuracy(5) ; ')
    print('screen("wcs","I",.9);')
    print('screen("wcs","I",0.05);')
    #    print('snapshot(0,2e-7,5e-11);')
    beamdyn.close();
    sys.stdout=oldstdout


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
    print "Running GPT in "+ FolderName

    os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100\
              \n gpt_track -o result.gdf beamdynTsim.in boundaries1.in  GPTLICENSE=1539721513 \
              \n gdftrans -o std1.gdf result.gdf position x y rxy z Bx By Bz G \
              \n gdf2a -o std1.txt std1.gdf")


    input_file='std1.txt'
    
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
    initialbx=[]
    for l in fin:
        
        
        if(len(l.split())==0):
            
            
            
            inx=float(thispartx[0])
            finx=float(thispartx[1])
            inxarray.append(inx)
            finxarray.append(finx)
            
            gammainitial=thispartGamma[0]
            gammafinal=thispartGamma[1]
            Tfinal=thispartT[1]
            Tinitial=thispartT[0]
            finalT.append(Tfinal)
            initialT.append(Tinitial)
            #
            #            print Tfinal,Tinitial
            initialbx.append(float(thispartbx[0]))
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
            thispartbx.append(bx)
            thispartby.append(by)
            thispartbz.append(bz)

    UI.append(np.mean(initialT))
    print np.std(initialbx)
    sigma.append(rms(finxarray))

    sigma2=np.std(finxarray)

    os.chdir("../../TSim");

#print sigma

sigmafile=open("OldSigma.txt",'w')
S= UI , sigma
np.savetxt(sigmafile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)

plt.figure(i)
plt.plot(np.asarray(UI),np.asarray(sigma),'k.')
plt.plot(UI,PlotModel(UI,10))
Tfit=optimization.curve_fit(Model, UI, sigma, maxfev=10000)

print Tfit[0]
plt.plot(UI,PlotModel(UI,Tfit[0]))
plt.show()
