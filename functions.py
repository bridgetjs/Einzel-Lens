'''
    Python Module for Development of Einzel Lens using GPT
    
'''
from __future__ import division
from matplotlib import cm
import numpy as np
import os
import math
import sys
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import numpy.random as npr
import matplotlib.pyplot as plt
import scipy.stats as st

Afile=''
Bfile=''

oldstdout = sys.stdout #Save old output
sigmax=2e-3 #Size of Intial Beam
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
PlotterCounter=1
GrandFolder="LensPosVar"
Pathname="SF_Files/"+GrandFolder
LensZRange=[20,30,40,50,60,70,80,90]

def data(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]  #read data assuming data is in the exepcted format
    return x,y

def ABSet(Afilename,Bfilename): #Set file names as global varriables
    global Afile,Bfile
    Afile=Afilename
    Bfile=Bfilename

def VarTest(): #Some random test code
    plt.figure(PlotterCounter)
    plt.plot(npr.rand(5,1),npr.rand(5,1))
    plt.show()
    print PlotterCounter
    print sigmax


def IsMax(Sample,Dist,N ): #Function to return an array of N values V that fall under probablity distribution dist
    V=[]
    P=[]
    PeakMax=max(Dist) #Find max of distribution
    RandomProbs=PeakMax*npr.rand(N ,1) #N random probabilities between 0 and PeakMax
    for indx,Vel in enumerate(Sample): #Loop through Sample
        #If the random probablity is below the distribution, keep the V correspodning to that index
        if (RandomProbs[indx]<Dist[indx]):
            V.append(Vel)
            P.append(RandomProbs[indx])

    return V,P

def Maxwell(N,T): #Produce a Boltzmann distribution
#    Import constants from scipy
    kb=sc.k ; me=sc.m_e
#   Calculate Normalisation and exponential factors
    Nfac=math.sqrt(me/(2*math.pi*kb*T))
    Expofac=me/(kb*T)

#   Produce random velocities within a 5 sigma range
    VSpread=5*math.sqrt(1/Expofac)
    VSampleX=npr.randint(-VSpread,VSpread,N)
    VSampleY=npr.randint(-VSpread,VSpread,N)
    VSampleZ=npr.randint(-VSpread,VSpread,N)


    DistX=[]
    DistY=[]
    DistZ=[]
    for V in VSampleX:
        DistX.append(Nfac*math.exp(-(Expofac*V**2)/2))
    for V in VSampleY:
        DistY.append(Nfac*math.exp(-(Expofac*V**2)/2))
    for V in VSampleZ:
        DistZ.append(Nfac*math.exp(-(Expofac*V**2)/2))

    #   Call is max function to get a boltzmann distribution
    Vx,Px=IsMax(VSampleX,DistX,N)
    Vy,Py=IsMax(VSampleY,DistY,N)
    Vz,Pz=IsMax(VSampleZ,DistZ,N)

    S=[Vx,Px,Vy,Py,Vz,Pz]
    return Vx,Px,Vy,Py,Vz,Pz

def BeamGenerator(N,T,sigma,centreZ):
    
    Ntot=10*N
    c=sc.c
    Vx,Px,Vy,Py,Vz,Pz=Maxwell(Ntot,T)
    

#    print Vx
    Bx=[v/c for v in Vx[0:N]]
    By=[v/c for v in Vy[0:N]]
    Bz=[v/c for v in Vz[0:N]]

 
    x=st.norm.rvs(loc=0, scale=sigma, size=N)
    y=st.norm.rvs(loc=0, scale=sigma, size=N)
    z=st.norm.rvs(loc=centreZ, scale=sigma, size=N)
    beam=[x,y,z,Bx,By,Bz]
    
    return beam

def ABFIT():
    
    U,A=data(Afile)
    Afit=si.UnivariateSpline(U,A)
    
    U,B=data(Bfile)
    Bfit=si.UnivariateSpline(U,B)
    
    return Afit,Bfit

def Model2(U,T):
    
    e=sc.e
#    T=abs(T)
    Afit,Bfit=ABFIT()

    sigma=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
    return sigma

def PlotModel(Energy,T,Afit,Bfit,sigmax):
    Sigmaarray=[]
    e=sc.e
    T=abs(T)
    for U in Energy:
        sigmaU=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
        Sigmaarray.append(sigmaU)
    return Sigmaarray

def rms(vec):
    ssum=0
    for i in vec:
        ssum+=i**2
    return math.sqrt(ssum/len(vec))

def BeamDynWriter(Name,Option,Outputstyle,beam,MCPPos=0.95,*args):
    oldstdout = sys.stdout
    beamdyn=open(Name,'w')
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

    print('nps  = 200;		        # Number of particles []')
    print('Qtot = -1.e-12 ;		# Total charge in bunch [C]')
    print('nmac = (Qtot/qe)/nps ;		# Initial nmacro')

    if Option=="B":
        for j in range(1,2000):
            xmom=-0.5e-4+j*(1e-4/2000)
            print 'setstartpar("beam",0,0,0.025,%1.3e,0,0,me,qe,1);' %xmom
    if Option=="A":
        print('setstartline("beam",nps,-radius,0 ,0.025,radius,0,0.025,0,0,0,me,qe,1);')
    if Option=="GenBunch":
#        Call Maxwell Generation?
        print('Tbunch=10;')
        print('GBwidth=G*sqrt(1.381e-23*Tbunch/(me*c*c));')
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
    if Option=="beam":
        
        for j in range(0,len(beam[0])-1):
            print 'setstartpar("beam",%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,%1.3e,me,qe,1);' %(beam[0][j] ,beam[1][j] ,beam[2][j] ,beam[3][j] ,beam[4][j] ,beam[5][j])


    print('map2D_E("wcs","z",0,"fieldmap.gdf","R","Z","Er","Ez",1) ;')
    print('Zcol=0;')
    print('copperscatter("wcs","z",0, "copper",0.01*nmac,2*nmac) ;')

    print('accuracy(5) ; ')

    if Outputstyle=="Screens":
        print 'screen("wcs","I",%1.3f);' %MCPPos
        print('screen("wcs","I",0.1);')
    if Outputstyle=="Snaps":
        print('snapshot(0,1e-7,5e-11);')
    if Outputstyle=="Touts":
        print ('tout(0,5e-8,5e-11) ;')

    beamdyn.close();
    sys.stdout=oldstdout

def Fisher():
    
    os.system("fishfile -o fieldmap.gdf -g boundaries.in EinzelLens.am")
    file=open("boundaries.in",'r')
    outfile=open("boundaries1.in",'w')
    lines = file.readlines()
    fin = lines[:-5]
    outfile.writelines(fin)
    file.close()
    outfile.close()

def GPTCall(DynamicFile,GroupBy,Outtxt,*args):
    Outstring = 'x y rxy z Bx By Bz G'
    for arg in args:
        if arg=="Efields":
            Outstring = "x y rxy z Bx By Bz G fEx fEy fEz"

    os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100")
    os.system('gpt_track -o result.gdf %s boundaries1.in  GPTLICENSE=1539721513' %DynamicFile)
    os.system('gdftrans -o Out.gdf result.gdf %s %s ' %(GroupBy, Outstring) )
    os.system('gdf2a -o %s Out.gdf' %Outtxt)


def line(x, a, b):
    return a*x + b

def Plotter(Infile,FolderName,i,*args):
    global PlotterCounter
    fin=open(Infile,'r')
    
    thispartx=[]; thisparty=[]; thispartrxy=[]; thispartz=[]; thisparttime=[];
    thispartGamma=[]; thispartG=[]; thispartT=[];
    
    thispartbx=[]; thispartby=[]; thispartbz=[];

    inxarray=[] ; finxarray=[];
    
    initialT=[] ; finalT=[] ;
    
    div0array=[] ;
    
    initialbx=[]
    
    Ex=[]; Ey=[]; Ez=[];
    newpart=0
    
    for indx,ar in enumerate(args):
        if ar=='xz':
            xzfig=PlotterCounter
            PlotterCounter+=1
        if ar=='yz':
            yzfig=PlotterCounter
            PlotterCounter+=1
        if ar=='Tz':
            TZfig=PlotterCounter
            PlotterCounter+=1
        if ar=='Efields':
            Efig=PlotterCounter
            PlotterCounter+=1

    for l in fin:
        if(len(l.split())==0):
            
#           Plotting Code here:
            for indx,ar in enumerate(args):
                
                if ar=='xz':
                    plt.figure(xzfig)
                    plt.plot(thispartz,thispartx)
                if ar=='yz':
                    plt.figure(yzfig)
                    plt.plot(thispartz,thisparty)
                if ar=='inx':
                    inx=float(thispartx[0])
                    inxarray.append(inx)
                    Tinitial=thispartT[0]
                    gammainitial=thispartGamma[0]
                    initialT.append(Tinitial)
                if ar=='finx':
                    finx=float(thispartx[1])
                    finxarray.append(finx)
                    Tfinal=thispartT[1]
                    finalT.append(Tfinal)
                    gammafinal=thispartGamma[1]
                if ar=='div0':
                    div0=float(thispartbx[0])/float(thispartbz[0])
                    div0array.append(div0)
                if ar=='Tz':
                    plt.figure(TZfig)
                    plt.plot(thispartz,thispartT)
                if ar=='Efields':
                    plt.figure(Efig)
                    plt.plot(thispartz,Ex)
        
        
            thispartx=[]; thisparty=[]; thispartrxy=[]; thispartz=[]; thisparttime=[];
            thispartGamma=[]; thispartG=[]; thispartT=[];
            
            thispartbx=[]; thispartby=[]; thispartbz=[];
            
            Ex=[]; Ey=[]; Ez=[];
          
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
        
        
            thispartx.append(x);  thisparty.append(y);
            thispartrxy.append(rxy) ; thispartz.append(z)
            
            thispartGamma.append(gamma)
            thispartT.append(T)
            thispartG.append(G)

            thispartbx.append(bx);  thispartby.append(by);   thispartbz.append(bz)

            for ar in args:
                if ar=='Efields':
                    Ex.append(l.split()[8])
                    Ey.append( l.split()[9])
                    Ez.append( l.split()[10])

    for indx,ar in enumerate(args):
        if ar=='xz':
            plt.figure(xzfig)
            plt.xlabel('z (m)')
            plt.ylabel('x (m)')
        if ar=='yz':
            plt.figure(yzfig)
            plt.xlabel('z (m)')
            plt.ylabel('y (m)')
        if ar=='rms':
            rmsx=rms(finxarray)
            return rmsx
        if ar=='B':
            plt.figure(PlotterCounter)
            PlotterCounter+=1
            plt.plot(div0array,finxarray,'r.')
            plt.xlabel('initial divergence (m/rad)')
            plt.ylabel('final x (m) ')
            a0=([0.0, 0.0])
            fit = optimization.curve_fit(line, div0array, finxarray, a0)
            Bvalue=fit[0][0]
            plt.plot(np.asarray(div0array),Bvalue*np.asarray(div0array))
            plotname="../../../../GPT_B/Archive/"+FolderName
            plt.savefig(plotname+ ".eps")
            plt.close()
            return np.mean(initialT),Bvalue
        if ar=='A':
            plt.figure(PlotterCounter)
            PlotterCounter+=1
            plt.plot(inxarray,finxarray,'r.')
            plt.xlabel('initial x (m)')
            plt.ylabel('final x (m) ')
            a0=([0.0, 0.0])
           
            fit = optimization.curve_fit(line, inxarray, finxarray,a0)
            Avalue=fit[0][0]
            plt.plot(np.asarray(inxarray),Avalue*np.asarray(inxarray))
            plotname="../../../../GPT_A/Archive/"+FolderName
            plt.savefig(plotname+ ".eps")
            return np.mean(initialT),Avalue
        if ar=='Tz':
            plt.figure(TZfig)
            plt.plot()
            plt.xlabel('z (m)')
            plt.ylabel('Kinetic Energy (keV)')
        if ar=='RealBunch':
            return finxarray,np.mean(initialT);
        if ar=='Efields':
            plt.figure(Efig)
            plt.xlabel('z (m)')
            plt.ylabel('E_x (V/m)')
def EinzelLensGen(Var,Range):
    return 0

def Eplotter():
    infile='SHORTOUTSF7.TXT'



