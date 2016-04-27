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
oldstdout = sys.stdout

Afile=''
Bfile=''

oldstdout = sys.stdout #Save old output
sigmax=2e-3 #Size of Intial Beam
colourstring='bgrcmybgrcmybgrcmybgrcmybgrcmybgrcmy'
PlotterCounter=1
#GrandFolder="LensPosVar"
GrandFolder="LensPropVar"
Pathname="SF_Files/"+GrandFolder
LensZRange=[15,20,25,30,35,40]
ScreenPosRange=[20,25,30,35,40,45,50]
Z0=0.025

ApertureRange=[2,2.5,3,3.5,4,4.5,5]
SeptRange=[0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5]
NewSeptRange=[1.2,1.4,1.6,1.8]
TotSeptRange=sorted(SeptRange+NewSeptRange)
VoltageRange=[3000,3500,4000,4500,5000]



def data(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]  #read data assuming data is in the exepcted format
    return x,y

def data3(filename): #function to return the data belonging to an A or B file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2]  #read data assuming data is in the exepcted format
    return x,y,z

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
#    print centreZ
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
    Afit=si.UnivariateSpline(U,A,s=0)
    
    U,B=data(Bfile)
    Bfit=si.UnivariateSpline(U,B,s=0)
    
    return Afit,Bfit

def ABCheck():

    Afit,Bfit=ABFIT()
    UA,A=data(Afile)
    UB,B=data(Bfile)
    
    plt.figure(51)
    plt.plot(UA,Afit(UA))
    plt.plot(UA,A,'.')
    plt.xlabel('U (keV)')
    plt.ylabel('A (no units)')

    plt.figure(52)
    plt.plot(UB,Bfit(UB))
    plt.plot(UB,B,'.')
    plt.xlabel('U (keV)')
    plt.ylabel('B m/rad')
    plt.show()

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
    Afit,Bfit=ABFIT()
    for U in Energy:
        sigmaU=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
        Sigmaarray.append(sigmaU)
    return Sigmaarray

def PlotModel2(Energy,T):
    Sigmaarray=[]
    e=sc.e
    T=abs(T)
    Afit,Bfit=ABFIT()
    for U in Energy:
        sigmaU=np.sqrt( (Afit(U))**2 * (sigmax**2)  +   ( (1.381e-23 * T ) / (2*e*U*1e3)  )  * (Bfit(U))**2   )
        Sigmaarray.append(sigmaU)
    return Sigmaarray



def rms(vec):
    ssum=0
    for i in vec:
        ssum+=i**2
    return math.sqrt(ssum/len(vec))

def BeamDynWriter(Name,Option,Outputstyle,beam=0,MCPPos=50,InitialPos=0.025,Number=2000):
#    print MCPPos,InitialPos
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

    print 'nps  = %d;		        # Number of particles []' %Number
    print('Qtot = -1.e-12 ;		# Total charge in bunch [C]')
    print('nmac = (Qtot/qe)/nps ;		# Initial nmacro')
    print 'InitialZ=%1.3f;' %InitialPos
    if Option=="B":
        for j in range(1,Number):
            xmom=-0.5e-4+j*(1e-4/Number)
            print 'setstartpar("beam",0,0,%1.3f,%1.3e,0,0,me,qe,1);' %(InitialPos,xmom)
    if Option=="A":
        print('setstartline("beam",nps,-radius,0 ,InitialZ,radius,0,InitialZ,0,0,0,me,qe,1);')
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
    if Option=="BNonlin":
        for j in range(1,Number):
            xmom=-0.65e-3+j*(1.3e-3/Number)
            print 'setstartpar("beam",0,0,InitialZ,%1.3e,0,0,me,qe,1);' %xmom
    if Option=="ANonlin":
        print('setstartline("beam",nps,-3*radius,0 ,InitialZ,3*radius,0,InitialZ,0,0,0,me,qe,1);')

    print('map2D_E("wcs","z",0,"fieldmap.gdf","R","Z","Er","Ez",1) ;')
    print('Zcol=0;')
    print('copperscatter("wcs","z",0, "copper",0.01*nmac,2*nmac) ;')

    print('accuracy(5) ; ')

    if Outputstyle=="Screens":
        print 'screen("wcs","I",%f);' %(MCPPos/100)
        print('screen("wcs","I",0.1);')
    if Outputstyle=="Snaps":
        print('snapshot(0,1e-7,5e-11);')
    if Outputstyle=="Touts":
        print ('tout(0,5e-8,5e-11) ;')

    beamdyn.close();
    sys.stdout=oldstdout

def Fisher():
#   Run the gpt fish file function
    os.system("fishfile -o fieldmap.gdf -g boundaries.in EinzelLens.am")
    #Open the file and remove the last five lines : avoid scattering off the ground plate
    file=open("boundaries.in",'r')
    outfile=open("boundaries1.in",'w')
    lines = file.readlines()
    fin = lines[:-5]
    outfile.writelines(fin)
#    Close and delete files
    file.close(); os.remove("boundaries.in")
    outfile.close()

def GPTCall(DynamicFile,GroupBy,Outtxt,*args):
    os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100")
    os.system('gpt_track -o result.gdf %s boundaries1.in  GPTLICENSE=1539721513' %DynamicFile)
    
    Outstring = 'x y rxy z Bx By Bz G'
    Output ='gdftrans -o Out.gdf result.gdf %s x y rxy z Bx By Bz G' %GroupBy
    
    for arg in args:
        if arg=="Efields":
            Outstring = "x y rxy z Bx By Bz G fEx fEy fEz"
            Output ='gdftrans -o Out.gdf result.gdf %s x y rxy z Bx By Bz G fEx fEy fEz' %GroupBy
        if arg=='std':
            Output ='gdfa -o Out.gdf result.gdf %s stdx stdy avgz avgBx avgBy avgBz' %GroupBy
    os.system(Output)
    os.system('gdf2a -o %s Out.gdf' %Outtxt)



def line(x, a, b):
    return a*x + b

def Plotter(Infile,FolderName,*args):
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

    new_part=0
    std_data=0

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
                
            if 'xz' in args:
                plt.figure(xzfig)
                plt.plot(thispartz,thispartx)
            if 'yz' in args:
                plt.figure(yzfig)
                plt.plot(thispartz,thisparty)
            if 'inx' in args:
                inx=float(thispartx[0])
                inxarray.append(inx)
                Tinitial=thispartT[0]
                gammainitial=thispartGamma[0]
                initialT.append(Tinitial)
            if 'finx' in args:
                if len(thispartx)!=2: continue
                finx=float(thispartx[1])
                finxarray.append(finx)
                Tfinal=thispartT[1]
                finalT.append(Tfinal)
                gammafinal=thispartGamma[1]
            if 'div0' in args:
                div0=float(thispartbx[0])/float(thispartbz[0])
                div0array.append(div0)
            if 'Tz' in args:
                plt.figure(TZfig)
                plt.plot(thispartz,thispartT)
            if 'Efields' in args:
                plt.figure(Efig)
                plt.plot(thispartz,Ex)
        
            if 'RealBunchStds' in args:
                ScreenStd=float(thispartx[1])
                InitialKin=float(thispartT[0])
        
            if 'RealBunchShape' in args:
                ScreenStdx=float(thispartx[1])
                InitialKin=float(thispartT[0])
                ScreenStdy=float(thisparty[1])
        
            thispartx=[]; thisparty=[]; thispartrxy=[]; thispartz=[]; thisparttime=[];
            thispartGamma=[]; thispartG=[]; thispartT=[];
            
            thispartbx=[]; thispartby=[]; thispartbz=[];
            
            Ex=[]; Ey=[]; Ez=[];
          
            new_part=0
            std_data=0
            continue
            
        if(l.split()[0] == 'ID'): continue
        if(l.split()[0] == 'x'):
            new_part=1
            continue
        if(l.split()[0] == 'position' or l.split()[0] == 'time'):
            std_data=1
            continue

        if (std_data==1):
            
            stdx=l.split()[1]
            stdy=l.split()[2]
            avgz=l.split()[3]
            avgBx=float(l.split()[4])
            avgBy=float(l.split()[5])
            avgBz=float(l.split()[6])
            b=math.sqrt(avgBx**2+avgBy**2+avgBz**2)
            gamma=1/math.sqrt(1-b**2)
            
            T=511*(gamma-1)

            thispartx.append(stdx);  thisparty.append(stdy);
            thispartz.append(avgz)
            thispartT.append(T)
            thispartbx.append(avgBx);  thispartby.append(avgBy);   thispartbz.append(avgBz)
        
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

            if 'Efields' in args:
                Ex.append(l.split()[8])
                Ey.append( l.split()[9])
                Ez.append( l.split()[10])

    if 'xz' in args:
        plt.figure(xzfig)
        plt.xlabel('z (m)')
        plt.ylabel('x (m)')
    if 'yz' in args:
        plt.figure(yzfig)
        plt.xlabel('z (m)')
        plt.ylabel('y (m)')
    if 'rms' in args:
        rmsx=rms(finxarray)
        return rmsx
    if 'B' in args:
        fit = optimization.curve_fit(line, div0array, finxarray)
        Bvalue=fit[0][0]
        eB=np.sqrt(np.diag(fit[1])[0])
        return np.mean(initialT),Bvalue,eB
        
#        plt.figure(PlotterCounter)
#        PlotterCounter+=1
#        plt.plot(div0array,finxarray,'r.')
#        plt.xlabel('initial divergence (rad)')
#        plt.ylabel('final x (m) ')
#         plt.plot(np.asarray(div0array),Bvalue*np.asarray(div0array))
#        plotname="../../../../GPT_B/Archive/"+FolderName
#        plt.savefig(plotname+ ".eps")
#        plt.close()


    if 'A' in args:
        a0=([0.0, 0.0])
        fit = optimization.curve_fit(line, inxarray, finxarray)
        Avalue=fit[0][0]
        #            print np.sqrt(np.diag(fit[1])[0])
        eA=np.sqrt(np.diag(np.cov(inxarray,finxarray))[0])
        
        return np.mean(initialT),Avalue,eA

#        plt.figure(PlotterCounter)
#        PlotterCounter+=1
#        
#        plt.plot(inxarray,finxarray,'r.')
#        plt.xlabel('initial x (m)')
#        plt.ylabel('final x (m) ')
#        
##            print np.cov(inxarray,finxarray)
#
#        plt.plot(np.asarray(inxarray),Avalue*np.asarray(inxarray))
#        plotname="../../../../GPT_A/Archive/"+FolderName
#        plt.savefig(plotname+ ".eps")
#        plt.close()

    if 'Tz' in args:
        plt.figure(TZfig)
        plt.plot()
        plt.xlabel('z (m)')
        plt.ylabel('Kinetic Energy (keV)')

    if 'RealBunch' in args:
        return finxarray,np.mean(initialT);
    if 'Efields' in args:
        plt.figure(Efig)
        plt.xlabel('z (m)')
        plt.ylabel('E_x (V/m)')
                
    if 'RealBunchStds' in args:
        return ScreenStd,InitialKin
    if 'RealBunchShape' in args:
        return ScreenStdx,ScreenStdy,InitialKin

    if 'avgs' in args:
        plt.figure(xzfig)
        plt.xlabel('z (m)')
        plt.ylabel('stdx (m)')
        
        plt.figure(yzfig)
        plt.xlabel('z (m)')
        plt.ylabel('stdy (m)')
        
        plt.figure(TZfig)
        plt.plot()
        plt.xlabel('z (m)')
        plt.ylabel('Mean beam energy (keV)')
    if 'Aplot' in args:
        fit = optimization.curve_fit(line, inxarray, finxarray)
        Avalue=fit[0][0]
        eA=np.sqrt(np.diag(np.cov(inxarray,finxarray))[0])
        
        plt.figure(1)
        
        plt.plot(inxarray,finxarray,'r.')
        plt.xlabel('initial x (m)')
        plt.ylabel('final x (m) ')
        plt.plot(np.asarray(inxarray),Avalue*np.asarray(inxarray))

    if 'Bplot' in args:
        fit = optimization.curve_fit(line, div0array, finxarray)
        Bvalue=fit[0][0]
        eB=np.sqrt(np.diag(fit[1])[0])
        
        plt.figure(2)
       
        plt.plot(div0array,finxarray,'r.')
        plt.xlabel('initial divergence (rad)')
        plt.ylabel('final x (m) ')
        plt.plot(np.asarray(div0array),Bvalue*np.asarray(div0array))

    if 'show' in args:
        plt.show()



def EinzelLensGen(Var,Range):
    return 0

def Eplotter():
    infile='SHORTOUTSF7.TXT'

def GPTrun(ParentFolder,Mode,Temp=10,N=2000,ScreenPos=50,InitialZ=0.025,BunchSize=0.002):
    global sigmax
    if sigmax!=BunchSize:
        
        sigmax=BunchSize
    
    print N,Temp,ScreenPos,ParentFolder,InitialZ,sigmax
    
    Ulist=[]
    if Mode=='T':
        stdlist=[]
        estdlist=[]
    if Mode=='A':
        A_list=[]
        eA_list=[]
    if Mode=='B':
        B_list=[]
        eB_list=[]

    for i in range(1,21):
        sys.stdout=oldstdout
        VBackPlate=0-250*i;
        FolderName='VBack=%d' %VBackPlate
        Path="../" +Pathname +"/"+ ParentFolder+"/"+FolderName
        
        if  os.path.exists(Path):
            os.chdir(Path);
        else:
            print "Shit's fucked with ",Path
            sys.exit()
        
        Fisher()

        if Mode=='T':
            print "Running real bunch simulation in "+ FolderName
            U=[]
            stdx=[]
            for j in range(0,10):
                dynfile="TSim_beamdyn.in"
                beam=BeamGenerator(N,Temp,BunchSize,InitialZ)#N,T,sigma,central Z position
                BeamDynWriter(dynfile,"beam","Screens",beam,ScreenPos)
                        #                DynamicFile,     GroupBy,    Outtxt
                GPTCall(dynfile,"position","std1.txt",'std')
                        
                ScreenSTD,E=Plotter("std1.txt",FolderName,'RealBunchStds')
                stdx.append(ScreenSTD)#Calculate std of final positions
                U.append(E) #Add intial energy to array
                os.remove(dynfile)
            
            Ulist.append(np.mean(U))
            stdlist.append(np.mean(stdx))
            estdlist.append(np.std(stdx))
        if Mode=='A':
            print "Running A study in "+ FolderName
            dynfile="A_beamdyn.in"
            BeamDynWriter(dynfile,"A","Screens",MCPPos=ScreenPos,InitialPos=InitialZ)
        #           DynamicFile,     GroupBy,    Outtxt
            GPTCall(dynfile,"position","stdA.txt")
    
            U_i,A,eA=Plotter("stdA.txt",FolderName,i,'inx','finx','A')
            A_list.append(A)
            eA_list.append(eA)
            Ulist.append(U_i)
#            os.remove(dynfile)
        if Mode=='B':
            print "Running B study in "+ FolderName
            dynfile="B_beamdyn.in"
            BeamDynWriter(dynfile,"B","Screens",MCPPos=ScreenPos,InitialPos=InitialZ)
        #           DynamicFile,     GroupBy,    Outtxt
            GPTCall(dynfile,"position","stdB.txt")
            U_i,B,eB=Plotter("stdB.txt",FolderName,i,'inx','div0','finx','B')
            B_list.append(B)
            eB_list.append(eB)
            Ulist.append(U_i)
#            os.remove(dynfile)

        os.chdir("../../../../VaryParameters")
            
            
    if Mode=='T':
        return TFit(Ulist,stdlist,estdlist,Temp)
    if Mode=='A':
        AfileName='Avals/Adata(%s).txt' %ParentFolder
        Afile=open(AfileName,'w')
        S= Ulist , A_list, eA_list
        np.savetxt(Afile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
        Afile.close()
    if Mode=='B':
        BfileName='BVals/Bdata(%s).txt' %ParentFolder
        Bfile=open(BfileName,'w')
        S= Ulist, B_list, eB_list
        np.savetxt(Bfile,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
        Bfile.close()

#np.mean(U),np.mean(stdx),np.std(stdx)

def TFit(Ulist,stdlist,estdlist,Temp,figureN=25):
    Afit,Bfit=ABFIT()
    plt.figure(figureN)
    plt.errorbar(np.asarray(Ulist),np.asarray(stdlist),np.asarray(estdlist),fmt='.',markersize=0)
    Tfit=optimization.curve_fit(Model2, Ulist, stdlist,sigma=estdlist, maxfev=100000,p0=Temp)
    Tfitted=Tfit[0][0]
    eTfitted=np.sqrt(np.diag(Tfit[1])[0])
    print 'Temperature = %2.3f +/- %2.3f' %(Tfitted,eTfitted)
    plt.plot(Ulist,PlotModel2(Ulist,Tfitted),label='T=(%2.1f $\pm$ %2.1f)K' %(Tfitted,eTfitted))
    return Tfitted,eTfitted

def SingleGPTRun(ParentFolder,Voltage=-1000,Mode='TTraj',Temp=10,N=50,ScreenPos=50,InitialZ=0.025):
    
    FolderName='VBack=%d' %Voltage
    Path="../" +Pathname +"/"+ ParentFolder+"/"+FolderName
        
    if  os.path.exists(Path):
        os.chdir(Path);
    else:
        print "Shit's fucked with ",Path
        sys.exit()
        
    Fisher()
    if Mode=='TTraj':
        dynfile="beamdyn_Traj.in"
        beam=BeamGenerator(N,Temp,sigmax,InitialZ)
        BeamDynWriter(dynfile,"beam","Snaps",beam)
        
        print "Running real bunch simulation in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"time","std1.txt")
        
        Plotter("std1.txt",FolderName,'xz','yz','Tz','show')
        os.remove(dynfile)
    if Mode=='E':
        dynfile="beamdyn_E.in"
        beam=BeamGenerator(N,Temp,sigmax,InitialZ)
        BeamDynWriter(dynfile,"beam","Snaps",beam)
        
        print "Running real bunch simulation in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"time","std1.txt",'std')
        
        Plotter("std1.txt",FolderName,'Tz','show')
        os.remove(dynfile)

    if Mode=='ATraj':
        dynfile="beamdyn_ATraj.in"

        BeamDynWriter(dynfile,"A","Snaps",Number=N,InitialPos=InitialZ)
            
        print "Running A study in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"time","std1.txt")
        
        Plotter("std1.txt",FolderName,'xz','show')
        os.remove(dynfile)
    if Mode=='BTraj':
        dynfile="beamdyn_BTraj.in"
        
        BeamDynWriter(dynfile,"B","Snaps",Number=N,InitialPos=InitialZ)
                
        print "Running Bstudy in "+ FolderName
        #                DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"time","std1.txt")
            
        Plotter("std1.txt",FolderName,'xz','show')
        os.remove(dynfile)
    if Mode=='A':
        print "Running A study in "+ FolderName
        dynfile="A_beamdyn.in"
        BeamDynWriter(dynfile,"A","Screens",MCPPos=ScreenPos,InitialPos=InitialZ,Number=N)
        #           DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"position","stdA.txt")

        Plotter("stdA.txt",FolderName,'inx','finx','Aplot','show')
        os.remove(dynfile)
    if Mode=='B':
        print "Running A study in "+ FolderName
        dynfile="B_beamdyn.in"
        BeamDynWriter(dynfile,"B","Screens",MCPPos=ScreenPos,InitialPos=InitialZ,Number=N)
        #           DynamicFile,     GroupBy,    Outtxt
        GPTCall(dynfile,"position","stdB.txt")
        
        Plotter("stdB.txt",FolderName,'inx','div0','finx','Bplot','show')
        os.remove(dynfile)

    os.chdir("../../../../VaryParameters")

def LensGenerator(my_test,var='test'):

    for Test in my_test:
        ParentFolder='Test=%1.2f' %Test
#        ParentFolder='Tuned2'
        for i in range(1,21):
            sys.stdout=oldstdout
            VBackPlate=0-250*i;
            FolderName='VBack=%d' %VBackPlate
            print "Generating in "+ FolderName
            Path="../" +Pathname +"/"+ ParentFolder+"/"+FolderName
            print Path
            if not os.path.exists(Path):
                os.makedirs(Path);
    
            os.chdir(Path);
            print os.getcwd()
            fid4=open('Batch.bat','w');
            sys.stdout = fid4
            print ('del %0')
            fid4.close();
            os.chdir("../../../../VaryParameters")

    os.system("open '../SFBatch - Shortcut.lnk'");

def Saver(filename,*args):
    S=args
    file=open(filename,'w')
    np.savetxt(file,zip(*S),fmt='%1.3e', delimiter='      ', newline='\n',)
    file.close()

