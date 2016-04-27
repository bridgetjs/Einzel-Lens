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
ApertureSize=3.5
V2=5000;
V1=0;
LensZRange=fc.LensZRange
PlateSept=1; # Z seperation between the ring
FlightTubeLength= 50; # length of the Cavity (cm)

GrandFolder=fc.GrandFolder
Pathname="SF_Files/"+GrandFolder

if not os.path.exists(Pathname):
    os.makedirs(Pathname)
os.chdir(Pathname)

LensZ=30
ApertureSizeRange=[3.5,4]

#fc.ApertureRange
#
for LensZ in fc.LensZRange:
    ParentFolderName='Pos=%d' %LensZ

#    ParentFolderName='Aperture=%1.2f' %ApertureSize

    ##################################### Consts
    GroundRingRad=0
    CavityRadius= 10.5; # Radius of Cavity (cm)
    PlateDepths=0.2
    CavityLength= FlightTubeLength; # length of the Cavity (cm)
    OuterPlateDepth = PlateDepths; # Depth of the outer plates in z cm
    InnerPlateDepth= PlateDepths;  # Depth of the inner plate in z
    BackPlateDepth=0.1;
    GroundPlateDepth=0.2;

    GridSizeR=100
    GridSizeZ=500

    OuterElectrodeRingRad=ApertureSize; #Radius of the Outer Electrode Ring
    InnerElectrodeRingRad=ApertureSize; #Radius of the Inner Electrode Ring

    coz2=LensZ-InnerPlateDepth/2; # Middle plate is centered
    coz1=coz2 - PlateSept - OuterPlateDepth; #position first plate
    coz3=coz2 + PlateSept + InnerPlateDepth;
    PosBackPlate=0;
    PosGroundPlate=10

    InnerCurvature=0#InnerPlateDepth/2;
    OuterCurvature=0#OuterPlateDepth/2;
    OuterPlateWidth= (CavityRadius-OuterElectrodeRingRad);
    InnerPlateWidth= (CavityRadius-InnerElectrodeRingRad);
    GroundPlateWidth=CavityRadius;
    GroundCurvature=0#GroundPlateDepth/2




    for i in range(loopmin,loopmax):
        
        sys.stdout=oldstdout
        VBackPlate=0-250*i;

        #Print a string with the Voltages of the Plates

        FolderName='VBack=%d' %VBackPlate


        #Make and change to the corresponding directory
        Path=ParentFolderName+"/"+FolderName
        if not os.path.exists(ParentFolderName):
            os.makedirs(ParentFolderName)
        if not os.path.exists(Path):
            os.makedirs(Path);
        os.chdir(Path);
        #Also make directory for any plots that are produced


        #Make Einzel Lenses Mesh File
        fid=open("EinzelLens.am", 'w');

        sys.stdout = fid
        print( 'Einzel Lens');
        print( '&reg kprob=0,    ! Poisson or Pandira problem \r');
        print( 'xjfact=0.0,      ! Electrostatic problem ');
        print( 'dx=0.1,         ! Mesh interval ');
        print( 'dy=0.1,         ! Mesh interval ');
        print( 'icylin=1,        ! Cylindrical coordinates ');
        print( 'nbsup=1,         ! N boundary condition at upper edge');
        print( 'nbslo=0,         ! D boundary condition at lower edge');
        print( 'nbsrt=1,         ! N boundary condition at right edge');
        print( 'nbslf=1,         ! N boundary condition at left edge');
        print( 'ltop=10,       ! Maximum row number for field interpolation');
        print( 'conv=1 &\n');


        print('&po x=0.,y=0. & ')
        print('&po x=',CavityRadius,',y=0. &');
        print('&po x=',CavityRadius,',y=',CavityLength,' &');
        print('&po x=0.,y=',CavityLength,' &');
        print('&po x=0.,y=0. &\n\n\n');


        #Plate 1
        print('&reg mat=0,voltage=',V1,',ibound=-1 & !Plate 1');
        print('&po x=',OuterElectrodeRingRad,',y=',coz1+OuterCurvature,' & ');
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz1,' &')
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz1+OuterPlateDepth,' &');
        print('&po x=',OuterElectrodeRingRad + OuterCurvature,',y=',coz1 + OuterPlateDepth,' &');
        print('&po x=',OuterElectrodeRingRad,',y=',coz1+OuterCurvature,' & ');


        #Plate 2
        print('&reg mat=0,voltage=',V2,',ibound=-1 & !Plate 2 ');
        print('&po x=',InnerElectrodeRingRad,',y=',coz2+InnerCurvature,' & ');
        print('&po x=',InnerElectrodeRingRad+InnerPlateWidth,',y=',coz2,' &')
        print('&po x=',InnerElectrodeRingRad+InnerPlateWidth,',y=',coz2+InnerPlateDepth,' &');
        print('&po x=',InnerElectrodeRingRad + InnerCurvature,',y=',coz2 + InnerPlateDepth,' &');
        print('&po x=',InnerElectrodeRingRad,',y=',coz2+InnerCurvature,' & ');

        #Plate 3
        print('&reg mat=0,voltage=',V1,',ibound=-1 & !Plate 3 ');
        print('&po x=',OuterElectrodeRingRad,',y=',coz3+OuterCurvature,' & ');
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz3,' &')
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz3+OuterPlateDepth,' &');
        print('&po x=',OuterElectrodeRingRad + OuterCurvature,',y=',coz3 + OuterPlateDepth,' &');
        print('&po x=',OuterElectrodeRingRad,',y=',coz3+OuterCurvature,' & ');

        #Back Plate
        print('&reg mat=0,voltage=',VBackPlate,',ibound=-1 &!Accelerator Plate ');
        print('&po x=0,y=',PosBackPlate,' & ');
        print('&po x=',CavityRadius,',y=',PosBackPlate,' &');
        print('&po x=',CavityRadius,',y=',PosBackPlate+BackPlateDepth,' &');
        print('&po x=0,y=',PosBackPlate+BackPlateDepth,' &');
        print('&po x=0,y=',PosBackPlate,' & \n \n');


        #'Grounding' Plate
        print('&reg mat=0,voltage=',V1,',ibound=-1 & !Grounding Plate ');
        print('&po x=',GroundRingRad,',y=',PosGroundPlate+GroundCurvature,' & ');
        print('&po x=',GroundRingRad+GroundPlateWidth,',y=',PosGroundPlate,' &')
        print('&po x=',GroundRingRad+GroundPlateWidth,',y=',PosGroundPlate+GroundPlateDepth,' &');
        print('&po x=',GroundRingRad + GroundCurvature,',y=',PosGroundPlate + GroundPlateDepth,' &');
        print('&po x=',GroundRingRad,',y=',PosGroundPlate+GroundCurvature,' & ');
        fid.close();

        #Make interpolate input file with a Gridsize x Gridsize grid
        fid3=open('EINZELLENS.IN7','w');
        sys.stdout = fid3
        print('rect \r\n');
        print('0 0 ',CavityRadius,' ',CavityLength,' \r\n');
        print(' ',GridSizeR,' ',GridSizeZ,' \r\nend ');
        fid3.close();

        #Make a Windows Batch File to run the files through AUTOMESH, POISSON and SUPERFISH.
        fid4=open('Batch.bat','w');
        sys.stdout = fid4
        print('START /wait C:\\LANL\\automesh.exe "Z:\\Documents\\MPhys Project\\Einzel Lens\\SF_files\\'+GrandFolder+'\\'+ParentFolderName+'\\'+FolderName+'\\EinzelLens.am" \r\n');
        print('START /wait C:\\LANL\\Poisson.exe  "Z:\\Documents\\MPhys Project\\Einzel Lens\\SF_files\\'+GrandFolder+'\\'+ParentFolderName+'\\'+FolderName+'\\EINZELLENS.T35" \r\n');
        print('START /wait  C:\\LANL\\SF7.EXE "Z:\\Documents\\MPhys Project\\Einzel Lens\\SF_files\\'+GrandFolder+'\\'+ParentFolderName+'\\'+FolderName+'\\EINZELLENS.IN7" \r\n ');
        print ('del %0')
        
        fid4.close();
        os.chdir("../../")

os.chdir("../../")
sys.stdout=oldstdout

os.system("open 'SFBatch - Shortcut.lnk'");
str = raw_input("Press Enter to Continue once superfish has finished running  ");



###############################################  Prepare OUTSF7 Files for GPT    ###############################################################
for V2 in fc.VoltageRange:
    ParentFolderName='Voltage=%d' %V2
#
#for ApertureSize in fc.ApertureRange:
#    
#    ParentFolderName='Aperture=%1.2f' %ApertureSize
#for PlateSept in [1]:
#    ParentFolderName='Sept=%1.2f' %PlateSept

#for k in range(0,1):
#    ParentFolderName='HalfLength'
#
for LensZ in fc.LensZRange:
    ParentFolderName='Pos=%d' %LensZ
    
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
        os.remove(infile)
sys.stdout=oldstdout
