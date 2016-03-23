from __future__ import print_function
from matplotlib import cm
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import os



OuterElectrodeRingRad=1.5; #Radius of the Outer Electrode Ring
InnerElectrodeRingRad=1.5; #Radius of the Inner Electrode Ring

CavityLength= 16; # length of the Cavity (cm)
CavityRadius= 8; # Radius of Cavity (cm)
PlateSept=2; # Z seperation between the ring
OuterPlateDepth = 0.25; # Depth of the outer plates in z cm
InnerPlateDepth= 0.25;  # Depth of the inner plate in z
BackPlateDepth=0.25;
V1=100;
V2=5000;
PosBackPlate=0;

OuterCurvature=OuterPlateDepth/2;

GridSize=501;

coz2=8-InnerPlateDepth/2; # Middle plate is centered
coz1=coz2-PlateSept - OuterPlateDepth; #position first plate
coz3=coz2 + PlateSept + InnerPlateDepth;
coz2=8-InnerPlateDepth/2; # Middle plate is centered
coz1=coz2-PlateSept - OuterPlateDepth; #position first plate
coz3=coz2 + PlateSept + InnerPlateDepth;

InnerCurvature=InnerPlateDepth/2;
OuterPlateWidth= (CavityRadius-OuterElectrodeRingRad-1);
InnerPlateWidth= (CavityRadius-InnerElectrodeRingRad-1);
ParentFolderName="AcceleratorPlates"


for i in range(0,5):
        VBackPlate=0-100*i;
        
        
        #Print a string with the Voltages of the Plates
        
        FolderName='VBack=%d' %VBackPlate
        
        
        #Make and change to the corresponding directory
        Path=ParentFolderName+ "/"+ FolderName
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
        print( 'xjfact=0.0,      ! Electrostatic problem \n');
        print( 'dx=0.05,         ! Mesh interval \n');
        print( 'dy=0.05,         ! Mesh interval \n');
        print( 'icylin=1,        ! Cylindrical coordinates \n');
        print( 'nbsup=0,         ! D boundary condition at upper edge\n');
        print( 'nbslo=0,         ! D boundary condition at lower edge\n');
        print( 'nbsrt=0,         ! D boundary condition at right edge\n');
        print( 'nbslf=1,         ! N boundary condition at left edge\n');
        print( 'ltop=10,       ! Maximum row number for field interpolation\n');
        print( 'conv=1 &\n\n');
        
        
        print('&po x=0.,y=0. & \n')
        print('&po x=',CavityRadius,',y=0. &\n');
        print('&po x=',CavityRadius,',y=',CavityLength,' &\n');
        print('&po x=0.,y=',CavityLength,' &\n');
        print('&po x=0.,y=0. &\n\n\n');
        
        
        #Plate 1
        print('&reg mat=0,voltage=',V1,',ibound=-1 & !Plate 1\n');
        print('&po x=',OuterElectrodeRingRad,',y=',coz1+OuterCurvature,' & \n');
        print('&po nt=2,x0=',OuterElectrodeRingRad + OuterCurvature,',y0=',coz1+OuterCurvature,',r=',OuterCurvature,',theta=270.& \n');
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz1,' &\n')
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz1+OuterPlateDepth,' &\n');
        print('&po x=',OuterElectrodeRingRad + OuterCurvature,',y=',coz1 + OuterPlateDepth,' &\n');
        print('&po nt=2,x0=',OuterElectrodeRingRad + OuterCurvature,',y0=',coz1+OuterPlateDepth-OuterCurvature,',r=',OuterCurvature,',theta=180. &\n\n\n');
        
        #Plate 2
        print('&reg mat=0,voltage=',V2,',ibound=-1 & !Plate 2 \n');
        
        
        print('&po x=',InnerElectrodeRingRad,',y=',coz2+InnerCurvature,' & \n');
        print('&po nt=2,x0=',InnerElectrodeRingRad + InnerCurvature,',y0=',coz2+InnerCurvature,',r=',InnerCurvature,',theta=270.& \n');
        print('&po x=',InnerElectrodeRingRad+InnerPlateWidth,',y=',coz2,' &\n')
        print('&po x=',InnerElectrodeRingRad+InnerPlateWidth,',y=',coz2+InnerPlateDepth,' &\n');
        print('&po x=',InnerElectrodeRingRad + InnerCurvature,',y=',coz2 + InnerPlateDepth,' &\n');
        print('&po nt=2,x0=',InnerElectrodeRingRad + InnerCurvature,',y0=',coz2+InnerPlateDepth-InnerCurvature,',r=',InnerCurvature,',theta=180. &\n\n\n');
        
        
        #Plate 3
        print('&reg mat=0,voltage=',V1,',ibound=-1 & !Plate 3 \n');
        
        print('&po x=',OuterElectrodeRingRad,',y=',coz3+OuterCurvature,' & \n');
        print('&po nt=2,x0=',OuterElectrodeRingRad + OuterCurvature,',y0=',coz3+OuterCurvature,',r=',OuterCurvature,',theta=270.& \n');
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz3,' &\n')
        print('&po x=',OuterElectrodeRingRad+OuterPlateWidth,',y=',coz3+OuterPlateDepth,' &\n');
        print('&po x=',OuterElectrodeRingRad + OuterCurvature,',y=',coz3 + OuterPlateDepth,' &\n');
        print('&po nt=2,x0=',OuterElectrodeRingRad + OuterCurvature,',y0=',coz3+OuterPlateDepth-OuterCurvature,',r=',OuterCurvature,',theta=180. &\n\n\n')
        
        
        #Back Plate
        print('&reg mat=0,voltage=',VBackPlate,',ibound=-1 &!Accelerator Plate \n');
        print('&po x=0,y=',PosBackPlate,' & \n');
        print('&po x=',CavityRadius,',y=',PosBackPlate,' &\n');
        print('&po x=',CavityRadius,',y=',PosBackPlate+BackPlateDepth,' &\n');
        print('&po x=0,y=',PosBackPlate+BackPlateDepth,' &\n');
        print('&po x=0,y=',PosBackPlate,' &\n \n \n');
        fid.close();
        
        
        #Make interpolate input file with a Gridsize x Gridsize grid
        fid3=open('EINZELLENS.IN7','w');
        sys.stdout = fid3
        print('grid \n');
        print('0 0 ',CavityRadius,' ',CavityLength,' \n');
        print(' ',GridSize-1,' ',GridSize-1,' \nend \n ');
        fid3.close();
        
        #Make a Windows Batch File to run the files through AUTOMESH, POISSON and SUPERFISH.
        fid4=open('Batch.bat','w');
        sys.stdout = fid4
        print('START /wait C:\\LANL\\automesh.exe "Z:\\Documents\\MPhys Project\\Einzel Lens\\'+ParentFolderName+'\\'+FolderName+'\\EinzelLens.am" \r\n');
        print('START /wait C:\\LANL\\Poisson.exe  "Z:\\Documents\\MPhys Project\\Einzel Lens\\'+ParentFolderName+'\\'+FolderName+'\\EINZELLENS.T35" \r\n');
        print('START "" C:\\LANL\\SF7.EXE "Z:\\Documents\\MPhys Project\\Einzel Lens\\'+ParentFolderName+'\\'+FolderName+'\\EINZELLENS.IN7" \r\n ');
        fid4.close()

        os.chdir("../../")

#os.system("open 'py - Shortcut.lnk'");

str = raw_input("Enter your input: ");
print ('Received input is : ', str)

#for i in range(0,5):
#
#    VBackPlate=0-100*i;
#
#
#    #Print a string with the Voltages of the Plates
#
#    FolderName='VBack=%d' %VBackPlate
#    Path=ParentFolderName+ "/"+ FolderName
#
#    #Make and change to the corresponding directory
#    infile=Path + "OUTSF7.TXT"
#    outfile=Path + "SHORTOUTSF7.TXT"
#    numline=33 #3 lines to skip
#    p=""
#    o=open(infile,"a")
#    f=open(outfile)
#    for i in range(numline):
#        f.next()
#    for line in f:
#        if p:
#            o.write(p)
#        p=line
#    f.close()
#    o.close()
