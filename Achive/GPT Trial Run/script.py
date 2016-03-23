from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os

#beamdyn=open("beamdyn.in")

if not os.path.exists("Plots"):
    os.makedirs("Plots")

os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100\
          \n gpt_track -o result.gdf beamdyn.in boundaries.in Efac=1. GPTLICENSE=1539721513 \
          \n gdftrans -o traj1.gdf result.gdf time x y rxy z G \
          \n gdf2a -o traj1.txt traj1.gdf")

#\n fishfile -o fieldmap.gdf -g boundaries.in EinzelLens.am\

input_file='traj1.txt'

fin=open(input_file,'r')

new_part=0


thispartx=[]
thisparty=[]
thispartrxy=[]
thispartz=[]
thispartG=[]
thisparttime=[]
thispartT=[]
for l in fin:
    
    
    if(len(l.split())==0):
        
        
        plt.figure(1)
        plt.plot(thispartx, thispartz)
        
        plt.figure(2)
        plt.plot(thisparty, thispartz)
        
        
        plt.figure(3)
        plt.plot(thispartz, thispartT)
        
        plt.figure(4)
        plt.plot(thispartrxy, thispartz)
        
        plt.figure(5)
        plt.plot(thisparttime,thispartz)
        
        thispartx=[]
        thisparty=[]
        thispartrxy=[]
        thispartz=[]
        thispartG=[]
        thisparttime=[]
        thispartT=[]
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
        G=l.split()[4]
        time=l.split()[5]
        
        T=511*(float(G)-1)
        thispartx.append(x)
        thisparty.append(y)
        thispartrxy.append(rxy)
        thispartz.append(z)
        thispartG.append(G)
        thisparttime.append(time)
        thispartT.append(T)

plt.figure(1)
plt.xlabel('x (m)')
plt.ylabel('z m ')
plt.savefig("plots/x.eps")
plt.savefig("plots/x.png")

plt.figure(2)
plt.xlabel('y (m)')
plt.ylabel('z (m) ')
plt.savefig("plots/y.eps")
plt.savefig("plots/y.png")

plt.figure(3)
plt.xlabel('z (m)')
plt.ylabel('T (eV)')
plt.savefig("plots/T.eps")
plt.savefig("plots/T.png")

plt.figure(4)
plt.xlabel('rxy (m)')
plt.ylabel('z (m)')
plt.savefig("plots/rxy.eps")
plt.savefig("plots/rxy.png")

plt.figure(5)
plt.xlabel('time (s)')
plt.ylabel('z (m)')
plt.savefig("plots/z.eps")
plt.savefig("plots/z.png")





fin.close()

