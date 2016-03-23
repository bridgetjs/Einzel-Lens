from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.optimize as optimization
import sys 
#beamdyn=open("beamdyn.in")

def func(x, a, b):
    return a*x + b


if not os.path.exists("Plots"):
    os.makedirs("Plots")

os.system("ASCI2GDF -o fieldmap.gdf SHORTOUTSF7.TXT R 1e-2 Z 1e-2 Er 100 Ez 100 absE 100 V 100\
          \n gpt_track -o result.gdf beamdyn.in boundaries.in Efac=1. GPTLICENSE=1539721513 \
          \n gdftrans -o std1.gdf result.gdf position x y rxy z G \
          \n gdf2a -o std1.txt std1.gdf")
#gdftrans -v  -o std.gdf result.gdf ID x y z rxy
#\n fishfile -o fieldmap.gdf -g boundaries.in EinzelLens.am\

input_file='std1.txt'

fin=open(input_file,'r')

new_part=0


thispartx=[]
thisparty=[]
thispartrxy=[]
thispartz=[]
thispartG=[]
thisparttime=[]
thispartT=[]
inxarray=[]
finxarray=[]
for l in fin:
    
    
    if(len(l.split())==0):
        inx=float(thispartx[0])
        finx=float(thispartx[1])
        #print inx , finx
        inxarray.append(inx)
        finxarray.append(finx)
        #print(float(thispartx[1])/float(thispartx[0]))
        
        plt.figure(2)
        plt.plot(inx   ,finx,'r.')
#        
#        plt.figure(2)
#        plt.plot(thisparty, thispartz)
#        
#        
#        plt.figure(3)
#        plt.plot(thispartz, thispartT)
#        
#        plt.figure(4)
#        plt.plot(thispartrxy, thispartz)
#        
#        plt.figure(5)
#        plt.plot(thisparttime,thispartz)

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
plt.plot(inxarray,finxarray,'r.')
plt.xlabel('initial x (m)')
plt.ylabel('final x (m) ')
a0=([0.0, 0.0])
b= optimization.curve_fit(func, inxarray, finxarray, a0)
print "The matrix element is",b[0][0]
plt.plot(np.asarray(inxarray),b[0][0]*np.asarray(inxarray),'k')
plt.savefig("plots/A(Vacc=-2.5kV).eps")
plt.savefig("plots/A(Vacc=-2.5kV).png")

#plt.show()
#
#plt.figure(1)
#plt.xlabel('x (m)')
#plt.ylabel('z m ')
#plt.savefig("plots/x.eps")
#plt.savefig("plots/x.png")
#
#plt.figure(2)
#plt.xlabel('y (m)')
#plt.ylabel('z (m) ')
#plt.savefig("plots/y.eps")
#plt.savefig("plots/y.png")
#
#plt.figure(3)
#plt.xlabel('z (m)')
#plt.ylabel('T (eV)')
#plt.savefig("plots/T.eps")
#plt.savefig("plots/T.png")
#
#plt.figure(4)
#plt.xlabel('rxy (m)')
#plt.ylabel('z (m)')
#plt.savefig("plots/rxy.eps")
#plt.savefig("plots/rxy.png")
#
#plt.figure(5)
#plt.xlabel('time (s)')
#plt.ylabel('z (m)')
#plt.savefig("plots/z.eps")
#plt.savefig("plots/z.png")
#




fin.close()

