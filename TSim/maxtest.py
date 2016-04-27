from __future__ import division
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import functions as fc
import scipy.constants as sc
import scipy.stats as st
import numpy as np
import math
c=sc.c
k=sc.k
m=sc.m_e
pi=math.pi
T=10
beam=fc.Maxwell(10000,T)# Vx,Px,Vy,Py,Vz,Pz

#print beam[0][0], beam[1][0][0],beam[2][0],beam[3][0][0],beam[4][0],beam[5][0][0],
#fc.BeamDynWriter("Test.in","beam","Screens",beam)
V=[]
P=[]
N=min(len(beam[0]),len(beam[2]),len(beam[4]))
Vx=[v for v in beam[0][0:N]]
Vy=[v for v in beam[2][0:N]]
Vz=[v for v in beam[4][0:N]]
Px=[v for v in beam[1][0:N]]
Py=[v for v in beam[3][0:N]]
Pz=[v for v in beam[5][0:N]]

for j in range(0,N):
    thisV=np.sqrt(Vx[j]**2 + Vy[j]**2 + Vz[j]**2) #V=sqrt(Vx^2+Vy^2+Vz^2)
    V.append(thisV)
    P.append(4*pi*thisV**2*Px[j][0]* Py[j][0]* Pz[j][0]) # P=4*pi^2 V^2 Px*Py*Pz

plt.figure(5)
n, Vels,patches=plt.hist(V, bins=100, color='green',alpha=0.75,normed=True)
Boltz=np.sqrt((m/(2*pi*k*T))**(3)) *4*pi*Vels**2 * np.exp((-m/(2*k*T))*Vels**2)
plt.plot(Vels,Boltz,'k',linewidth=1.5,label='10 K Maxwell-Boltzmann distribution')
plt.xlabel('Velocity (m/s)',fontsize=16)
plt.ylabel('Probability',fontsize=16)
plt.legend(fontsize=16)
#print Boltz


beam=fc.BeamGenerator(10000,T,2e-3,0.05)
plt.figure(6)
plt.plot(beam[0],beam[2],'.')
plt.show()





plt.show()