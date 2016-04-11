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
T=5
beam=fc.Maxwell(100000,T)# Vx,Px,Vy,Py,Vz,Pz

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

#for j in range(0,350000):
#    thisV=np.sqrt(beam[0][j]**2 + beam[2][j]**2 + beam[4][j]**2) #V=sqrt(Vx^2+Vy^2+Vz^2)
#    V.append(thisV)
#    P.append(4*pi*thisV**2*beam[1][j][0]* beam[3][j][0]* beam[5][j][0]) # P=4*pi^2 V^2 Px*Py*Pz
plt.figure(1)
plt.plot(V,P,'b.')
Vels=np.linspace(0,max(V),10000)
Boltz=np.sqrt((m/(2*pi*k*T))**(3)) *4*pi*Vels**2 * np.exp((-m/(2*k*T))*Vels**2)
plt.plot(Vels,Boltz,'k')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Probability')
print np.sum(P)
plt.show()