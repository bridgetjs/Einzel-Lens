from __future__ import division
import matplotlib.pyplot as plt
import functions as fc
import scipy.constants as sc
import scipy.stats as st
import numpy as np
import math
c=sc.c
k=sc.k
m=sc.m_e
pi=math.pi
beam=fc.Maxwell(100000,10)# Vx,Px,Vy,Py,Vz,Pz

#print beam[0][0], beam[1][0][0],beam[2][0],beam[3][0][0],beam[4][0],beam[5][0][0],
#fc.BeamDynWriter("Test.in","beam","Screens",beam)
V=[]
P=[]
for j in range(0,35000):
    thisV=np.sqrt(beam[0][j]**2 + beam[2][j]**2 + beam[4][j]**2) #V=√(Vx^2+Vy^2+Vz^2)
    V.append(thisV)
    P.append(4*pi*thisV**2*beam[1][j][0]* beam[3][j][0]* beam[5][j][0]) # P=4π^2 V^2 Px*Py*Pz
plt.figure(1)
plt.plot(V,P,'b.')
Vels=np.linspace(0,60000)
Boltz=np.sqrt((m/(2*pi*k*10))**(3)) *4*pi*Vels**2 * np.exp((-m/(2*k*10))*Vels**2)
plt.plot(Vels,Boltz,'k')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Probability')
plt.show()