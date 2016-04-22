from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import functions as fc
import scipy.optimize as optimization
import re

def data4(filename): #function to return the data file
    file=open(filename,'r') #Open file name
    data=np.loadtxt(file) #Load data
    file.close() #close file
    x=data[:,0]; y=data[:,1]; z=data[:,2];   #read data assuming data is in the exepcted format
    return x,y,z


def solve_linear_equation ( equ ):
    match = re.match(r"(\d+)x\+(\d+)=(\d+)", equ)
    m, c, y = match.groups()
    m, c, y = float(m), float(c), float(y) # Convert from strings to numbers
    x = (y-c)/m
    print ("x = %f" % x)
TempRange=[10,20,30]

S,T,eT=data4('PlateSeptT.txt')

Sep=fc.TotSeptRange
l=len(Sep)

T10=T[0:l]
T20=T[l:2*l]
T30=T[2*l:3*l]

eT10=eT[0:l]
eT20=eT[l:2*l]
eT30=eT[2*l:3*l]

bestsept=[]

plt.figure(1)
plt.errorbar(Sep,T10,eT10,fmt='r.',markersize=0)
line = optimization.curve_fit(fc.line, Sep, T10,sigma=eT10)
m=line[0][0]
c=line[0][1]
plt.plot(np.asarray(Sep),np.asarray(Sep)*line[0][0]+line[0][1],'c',label=r'U=(%1.2f x + %1.2f )' %(line[0][0],line[0][1]))
bestsept.append( (10-c)/m)


line2 = optimization.curve_fit(fc.line, Sep, T20,sigma=eT20)
m=line2[0][0]
c=line2[0][1]
bestsept.append( (20-c)/m)

plt.plot(np.asarray(Sep),np.asarray(Sep)*line2[0][0]+line2[0][1],'m',label=r'U=(%1.2f x + %1.2f )' %(line2[0][0],line2[0][1]))
plt.errorbar(Sep,T20,eT20,fmt='g.',markersize=0)



line3 = optimization.curve_fit(fc.line, Sep, T30,sigma=eT30)
m=line3[0][0]
c=line3[0][1]
bestsept.append( (30-c)/m)
plt.plot(np.asarray(Sep),np.asarray(Sep)*line3[0][0]+line3[0][1],'y',label=r'U=(%1.2f x + %1.2f )' %(line3[0][0],line3[0][1]))
plt.errorbar(Sep,T30,eT30,fmt='b.',markersize=0)

for Temp in TempRange:
    plt.figure(1)
    plt.axhline(y=Temp,color='k')

print bestsept,np.average(bestsept)


plt.figure(1)
plt.axis([0, 3, 0, 40])
plt.xlabel('Plate Seperation (cm)')
plt.ylabel('Temperature (K)')
plt.show()
