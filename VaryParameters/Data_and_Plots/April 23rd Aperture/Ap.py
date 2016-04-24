import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

import os
print os.getcwd()

Tlist=10,20,30
A,T,eT=fc.data3('Data and Plots/April 23rd Aperture/ApertureData.txt')

plt.figure(1)

plt.errorbar(A,T,eT,fmt='.')
for T in Tlist: plt.axhline(y=T)
plt.xlabel('ApertureSize (cm)')
plt.ylabel('Fitted Temperature (K)')
plt.show()