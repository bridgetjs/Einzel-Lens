from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../')
import scipy.optimize as optimization
import scipy.interpolate as si
import scipy.constants as sc
import functions as fc
import timeit
oldstdout = sys.stdout

N=2000

AfileName='AVals/Adata(Tuned2).txt'
BfileName='BVals/Bdata(Tuned2).txt'

fc.ABSet(AfileName,BfileName)
print fc.Afile,fc.Bfile
Afit,Bfit=fc.ABFIT()

Ulist=[]
stdlist=[]
estdlist=[]

Pathname=fc.Pathname
ParentFolder='Tuned2'

print fc.STDrun(10,'Tuned2',2000,100,'T')
