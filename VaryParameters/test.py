import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
fc.ABSet('AVals/Adata(pos=50).txt','BVals/Bdata(pos=50).txt')
Afile=fc.Afile
print Afile
UI,Avals=fc.data(Afile)
A,B=fc.ABFIT()
plt.plot(A(UI))
plt.show()
