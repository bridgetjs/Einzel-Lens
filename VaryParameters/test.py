import sys
sys.path.append('../')
import functions as fc
import matplotlib.pyplot as plt
import numpy as np

z=0.05
Temp=[10,10]
ParentFolder='Opt2'
fc.GPTrun(ParentFolder,'A',ScreenPos=50,InitialZ=z)
fc.GPTrun(ParentFolder,'B',ScreenPos=50,InitialZ=z)

AfileName='AVals/Adata(%s).txt' %ParentFolder
BfileName='BVals/Bdata(%s).txt' %ParentFolder
fc.ABSet(AfileName,BfileName)
Tlist=[]
wTlist=[]
x=[0,25]
y=[0,25]
for T in Temp:
    Tfit,eTfit=fc.GPTrun(ParentFolder,'T',ScreenPos=50,Temp=T,N=2000,InitialZ=z)
    plt.figure(1)
    plt.errorbar(T,Tfit,eTfit,fmt='.')

#    Tlist=[]
#    wTlist=[]
#    for k in range(0,3):
#
#        Tlist.append(Tfit)
#        wTlist.append(eTfit**(-2))
#
#        plt.figure(1)
#        plt.errorbar(T,Tfit,eTfit,fmt='.')
#
#    AvgT,sumofw=np.average(Tlist,weights=wTlist,returned=True)
#    eAvgT=1/np.sqrt(sumofw)
#
#    plt.figure(2)
#    plt.errorbar(T,AvgT,eAvgT,fmt='.')
#
#plt.figure(1)
#plt.plot(x,y)
#plt.xlabel('Input Temperature',fontsize=16)
#plt.ylabel('Fitted temperature (K)')
#plt.savefig('Shit/Tfits(ooops).eps')
#plt.axis([0,25,0,30])
#
#plt.figure(2)
#plt.plot(x,y)
#plt.xlabel('Input Temperature',fontsize=16)
#plt.ylabel('Fitted temperature (K)')
#plt.savefig('Tuning/Tfits(averaged).eps')
#plt.axis([0,25,0,30])


plt.figure(25)
plt.xlabel('Kinetic Energy after acceleration (keV)')
plt.ylabel('Beam size (standard deviation) (m)')
plt.legend()
plt.savefig('Shit/Curves.eps')

plt.show()
