#plot many quantities
#Juan Duarte

import os
import numpy as np
import pylab

fignumber = 0
########################################################SAT
filenameaux = 'montecarloresultsallfilevdsat'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindex =  header.index('Ioff')
Ionindex =  header.index('Ion')
Vthindex =  header.index('Vth')
gmmaxindex =  header.index('gmax')
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex-1,1)

datalist = np.loadtxt(filenameaux,skiprows = 1, usecols =columns)

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindex],factorIDSoff*datalist[:,Ioffindex],s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("IOFF,SAT (nA/um)", fontsize=18)
ax = pylab.gca()
ax.set_yscale('log')
pylab.xlim([500,1800])
pylab.ylim([1e1,1e4])
#pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSoff*datalist[:,Ioffindex],50)
pylab.xlim([0,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("IOFF,SAT (nA/um)", fontsize=18)
#pylab.savefig('Ioffsathistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,Vthindex],15)
pylab.xlim([0.10,0.25])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Vth,SAT (V)", fontsize=18)
pylab.savefig('vthsathistogram', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,gmmaxindex],15)
pylab.xlim([500,3500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("gmax,SAT (uA/um V)", fontsize=18)
pylab.savefig('gmsathistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,Ionindex],50)
pylab.xlim([600,1800])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
#pylab.savefig('Ionsathistogram', dpi=600, bbox_inches='tight')

target.close()
#pylab.savefig('IonIoffVDSAT', dpi=300, bbox_inches='tight')
########################################################LIN
filenameaux = 'montecarloresultsallfilevdlin'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindex =  header.index('Ioff')
Ionindex =  header.index('Ion')
Vthindex =  header.index('Vth')
gmmaxindex =  header.index('gmax')
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex-1,1)

datalist = np.loadtxt(filenameaux,skiprows = 1, usecols =columns)

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindex],factorIDSoff*datalist[:,Ioffindex],s=80, facecolors='none', edgecolors='r') 
pylab.xlabel("ION,LIN (nA/um)", fontsize=18)
pylab.ylabel("IOFF,LIN (nA/um)", fontsize=18)
ax = pylab.gca()
ax.set_yscale('log')
#pylab.xlim([500,1800])
#pylab.ylim([1e1,1e4])
#pylab.savefig('IonIoffVDLIN', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSoff*datalist[:,Ioffindex],50)
pylab.xlim([0,60])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("IOFF,LIN (nA/um)", fontsize=18)
#pylab.savefig('Iofflinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,Vthindex],15)
pylab.xlim([0.19,0.25])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Vth,LIN (V)", fontsize=18)
pylab.savefig('vthlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,gmmaxindex],15)
pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("gmax,SAT (uA/um V)", fontsize=18)
pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,Ionindex],50)
pylab.xlim([500,1000])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("ION,LIN (uA/um)", fontsize=18)
#pylab.savefig('Ionlinhistogram', dpi=600, bbox_inches='tight')

target.close()


##################################
pylab.show() 

