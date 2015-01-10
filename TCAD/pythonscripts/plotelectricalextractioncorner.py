#plot many quantities
#Juan Duarte

import os
import numpy as np
import pylab

eins =  3.453060e-13 
#'Ioff'+' '+'Ion'+' '+'Wt'+' '+'Lg'+' '+'Hfin'+' '+'tox'+' '+'WbR'+' '+'WbL' +' '+'Nfin'+' '+'vd'+' '+'Vth'+' gmax'+' SS'+' filename'+'\n'

def UFCMparameters(Wt,Lg,Hfin,tox,WbR,WbL,Nfin,eins):
  Cins = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)*eins/tox
  Ach=Hfin*(WbR+WbL+2*Wt)/2
  Weff = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)
  return Cins,Ach,Weff

fignumber = 0
########################################################SAT
filenameaux = 'montecarloresultsallfilevdsatcorner'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindex =  header.index('Ioff')
Ionindex =  header.index('Ion')
Vthindex =  header.index('Vth')
gmmaxindex =  header.index('gmax')
Hfinindex = header.index('Hfin') 
SSindex = header.index('SS') 
Wtindex = header.index('Wt') 
Lgindex = header.index('Lg') 
toxindex = header.index('tox') 
WbRindex = header.index('WbR') 
WbLindex = header.index('WbL') 
Nfinindex = header.index('Nfin') 
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex,1)

datalist = np.loadtxt(filenameaux,skiprows = 1, usecols =columns)

#Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindex],datalist[:,Lgindex],datalist[:,Hfinindex],datalist[:,toxindex],datalist[:,WbRindex],datalist[:,WbLindex],datalist[:,Nfinindex],eins)

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((datalist[:,Vthindex]),factorIDSon*datalist[:,Ionindex],s=80, facecolors='none', edgecolors='r') 
pylab.ylabel("ION,SAT (uA/um)", fontsize=18)
pylab.xlabel("Vth,SAT (V)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(datalist[:,Vthindex],1000*datalist[:,SSindex],s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("SS (mV/dec)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([0.1,0.25])
#pylab.ylim([60e-3,300e-3])
#pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')


fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindex],np.log(factorIDSoff*datalist[:,Ioffindex]),s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([500,1800])
#pylab.ylim([2,10])
#pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')


##################################
########################################################SAT
filenameaux = 'montecarloresultsallfilevdlincorner'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindex =  header.index('Ioff')
Ionindex =  header.index('Ion')
Vthindex =  header.index('Vth')
gmmaxindex =  header.index('gmax')
Hfinindex = header.index('Hfin') 
SSindex = header.index('SS') 
Wtindex = header.index('Wt') 
Lgindex = header.index('Lg') 
toxindex = header.index('tox') 
WbRindex = header.index('WbR') 
WbLindex = header.index('WbL') 
Nfinindex = header.index('Nfin') 
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex,1)

datalist = np.loadtxt(filenameaux,skiprows = 1, usecols =columns)


fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindex],(datalist[:,Vthindex]),s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,LIN (uA/um)", fontsize=18)
pylab.ylabel("Vth,LIN (V)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((datalist[:,Vthindex]),factorIDSon*datalist[:,Ionindex],s=80, facecolors='none', edgecolors='k') 
pylab.ylabel("ION,LIN (uA/um)", fontsize=18)
pylab.xlabel("Vth,LIN (V)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')

pylab.show() 

