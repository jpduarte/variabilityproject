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
filenameaux = 'summarydevice'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindexlin =  header.index('Iofflin') 
Ionindexlin =  header.index('Ionlin')
Vthindexlin =  header.index('Vthlin')
gmmaxindexlin =  header.index('gmaxlin')
Hfinindexlin = header.index('Hfinlin') 
SSindexlin = header.index('SSlin') 
Wtindexlin = header.index('Wtlin') 
Lgindexlin = header.index('Lglin') 
toxindexlin = header.index('toxlin') 
WbRindexlin = header.index('WbRlin') 
WbLindexlin = header.index('WbLlin') 
Nfinindexlin = header.index('Nfinlin') 

Ioffindexsat =  header.index('Ioffsat')
Ionindexsat =  header.index('Ionsat')
Vthindexsat =  header.index('Vthsat')
gmmaxindexsat =  header.index('gmaxsat')
Hfinindexsat = header.index('Hfinsat') 
SSindexsat = header.index('SSsat') 
Wtindexsat = header.index('Wtsat') 
Lgindexsat = header.index('Lgsat') 
toxindexsat = header.index('toxsat') 
WbRindexsat = header.index('WbRsat') 
WbLindexsat = header.index('WbLsat') 
Nfinindexsat = header.index('Nfinsat') 

DIBLindexsat = header.index('DIBL') 
datalist = np.loadtxt(filenameaux,skiprows = 1)

#Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindex],datalist[:,Lgindex],datalist[:,Hfinindex],datalist[:,toxindex],datalist[:,WbRindex],datalist[:,WbLindex],datalist[:,Nfinindex],eins)

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=80, facecolors='none', edgecolors='k') 
pylab.ylabel("ION,SAT (uA/um)", fontsize=18)
pylab.xlabel("Vth,SAT (V)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.savefig('IonVthVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("SS (mV/dec)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([0.1,0.25])
#pylab.ylim([60e-3,300e-3])
##pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
#pylab.savefig('IonSSVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([500,1800])
#pylab.ylim([2,10])
#pylab.savefig('IonIoffVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("DIBL (mV/V)", fontsize=18)
pylab.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = pylab.gca()
#pylab.xlim([20,100])
#pylab.savefig('DIBLIoffVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("DIBL (mV/V)", fontsize=18)
pylab.ylabel("ION,SAT (uA/um)", fontsize=18)
#pylab.xlim([20,100])
#pylab.savefig('DIBLIonVDSATcorner', dpi=300, bbox_inches='tight')
########

fignumber = 0
########################################################SAT
filenameaux = 'summarydevicecornerhspice'
factorIDSon = 1e6/(42e-3*2+7.6e-3)
factorIDSoff = 1e9/(42e-3*2+7.6e-3)

target = open(filenameaux, 'r') 
header = str.split(target.readline())

Ioffindexlin =  header.index('Iofflin') 
Ionindexlin =  header.index('Ionlin')
Vthindexlin =  header.index('Vthlin')
gmmaxindexlin =  header.index('gmaxlin')
Hfinindexlin = header.index('Hfinlin') 
SSindexlin = header.index('SSlin') 
Wtindexlin = header.index('Wtlin') 
Lgindexlin = header.index('Lglin') 
toxindexlin = header.index('toxlin') 
WbRindexlin = header.index('WbRlin') 
WbLindexlin = header.index('WbLlin') 
Nfinindexlin = header.index('Nfinlin') 

Ioffindexsat =  header.index('Ioffsat')
Ionindexsat =  header.index('Ionsat')
Vthindexsat =  header.index('Vthsat')
gmmaxindexsat =  header.index('gmaxsat')
Hfinindexsat = header.index('Hfinsat') 
SSindexsat = header.index('SSsat') 
Wtindexsat = header.index('Wtsat') 
Lgindexsat = header.index('Lgsat') 
toxindexsat = header.index('toxsat') 
WbRindexsat = header.index('WbRsat') 
WbLindexsat = header.index('WbLsat') 
Nfinindexsat = header.index('Nfinsat') 

DIBLindexsat = header.index('DIBL') 
datalist = np.loadtxt(filenameaux,skiprows = 1)

#Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindex],datalist[:,Lgindex],datalist[:,Hfinindex],datalist[:,toxindex],datalist[:,WbRindex],datalist[:,WbLindex],datalist[:,Nfinindex],eins)

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=80, facecolors='none', edgecolors='r') 
pylab.ylabel("ION,SAT (uA/um)", fontsize=18)
pylab.xlabel("Vth,SAT (V)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.savefig('IonVthVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=80, facecolors='none', edgecolors='r') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("SS (mV/dec)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([0.1,0.25])
#pylab.ylim([60e-3,300e-3])
##pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
#pylab.savefig('IonSSVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=80, facecolors='none', edgecolors='r') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
#pylab.xlim([500,1800])
#pylab.ylim([2,10])
#pylab.savefig('IonIoffVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=80, facecolors='none', edgecolors='r') 
pylab.xlabel("DIBL (mV/V)", fontsize=18)
pylab.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = pylab.gca()
#pylab.xlim([20,100])
#pylab.savefig('DIBLIoffVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=80, facecolors='none', edgecolors='r') 
pylab.xlabel("DIBL (mV/V)", fontsize=18)
pylab.ylabel("ION,SAT (uA/um)", fontsize=18)
#pylab.xlim([20,100])
#pylab.savefig('DIBLIonVDSATcorner', dpi=300, bbox_inches='tight')
########

pylab.show() 

