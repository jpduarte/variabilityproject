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
filenameaux = 'montecarloresultsallfilevdsat'
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

Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindex],datalist[:,Lgindex],datalist[:,Hfinindex],datalist[:,toxindex],datalist[:,WbRindex],datalist[:,WbLindex],datalist[:,Nfinindex],eins)

fignumber+=1
pylab.figure(fignumber)
pylab.hist(Cins,20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Cins", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(Weff,20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Weff", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(Ach,20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Ach", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,SSindex],20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("SS", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(np.log(abs(factorIDSoff*datalist[:,Ioffindex])),50)
#pylab.xlim([0,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("IOFF,SAT (nA/um)", fontsize=18)
ax = pylab.gca()
#ax.set_xscale('log')
#pylab.savefig('Ioffsathistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.scatter(datalist[:,Vthindex],datalist[:,SSindex],s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("SS", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
pylab.xlim([0.1,0.25])
pylab.ylim([60e-3,300e-3])
#pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')


fignumber+=1
pylab.figure(fignumber)
pylab.scatter(factorIDSon*datalist[:,Ionindex],np.log(factorIDSoff*datalist[:,Ioffindex]),s=80, facecolors='none', edgecolors='k') 
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
pylab.ylabel("IOFF,SAT (nA/um)", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')
pylab.xlim([500,1800])
pylab.ylim([2,10])
#pylab.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')

"""fignumber+=1
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
#pylab.savefig('vthsathistogram', dpi=300, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,gmmaxindex],15)
pylab.xlim([500,3500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("gmax,SAT (uA/um V)", fontsize=18)
#pylab.savefig('gmsathistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,Ionindex],50)
pylab.xlim([600,1800])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("ION,SAT (uA/um)", fontsize=18)
#pylab.savefig('Ionsathistogram', dpi=600, bbox_inches='tight')


fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,SSindex],20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("SS", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

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
Hfinindex = header.index('Hfin') 
SSindex = header.index('SS') 
Wtindex = header.index('Wt') 
Lgindex = header.index('Lg') 
toxindex = header.index('tox') 
WbRindex = header.index('WbR') 
WbLindex = header.index('WbL') 
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex,1)

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
#pylab.savefig('vthlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,gmmaxindex],15)
pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("gmax,SAT (uA/um V)", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(factorIDSon*datalist[:,Ionindex],50)
pylab.xlim([500,1000])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("ION,LIN (uA/um)", fontsize=18)
#pylab.savefig('Ionlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,Hfinindex],50)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("Hfin", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

fignumber+=1
pylab.figure(fignumber)
pylab.hist(datalist[:,SSindex],20)
#pylab.xlim([500,1500])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("SS", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

target.close()
"""
filenameaux = 'montecarloresultsallfilevdlin'
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
filenameindex =  header.index('filename')
columns = np.arange(0,filenameindex,1)

datalist = np.loadtxt(filenameaux,skiprows = 1, usecols =columns)


fignumber+=1
pylab.figure(fignumber)
pylab.hist(np.log(factorIDSoff*datalist[:,Ioffindex]),50)
#pylab.xlim([0,60])
pylab.ylabel("Frequency", fontsize=18)
pylab.xlabel("IOFF,LIN (nA/um)", fontsize=18)
#pylab.savefig('Iofflinhistogram', dpi=600, bbox_inches='tight')



target.close()
##################################
pylab.show() 

