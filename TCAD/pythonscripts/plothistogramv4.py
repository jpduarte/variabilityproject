#plot many quantities
#Juan Duarte

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats

eins =  3.453060e-11 
#'Ioff'+' '+'Ion'+' '+'Wt'+' '+'Lg'+' '+'Hfin'+' '+'tox'+' '+'WbR'+' '+'WbL' +' '+'Nfin'+' '+'vd'+' '+'Vth'+' gmax'+' SS'+' filename'+'\n'

def UFCMparameters(Wt,Lg,Hfin,tox,WbR,WbL,Nfin,eins):
  Cins = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)*eins/tox
  Ach=Hfin*(WbR+WbL+2*Wt)/2
  Weff = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)
  return Cins,Ach,Weff
fignumber = 0

#############################################
filenameaux = 'summarydevicemchspice14nmv3'
factorIDSon = 1e3/(42e-3*2+7.6e-3)
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

Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindexsat],datalist[:,Lgindexsat],datalist[:,Hfinindexsat],datalist[:,toxindexsat],datalist[:,WbRindexsat],datalist[:,WbLindexsat],datalist[:,Nfinindexsat],eins)

fignumber+=1
plt.figure(fignumber)
modelref1 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=500, marker="s",facecolors='none', edgecolors='k')
VthSScorr_modelref =  np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])


fignumber+=1
plt.figure(fignumber)
modelref2 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="s",facecolors='none', edgecolors='k') 
IonIoffcorr_modelref =  np.corrcoef(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
modelref2 = plt.scatter(factorIDSon*datalist[:,Ionindexlin],(factorIDSoff*datalist[:,Ioffindexlin]),s=500, marker="s",facecolors='none', edgecolors='k') 
IonIoffcorr_modelref =  np.corrcoef(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
modelref3 = plt.scatter(1000*datalist[:,DIBLindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="s",facecolors='none', edgecolors='k')
DIBLIoffcorr_modelref =  np.corrcoef(datalist[:,DIBLindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
stats.probplot(1000*datalist[:,DIBLindexsat], dist="norm", plot=plt)

fignumber+=1
plt.figure(fignumber)
modelref4 = plt.scatter(1000*datalist[:,DIBLindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="o", color='r', edgecolors='k')
DIBLIoffcorr_modelref =  np.corrcoef(datalist[:,DIBLindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))


target.close()

fignumber = 0
fignumber = 0

################################################################################################################legends

##########
#############################################
filenameaux = 'summarydevicemchspice10nmv3'
factorIDSon = 1e3/(42e-3*2+7.6e-3)
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

Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindexsat],datalist[:,Lgindexsat],datalist[:,Hfinindexsat],datalist[:,toxindexsat],datalist[:,WbRindexsat],datalist[:,WbLindexsat],datalist[:,Nfinindexsat],eins)

fignumber+=1
plt.figure(fignumber)
modelred1 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=500, marker="o",facecolors='none', edgecolors='b')
VthSScorr_modelreduced =  np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])


fignumber+=1
plt.figure(fignumber)
modelred2 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="o",facecolors='none', edgecolors='b') 
IonIoffcorr_modelreduced =  np.corrcoef(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
modelred2 = plt.scatter(factorIDSon*datalist[:,Ionindexlin],(factorIDSoff*datalist[:,Ioffindexlin]),s=500, marker="o",facecolors='none', edgecolors='b') 
IonIoffcorr_modelreduced =  np.corrcoef(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
modelred3 = plt.scatter(1000*datalist[:,DIBLindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="o",facecolors='none', edgecolors='b')
DIBLIoffcorr_modelreduced =  np.corrcoef(datalist[:,DIBLindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))


fignumber+=1
plt.figure(fignumber)
stats.probplot(1000*datalist[:,DIBLindexsat], dist="norm", plot=plt)


fignumber+=1
plt.figure(fignumber)
modelred4 = plt.scatter(1000*datalist[:,DIBLindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=500, marker="o", color='b', edgecolors='k')
DIBLIoffcorr_modelreduced =  np.corrcoef(datalist[:,DIBLindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]))


target.close()





##################################
fignumber = 0

fignumber+=1
plt.figure(fignumber)
plt.ylabel("SS (mV/dec)", fontsize=21)
plt.xlabel("Vth,SAT (V)", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ modelref1,modelred1],['14 nm','10 nm'],loc=1,prop={'size':20})
#plt.legend([ modelref1,modelred1],['14 nm, p='+("%.2f" % VthSScorr_modelref [0][1]),'10 nm, p='+("%.2f" % VthSScorr_modelreduced[0][1])],loc=1,prop={'size':20})

fignumber+=1
plt.figure(fignumber)
plt.xlabel("ION,SAT (mA/um)", fontsize=21)
plt.ylabel("IOFF,SAT (nA/um)", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ modelref2,modelred2],['14 nm','10 nm'],loc=1,prop={'size':20})
#plt.legend([ modelref2,modelred2],['14 nm, p='+("%.2f" % IonIoffcorr_modelref [0][1]),'10 nm, p='+("%.2f" % IonIoffcorr_modelreduced[0][1])],loc=1,prop={'size':20})
ax.set_yscale('log')
#plt.xlim([500,1800])
plt.ylim([1,1e4])
#plt.xlim([0.8,1.7])

fignumber+=1
plt.figure(fignumber)
plt.xlabel("ION,LIN (mA/um)", fontsize=21)
plt.ylabel("IOFF,LIN (nA/um)", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ modelref2,modelred2],['14 nm','10 nm'],loc=1,prop={'size':20})
#plt.legend([ modelref2,modelred2],['14 nm, p='+("%.2f" % IonIoffcorr_modelref [0][1]),'10 nm, p='+("%.2f" % IonIoffcorr_modelreduced[0][1])],loc=1,prop={'size':20})
ax.set_yscale('log')
#plt.xlim([500,1800])
plt.ylim([1,1e4])
#plt.xlim([0.8,1.7])

fignumber+=1
plt.figure(fignumber)
plt.xlabel("DIBL (mV/V)", fontsize=21)
plt.ylabel("IOFF,SAT (nA/um)", fontsize=21)
ax = plt.gca()
ax.set_yscale('log')
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ modelref3,modelred3],['14 nm','10 nm'],loc=1,prop={'size':20})
plt.ylim([1,1e4])

fignumber+=1
plt.figure(fignumber)
plt.ylabel("DIBL (mV/V)", fontsize=21)
plt.xlabel("Quantiles", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
ax.set_title("")
plt.legend([ modelref4,modelred4],['14 nm','10 nm'],loc=2,prop={'size':20})
#plt.legend([ modelref3,modelred3],['14 nm','10 nm'],loc=1,prop={'size':20})


target.close()



plt.show() 

