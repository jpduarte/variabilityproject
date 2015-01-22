#plot many quantities
#Juan Duarte

import os
import numpy as np
import matplotlib.pyplot as plt

eins =  3.453060e-13 
#'Ioff'+' '+'Ion'+' '+'Wt'+' '+'Lg'+' '+'Hfin'+' '+'tox'+' '+'WbR'+' '+'WbL' +' '+'Nfin'+' '+'vd'+' '+'Vth'+' gmax'+' SS'+' filename'+'\n'

def UFCMparameters(Wt,Lg,Hfin,tox,WbR,WbL,Nfin,eins):
  Cins = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)*eins/tox
  Ach=Hfin*(WbR+WbL+2*Wt)/2
  Weff = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)
  return Cins,Ach,Weff
fignumber = 0
########################################################SAT
filenameaux = 'summarydevicecorner'
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
plt.figure(fignumber)
tcad1 = plt.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, facecolors='none', edgecolors='k') 
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.savefig('IonVthVDSATcorner', dpi=300, bbox_inches='tight')
VthIoncorr_tcad =  np.corrcoef(datalist[:,Vthindexsat],factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
tcad2 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=800, facecolors='none', edgecolors='k') 
plt.xlabel("Vth,SAT (V)", fontsize=18)
plt.ylabel("SS (mV/dec)", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([0.1,0.25])
#plt.ylim([60e-3,300e-3])
##plt.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
#plt.savefig('IonSSVDSATcorner', dpi=300, bbox_inches='tight')
print 'Correlation Vth-SS'
VthSScorr_tcad =  np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])

fignumber+=1
plt.figure(fignumber)
tcad3 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=800, facecolors='none', edgecolors='k') 
plt.xlabel("ION,SAT (uA/um)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([500,1800])
#plt.ylim([2,10])
#plt.savefig('IonIoffVDSATcorner', dpi=300, bbox_inches='tight')
IoffIoncorr_tcad =  np.corrcoef(np.log(factorIDSoff*datalist[:,Ioffindexsat]),factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
tcad4 = plt.scatter((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=800, facecolors='none', edgecolors='k') 
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#plt.xlim([20,100])
#plt.savefig('DIBLIoffVDSATcorner', dpi=300, bbox_inches='tight')

print 'Correlation DIBL-Ioff'
DIBLIOFFcorr_tcad =  np.corrcoef((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
tcad5 = plt.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, facecolors='none', edgecolors='k') 
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
#plt.xlim([20,100])
#plt.savefig('DIBLIonVDSATcorner', dpi=300, bbox_inches='tight')
########
IonDIBLcorr_tcad =  np.corrcoef(1e3*datalist[:,DIBLindexsat],factorIDSon*datalist[:,Ionindexsat])
fignumber = 0


#########################################################################################################################################################MODEL
filenameaux = 'summarydevicecornerhspice4'
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
print Vthindexsat,Ionindexsat
fignumber+=1
plt.figure(fignumber)
model1 = plt.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
ax = plt.gca()
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.legend([ tcad1,model1],['TCAD','Model'],loc=1,prop={'size':20})
#ax.set_yscale('log')
#plt.savefig('IonVthVDSATcornerhspice7', dpi=600, bbox_inches='tight')
VthIoncorr_model =  np.corrcoef(datalist[:,Vthindexsat],factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
model2 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
plt.xlabel("Vth,SAT (V)", fontsize=18)
plt.ylabel("SS (mV/dec)", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([0.1,0.25])
#plt.ylim([60e-3,300e-3])
##plt.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.legend([ tcad2,model2],['TCAD','Model'],loc=1,prop={'size':20})
#plt.savefig('VthSSVDSATcornerhspice7', dpi=600, bbox_inches='tight')

print 'Correlation Vth-SS'
VthSScorr_model = np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])

fignumber+=1
plt.figure(fignumber)
model3 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=800, marker="h",facecolors='none', edgecolors='r') 
plt.xlabel("ION,SAT (uA/um)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([500,1800])
#plt.ylim([2,10])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.legend([ tcad3,model3],['TCAD','Model'],loc=2,prop={'size':20})
#plt.savefig('IonIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')
IoffIoncorr_model =  np.corrcoef(np.log(factorIDSoff*datalist[:,Ioffindexsat]),factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
model4 = plt.scatter((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]),s=800, marker="h",facecolors='none', edgecolors='r') 
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
plt.xlim([10,100])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.legend([ tcad4,model4],['TCAD','Model'],loc=2,prop={'size':20})
#plt.savefig('DIBLIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')
print 'Correlation DIBL-Ioff'
DIBLIOFFcorr_model =  np.corrcoef((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]))


fignumber+=1
plt.figure(fignumber)
model5 = plt.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlim([10,100])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.legend([ tcad5,model5],['TCAD','Model'],loc=2,prop={'size':20})
#plt.savefig('DIBLIonVDSATcornerhspice7', dpi=600)
IonDIBLcorr_model =  np.corrcoef(1e3*datalist[:,DIBLindexsat],factorIDSon*datalist[:,Ionindexsat])
########
fignumber = 0
###################################################################################################################################################Nominal
filenameaux = 'summarydevicenominal'
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
plt.figure(fignumber)
nominal1 = plt.scatter((datalist[1,Vthindexsat]),factorIDSon*datalist[1,Ionindexsat], s=800,marker="<",color='y',edgecolors='k')
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.savefig('IonVthVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
plt.figure(fignumber)
nominal2 = plt.scatter(datalist[1,Vthindexsat],1000*datalist[1,SSindexsat], s=800,marker="<",color='y',edgecolors='k')
plt.xlabel("Vth,SAT (V)", fontsize=18)
plt.ylabel("SS (mV/dec)", fontsize=18)
ax = plt.gca()



#ax.set_yscale('log')
#plt.xlim([0.1,0.25])
#plt.ylim([60e-3,300e-3])
##plt.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
#plt.savefig('IonSSVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
plt.figure(fignumber)
nominal3 = plt.scatter(factorIDSon*datalist[1,Ionindexsat],np.log(factorIDSoff*datalist[1,Ioffindexsat]), s=800,marker="<",color='y',edgecolors='k') 
plt.xlabel("ION,SAT (uA/um)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([500,1800])
#plt.ylim([2,10])
#plt.savefig('IonIoffVDSATcorner', dpi=300, bbox_inches='tight')

print DIBLindexsat,Ioffindexsat

fignumber+=1
plt.figure(fignumber)
nominal4 = plt.scatter((1e3*datalist[1,DIBLindexsat]),np.log(factorIDSoff*datalist[1,Ioffindexsat]), s=800,marker="<",color='y',edgecolors='k')
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#plt.xlim([20,100])
#plt.savefig('DIBLIoffVDSATcorner', dpi=300, bbox_inches='tight')

fignumber+=1
plt.figure(fignumber)
nominal5 = plt.scatter((1e3*datalist[1,DIBLindexsat]),factorIDSon*datalist[1,Ionindexsat], s=800,marker="<",color='y',edgecolors='k')
plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
#plt.xlim([20,100])
#plt.savefig('DIBLIonVDSATcorner', dpi=300, bbox_inches='tight')
########
#########################################################################################################################################################Legends
fignumber = 0


fignumber+=1
plt.figure(fignumber)

plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
ax = plt.gca()
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
plt.legend([ tcad1,model1,nominal1],['TCAD, p='+("%.2f" % VthIoncorr_tcad[0][1]),'Model, p='+("%.2f" % VthIoncorr_model[0][1]),'Nominal'],loc=1,prop={'size':20})
#ax.set_yscale('log')
#plt.savefig('IonVthVDSATcornerhspice7', dpi=600, bbox_inches='tight')


fignumber+=1
plt.figure(fignumber)

plt.xlabel("Vth,SAT (V)", fontsize=18)
plt.ylabel("SS (mV/dec)", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([0.1,0.25])
#plt.ylim([60e-3,300e-3])
##plt.savefig('IonIoffVDSAT', dpi=600, bbox_inches='tight')
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
plt.legend([ tcad2,model2,nominal2],['TCAD, p='+("%.2f" % VthSScorr_tcad[0][1]),'Model, p='+("%.2f" % VthSScorr_model[0][1]),'Nominal'],loc=1,prop={'size':20})
#plt.savefig('VthSSVDSATcornerhspice7', dpi=600, bbox_inches='tight')


fignumber+=1
plt.figure(fignumber)

plt.xlabel("ION,SAT (uA/um)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
#ax.set_yscale('log')
#plt.xlim([500,1800])
#plt.ylim([2,10])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
plt.legend([ tcad3,model3,nominal3],['TCAD, p='+("%.2f" % IoffIoncorr_tcad[0][1]),'Model, p='+("%.2f" % IoffIoncorr_model[0][1]),'Nominal'],loc=2,prop={'size':20})
#plt.savefig('IonIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')


fignumber+=1
plt.figure(fignumber)

plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("LOG(IOFF,SAT (nA/um))", fontsize=18)
ax = plt.gca()
plt.xlim([10,100])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
plt.legend([ tcad4,model4,nominal4],['TCAD, p='+("%.2f" % DIBLIOFFcorr_tcad[0][1]) ,'Model, p='+("%.2f" % DIBLIOFFcorr_model[0][1]) ,'Nominal'],loc=2,prop={'size':20})
#plt.savefig('DIBLIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')

fignumber+=1
plt.figure(fignumber)

plt.xlabel("DIBL (mV/V)", fontsize=18)
plt.ylabel("ION,SAT (uA/um)", fontsize=18)
plt.xlim([10,100])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
plt.legend([ tcad5,model5,nominal5],['TCAD, p='+("%.2f" % IonDIBLcorr_tcad[0][1]),'Model, p='+("%.2f" % IonDIBLcorr_model[0][1]),'Nominal'],loc=2,prop={'size':20})
#plt.savefig('DIBLIonVDSATcornerhspice8', dpi=600)

plt.show() 

