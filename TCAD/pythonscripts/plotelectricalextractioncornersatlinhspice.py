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

fignumber+=1
plt.figure(fignumber)
tcad1 = plt.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, facecolors='none', edgecolors='k') 
VthIoncorr_tcad =  np.corrcoef(datalist[:,Vthindexsat],factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
tcad2 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=800, facecolors='none', edgecolors='k') 
VthSScorr_tcad =  np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])

fignumber+=1
plt.figure(fignumber)
tcad3 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=800, facecolors='none', edgecolors='k') 
IoffIoncorr_tcad =  np.corrcoef(np.log(factorIDSoff*datalist[:,Ioffindexsat]),factorIDSon*datalist[:,Ionindexsat])

fignumber+=1
plt.figure(fignumber)
tcad4 = plt.scatter((1e3*datalist[:,DIBLindexsat]),(factorIDSoff*datalist[:,Ioffindexsat]),s=800, facecolors='none', edgecolors='k') 
DIBLIOFFcorr_tcad =  np.corrcoef((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]))

fignumber+=1
plt.figure(fignumber)
tcad5 = plt.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, facecolors='none', edgecolors='k') 
IonDIBLcorr_tcad =  np.corrcoef(1e3*datalist[:,DIBLindexsat],factorIDSon*datalist[:,Ionindexsat])
fignumber = 0

print np.amin(datalist[:,Ioffindexsat])
print np.amax(datalist[:,Ioffindexsat])

#########################################################################################################################################################MODEL
filenameaux = 'summarydevicecornerhspice4'
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

print Vthindexsat,Ionindexsat

#1 
fignumber+=1
plt.figure(fignumber)
model1 = plt.scatter((datalist[:,Vthindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
VthIoncorr_model =  np.corrcoef(datalist[:,Vthindexsat],factorIDSon*datalist[:,Ionindexsat])

#2 
fignumber+=1
plt.figure(fignumber)
model2 = plt.scatter(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
VthSScorr_model = np.corrcoef(datalist[:,Vthindexsat],1000*datalist[:,SSindexsat])

#3 
fignumber+=1
plt.figure(fignumber)
model3 = plt.scatter(factorIDSon*datalist[:,Ionindexsat],(factorIDSoff*datalist[:,Ioffindexsat]),s=800, marker="h",facecolors='none', edgecolors='r') 
IoffIoncorr_model =  np.corrcoef(np.log(factorIDSoff*datalist[:,Ioffindexsat]),factorIDSon*datalist[:,Ionindexsat])

#4 
fignumber+=1
plt.figure(fignumber)
model4 = plt.scatter((1e3*datalist[:,DIBLindexsat]),(factorIDSoff*datalist[:,Ioffindexsat]),s=800, marker="h",facecolors='none', edgecolors='r') 
DIBLIOFFcorr_model =  np.corrcoef((1e3*datalist[:,DIBLindexsat]),np.log(factorIDSoff*datalist[:,Ioffindexsat]))

#5
fignumber+=1
plt.figure(fignumber)
model5 = plt.scatter((1e3*datalist[:,DIBLindexsat]),factorIDSon*datalist[:,Ionindexsat],s=800, marker="h",facecolors='none', edgecolors='r') 
IonDIBLcorr_model =  np.corrcoef(1e3*datalist[:,DIBLindexsat],factorIDSon*datalist[:,Ionindexsat])
########
fignumber = 0
###################################################################################################################################################Nominal
filenameaux = 'summarydevicenominal'
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

#Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindex],datalist[:,Lgindex],datalist[:,Hfinindex],datalist[:,toxindex],datalist[:,WbRindex],datalist[:,WbLindex],datalist[:,Nfinindex],eins)

fignumber+=1
plt.figure(fignumber)
nominal1 = plt.scatter((datalist[1,Vthindexsat]),factorIDSon*datalist[1,Ionindexsat], s=800,marker="<",color='y',edgecolors='k')

fignumber+=1
plt.figure(fignumber)
nominal2 = plt.scatter(datalist[1,Vthindexsat],1000*datalist[1,SSindexsat], s=800,marker="<",color='y',edgecolors='k')

fignumber+=1
plt.figure(fignumber)
nominal3 = plt.scatter(factorIDSon*datalist[1,Ionindexsat],(factorIDSoff*datalist[1,Ioffindexsat]), s=800,marker="<",color='y',edgecolors='k') 

fignumber+=1
plt.figure(fignumber)
nominal4 = plt.scatter((1e3*datalist[1,DIBLindexsat]),(factorIDSoff*datalist[1,Ioffindexsat]), s=800,marker="<",color='y',edgecolors='k')

fignumber+=1
plt.figure(fignumber)
nominal5 = plt.scatter((1e3*datalist[1,DIBLindexsat]),factorIDSon*datalist[1,Ionindexsat], s=800,marker="<",color='y',edgecolors='k')



"""
Ioffsat Ionsat Wtsat Lgsat Hfinsat toxsat WbRsat WbLsat Nfinsat vdsat Vthsat gmaxsat SSsat Iofflin Ionlin Wtlin Lglin Hfinlin toxlin WbRlin WbLlin Nfinlin vdlin Vthlin gmaxlin SSlin DIBL

worst Ioff
9.34639289808e-08 0.000126987631237 0.00836 0.018 0.0462 0.00088 0.00418 0.00418 5.4e+18 0.86 0.124585965699 0.000221990263647 0.232313276599  5.35344988756e-09 8.16266990526e-05 0.00836 0.018 0.0462 0.00088 0.00418 0.00418 5.4e+18 0.05 0.196405664067 0.000121184854514 0.0893310836029  0.0886662942815

best Ioff
3.29770644893e-09 0.000114173713351 0.00684 0.022 0.0378 0.00072 0.00342 0.00342 6.6e+18 0.86 0.190288418916 0.000213370620407 0.0767817742939  1.26802954244e-09 6.82859841306e-05 0.00684 0.022 0.0378 0.00072 0.00342 0.00342 6.6e+18 0.05 0.214969232553 9.06403859657e-05 0.0708358899816  0.0304701402926

2.56258523377e-09 9.27272264654e-05 0.00684 0.022 0.0378 0.00088 0.00418 0.00418 6.6e+18 0.86 0.206710408488 0.00017745101957 0.0782830072375  7.84575267542e-10 5.78785203038e-05 0.00684 0.022 0.0378 0.00088 0.00418 0.00418 6.6e+18 0.05 0.237820277646 8.77044826598e-05 0.0712531833845  0.0384072458741

"""
Vthsat = 0.124585965699
SSsat = 0.232313276599
Ionsat = 0.000126987631237
Ioffsat = 9.34639289808e-08
DIBL = 0.0886662942815


fignumber = 0
markerdata = 's'
colordata = 'g'

fignumber+=1
plt.figure(fignumber)
worstIoff1 = plt.scatter(Vthsat,factorIDSon*Ionsat, s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
worstIoff2 = plt.scatter(Vthsat,1000*SSsat, s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
worstIoff3 = plt.scatter(factorIDSon*Ionsat,(factorIDSoff*Ioffsat), s=800,marker=markerdata,color=colordata,edgecolors='k') 

fignumber+=1
plt.figure(fignumber)
worstIoff4 = plt.scatter((1e3*DIBL),(factorIDSoff*Ioffsat), s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
worstIoff5 = plt.scatter((1e3*DIBL),factorIDSon*Ionsat, s=800,marker=markerdata,color=colordata,edgecolors='k')

Vthsat = 0.206710408488
SSsat = 0.0782830072375 
Ionsat = 9.27272264654e-05
Ioffsat = 2.56258523377e-09
DIBL = 0.0384072458741

fignumber = 0
markerdata = 'o'
colordata = 'c'

fignumber+=1
plt.figure(fignumber)
bestIoff1 = plt.scatter(Vthsat,factorIDSon*Ionsat, s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
bestIoff2 = plt.scatter(Vthsat,1000*SSsat, s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
bestIoff3 = plt.scatter(factorIDSon*Ionsat,(factorIDSoff*Ioffsat), s=800,marker=markerdata,color=colordata,edgecolors='k') 

fignumber+=1
plt.figure(fignumber)
bestIoff4 = plt.scatter((1e3*DIBL),(factorIDSoff*Ioffsat), s=800,marker=markerdata,color=colordata,edgecolors='k')

fignumber+=1
plt.figure(fignumber)
bestIoff5 = plt.scatter((1e3*DIBL),factorIDSon*Ionsat, s=800,marker=markerdata,color=colordata,edgecolors='k')


########
#########################################################################################################################################################Legends
fignumber = 0


fignumber+=1
plt.figure(fignumber)

plt.ylabel("ION,SAT (mA/um)", fontsize=21)
plt.xlabel("Vth,SAT (V)", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ tcad1,model1,nominal1,worstIoff1,bestIoff1],['TCAD','Model','Nominal','Highest Ioff','Lowest Ioff'],loc=1,prop={'size':20})
#plt.legend([ tcad1,model1,nominal1],['TCAD, p='+("%.2f" % VthIoncorr_tcad[0][1]),'Model, p='+("%.2f" % VthIoncorr_model[0][1]),'Nominal'],loc=1,prop={'size':20})
#IonVthVDSATcornerhspice

fignumber+=1
plt.figure(fignumber)

plt.xlabel("Vth,SAT (V)", fontsize=21)
plt.ylabel("SS (mV/dec)", fontsize=21)
ax = plt.gca()
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ tcad1,model1,nominal1,worstIoff1,bestIoff1],['TCAD','Model','Nominal','Highest Ioff','Lowest Ioff'],loc=1,prop={'size':20})
#plt.legend([ tcad2,model2,nominal2],['TCAD, p='+("%.2f" % VthSScorr_tcad[0][1]),'Model, p='+("%.2f" % VthSScorr_model[0][1]),'Nominal'],loc=1,prop={'size':20})
#plt.savefig('VthSSVDSATcornerhspice7', dpi=600, bbox_inches='tight')


fignumber+=1
plt.figure(fignumber)

plt.xlabel("ION,SAT (mA/um)", fontsize=21)
plt.ylabel("IOFF,SAT (nA/um)", fontsize=21)
ax = plt.gca()
ax.set_yscale('log')
#plt.xlim([500,1800])
#plt.ylim([2,10])
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ tcad1,model1,nominal1,worstIoff1,bestIoff1],['TCAD','Model','Nominal','Highest Ioff','Lowest Ioff'],loc=2,prop={'size':20})
#plt.legend([ tcad3,model3,nominal3],['TCAD, p='+("%.2f" % IoffIoncorr_tcad[0][1]),'Model, p='+("%.2f" % IoffIoncorr_model[0][1]),'Nominal'],loc=2,prop={'size':20})
#plt.savefig('IonIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')


fignumber+=1
plt.figure(fignumber)

plt.xlabel("DIBL (mV/V)", fontsize=21)
plt.ylabel("IOFF,SAT (nA/um)", fontsize=21)
ax = plt.gca()
ax.set_yscale('log')
plt.xlim([10,100])
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ tcad1,model1,nominal1,worstIoff1,bestIoff1],['TCAD','Model','Nominal','Highest Ioff','Lowest Ioff'],loc=2,prop={'size':20})
#plt.legend([ tcad4,model4,nominal4],['TCAD, p='+("%.2f" % DIBLIOFFcorr_tcad[0][1]) ,'Model, p='+("%.2f" % DIBLIOFFcorr_model[0][1]) ,'Nominal'],loc=2,prop={'size':20})
#plt.savefig('DIBLIoffVDSATcornerhspice7', dpi=600, bbox_inches='tight')

fignumber+=1
plt.figure(fignumber)

plt.xlabel("DIBL (mV/V)", fontsize=21)
plt.ylabel("ION,SAT (mA/um)", fontsize=21)
plt.xlim([10,100])
plt.xticks(fontsize=21) 
plt.yticks(fontsize=21) 
plt.legend([ tcad1,model1,nominal1,worstIoff1,bestIoff1],['TCAD','Model','Nominal','Highest Ioff','Lowest Ioff'],loc=2,prop={'size':20})
#plt.legend([ tcad5,model5,nominal5],['TCAD, p='+("%.2f" % IonDIBLcorr_tcad[0][1]),'Model, p='+("%.2f" % IonDIBLcorr_model[0][1]),'Nominal'],loc=2,prop={'size':20})
#plt.savefig('DIBLIonVDSATcornerhspice8', dpi=600)

plt.show() 

