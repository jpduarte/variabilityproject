#plot many quantities
#Juan Duarte

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats
eins =  3.453060e-13 
#'Ioff'+' '+'Ion'+' '+'Wt'+' '+'Lg'+' '+'Hfin'+' '+'tox'+' '+'WbR'+' '+'WbL' +' '+'Nfin'+' '+'vd'+' '+'Vth'+' gmax'+' SS'+' filename'+'\n'

def UFCMparameters(Wt,Lg,Hfin,tox,WbR,WbL,Nfin,eins):
  Cins = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)*eins/tox
  Ach=Hfin*(WbR+WbL+2*Wt)/2
  Weff = (np.sqrt(Hfin**2+WbR**2)+np.sqrt(Hfin**2+WbL**2)+Wt)
  return Cins,Ach,Weff
fignumber = 0
########################################################SAT
filenameaux = 'summarydevicetcadmc'
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
plt.hist(datalist[:,Ioffindexsat],50,alpha=0.9,label='TCAD')
#pylab.xlim([500,1500])
plt.ylabel("Frequency", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
#pylab.savefig('gmlinhistogram', dpi=600, bbox_inches='tight')

target.close()
fignumber = 0
#############################################
filenameaux = 'summarydevicemchspice'
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

fignumber+=1
plt.figure(fignumber)
plt.hist(datalist[:,Ioffindexsat],50,alpha=0.7,label='Model')
#pylab.xlim([500,1500])
plt.ylabel("Frequency", fontsize=18)
plt.xlabel("Vth,SAT (V)", fontsize=18)
plt.legend(loc='upper right')



Cins,Ach,Weff = UFCMparameters(datalist[:,Wtindexsat],datalist[:,Lgindexsat],datalist[:,Hfinindexsat],datalist[:,toxindexsat],datalist[:,WbRindexsat],datalist[:,WbLindexsat],datalist[:,Nfinindexsat],eins)


fignumber+=1
plt.figure(fignumber)
model3 = plt.scatter(Weff,Cins,s=300, marker="o",facecolors='none', edgecolors='k') 
ax = plt.gca()
plt.ylabel("Cins (F/m)", fontsize=18)
plt.xlabel("Weff (m)", fontsize=18)

fignumber+=1
plt.figure(fignumber)
model3 = plt.scatter(Weff,Ach,s=300, marker="o",facecolors='none', edgecolors='k') 
ax = plt.gca()
plt.ylabel("Ach (m^2)", fontsize=18)
plt.xlabel("Weff (m)", fontsize=18)

print 'Cins'
Cinssigma = np.std(Cins)
Cinsmean = np.mean(Cins)
print Cinssigma
print Cinsmean

print 'Weff'
Weffsigma = np.std(Weff)
Weffmean = np.mean(Weff)
print Weffsigma
print Weffmean

print 'Ach'
Achsigma = np.std(Ach)
Achmean = np.mean(Ach)
print Achsigma
print Achmean


fignumber+=1
plt.figure(fignumber)
plt.hist(Cins,50,alpha=0.9,label='Model')
mu, sigma = Cinsmean, Cinssigma # mean and standard deviation
Cinsnew = np.random.normal(mu, sigma, 2000)
plt.hist(Cinsnew,50,alpha=0.7,label='Model')
plt.ylabel("Frequency", fontsize=18)
plt.xlabel("Cins (F/m)", fontsize=18)

fignumber+=1
plt.figure(fignumber)
plt.hist(Ach,50,alpha=0.9,label='Model')
mu, sigma = Achmean, Achsigma # mean and standard deviation
Achnew = np.random.normal(mu, sigma, 2000)
plt.hist(Achnew,50,alpha=0.7,label='Model')
plt.ylabel("Frequency", fontsize=18)
plt.xlabel("Ach (m^2)", fontsize=18)

fignumber+=1
plt.figure(fignumber)
plt.hist(Weff,50,alpha=0.9,label='Model')
mu, sigma = Weffmean, Weffsigma # mean and standard deviation
Weffnew = np.random.normal(mu, sigma, 2000)
plt.hist(Weffnew,50,alpha=0.7,label='Model')
plt.ylabel("Frequency", fontsize=18)
plt.xlabel("Weff (m)", fontsize=18)

print np.corrcoef(Weff,Cins)
print np.corrcoef(Weff,Ach)

"""fignumber=2
plt.figure(fignumber)
model3 = plt.scatter(Weffnew,Cinsnew,s=300, marker="s",facecolors='none', edgecolors='r') 
ax = plt.gca()
plt.ylabel("Cins (F/m)", fontsize=18)
plt.xlabel("Weff (m)", fontsize=18)

fignumber=3
plt.figure(fignumber)
model3 = plt.scatter(Weffnew,Achnew,s=300, marker="s",facecolors='none', edgecolors='r') 
ax = plt.gca()
plt.ylabel("Ach (m^2)", fontsize=18)
plt.xlabel("Weff (m)", fontsize=18)"""

slope, intercept, r_value, p_value, std_err = stats.linregress(Weff,Cins)
print slope, intercept, r_value, p_value, std_err
mu, sigma = 0, np.std(slope*Weff+intercept-Cins) # mean and standard deviation
Cinsadd = np.random.normal(mu, sigma, 2000)
fignumber=2
plt.figure(fignumber)
line = slope*Weffnew+intercept+Cinsadd
model3 = plt.scatter(Weffnew,line,s=300, marker="s",facecolors='none', edgecolors='r') 

fignumber+=7
plt.figure(fignumber)
plt.hist(slope*Weff+intercept-Cins,50,alpha=0.9,label='Model') 

slope, intercept, r_value, p_value, std_err = stats.linregress(Weff,Ach)
print slope, intercept, r_value, p_value, std_err
mu, sigma = 0, np.std(slope*Weff+intercept-Ach) # mean and standard deviation
Achadd = np.random.normal(mu, sigma, 2000)
fignumber=3
plt.figure(fignumber)
line = slope*Weffnew+intercept+Achadd
model3 = plt.scatter(Weffnew,line,s=300, marker="s",facecolors='none', edgecolors='r') 

plt.show() 

