#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parseinverterout
import plotgeneral
import pylab
import scipy.stats as stats
###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/SRAM/'
###################################
#prepare hspice
#for Lparam in Lparam_array:
#preparehspice.preparehspiceidvg(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,NFINparam)

numbertrans = ['1','2','3','4','5','6']
PHIarra = [4.20,4.93,4.20,4.93,4.20,4.20]

totalsim=2000
count=  0
SNR = []
os.system('hspice -C')
while count<totalsim:
  print count
  shutil.copyfile('../SRAM/paramsram.include','../SRAM/paramsramaux.include')
  for number in numbertrans:
    mu, sigma = 0.016e-6, 0.020e-6*0.12/2 # mean and standard deviation
    Lg = np.random.normal(mu, sigma, 1)

    mu, sigma = 0.006e-6, 1e-9 # mean and standard deviation
    Wt = np.random.normal(mu, sigma, 1)

    mu, sigma = 0.001e-6/2, 0.0005e-6 # mean and standard deviation
    WbR = np.random.normal(mu, sigma, 1)

    mu, sigma = 0.001e-6/2, 0.0005e-6 # mean and standard deviation
    WbL = np.random.normal(mu, sigma, 1)

    mu, sigma = 0.042e-6, 1e-9 # mean and standard deviation
    Hfin = np.random.normal(mu, sigma, 1)

    mu, sigma = 0.0007e-6, 0.0008e-6*0.08/2 # mean and standard deviation
    tox = np.random.normal(mu, sigma, 1)

    mu, sigma = 7e24, 7e24*0.1/2 # mean and standard deviation
    Doping_FIN = np.random.normal(mu, sigma, 1)

    if number == '2' or number == '4':
      mu = 4.93
    else:
      mu = 4.2
    
    sigma = 0.07/2 # mean and standard deviation,v1 2sigma=0.03, v2 0.045/2, v3 0.07, 
    PHIG = np.random.normal(mu, sigma, 1)

    mu, sigma = 300, 50 # mean and standard deviation
    RSHSparam = np.random.normal(mu, sigma, 1)

    mu, sigma = 300, 50 # mean and standard deviation
    RSHDparam = np.random.normal(mu, sigma, 1)
    
    sf.inplace_change('../SRAM/paramsramaux.include', 'Lparam'+number+'value', str(Lg[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'HFINparam'+number+'value', str(Hfin[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'TFIN_TOPparam'+number+'value', str(Wt[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'TFIN_BASEparam'+number+'value', str(Wt[0]+WbR[0]+WbL[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'EOTparam'+number+'value', str(tox[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'NBODYparam'+number+'value', str(Doping_FIN[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'NFINparam'+number+'value', str(1))
    sf.inplace_change('../SRAM/paramsramaux.include', 'PHIGparam'+number+'value', str(PHIG[0]))    
    sf.inplace_change('../SRAM/paramsramaux.include', 'RSHSparam'+number+'value', str(RSHSparam[0]))
    sf.inplace_change('../SRAM/paramsramaux.include', 'RSHDparam'+number+'value', str(RSHDparam[0]))
  count+=1



#run hspise
  filetorun = 'sram2'
  os.system('hspice -C ' + wheretosimpath +filetorun+'.sp -o ' + wheretosimpath+filetorun)

#parse results
  finalnamedata = 'sram2'+str(count)
  outputnamefiles = parseinverterout.parseinvv1(wheretosimpath,filetorun,finalnamedata,"vl2")#return list with names of files 
  #print outputnamefiles
  figurenumber=1
  
  for namefile in outputnamefiles:
    SRNaux = plotgeneral.plotsram(wheretosimpath,namefile,'Vin','Vout',figurenumber,0)
    SNR.append(SRNaux)
    figurenumber+=1
  
pylab.figure(1)

pylab.xlabel("VL (V)", fontsize=21)
pylab.ylabel("VR (V)", fontsize=21)
ax = pylab.gca()
pylab.xticks(fontsize=21) 
pylab.yticks(fontsize=21) 

#ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
#ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
pylab.figure(2)
stats.probplot(SNR,  plot=pylab)
pylab.xlabel("Quantalies", fontsize=21)
pylab.ylabel("SNM (V)", fontsize=21)
ax = pylab.gca()
pylab.xticks(fontsize=21) 
pylab.yticks(fontsize=21) 



os.system('hspice -C -K')
pylab.show()  
  
  


      
