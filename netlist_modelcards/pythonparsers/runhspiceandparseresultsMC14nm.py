#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parsehspice

###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/mchspice10nmv3/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplateGEO114nm.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo4template_mc.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)

NFINparam = '1'
###################################
mu, sigma = 0.016e-6, 0.020e-6*0.12/2 # mean and standard deviation
Lg = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.006e-6, 1e-9 # mean and standard deviation
Wt = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.001e-6/2, 0.0005e-6 # mean and standard deviation
WbR = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.001e-6/2, 0.0005e-6 # mean and standard deviation
WbL = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.042e-6, 1e-9 # mean and standard deviation
Hfin = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.0007e-6, 0.0008e-6*0.08/2 # mean and standard deviation
tox = np.random.normal(mu, sigma, 2000)

mu, sigma = 7e24, 7e24*0.1/2 # mean and standard deviation
Doping_FIN = np.random.normal(mu, sigma, 2000)

mu, sigma = 4.20, 0.07/2 # mean and standard deviation,v1 2sigma=0.03, v2 0.045/2, v3 0.07
PHIG = np.random.normal(mu, sigma, 2000)

mu, sigma = 300, 50 # mean and standard deviation
RSHSparam = np.random.normal(mu, sigma, 2000)

mu, sigma = 300, 50 # mean and standard deviation
RSHDparam = np.random.normal(mu, sigma, 2000)


#
#prepare hspice
os.system('hspice -C')

index=0
for L in Lg:
  preparehspice.preparehspiceidvgGEO1v2(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,str(L),str(Hfin[index]),str(Wt[index]),str(WbR[index]+WbL[index]+Wt[index]),str(tox[index]),str(Doping_FIN[index]),NFINparam,str(PHIG[index]),str(RSHSparam[index]),str(RSHDparam[index]))

  #run hspise
  filenameoutput = 'idvgaux'
  os.system('hspice -C ' + wheretosimpath +'idvgaux.sp -o ' + wheretosimpath+filenameoutput)

  #parse results
  parsehspice.parsehspicev2(wheretosimpath,filenameoutput,vds,str(L),str(Wt[index]),str(WbR[index]),str(WbL[index]),str(Hfin[index]),str(tox[index]),str(Doping_FIN[index]))
  index +=1
  print index

 
os.system('hspice -C -K')
      
