#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parsehspice

###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/mchspice10nm/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplateGEO1.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo4template_extreme.nmos"#modelcard-109template.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)
Lparam_array = ['20e-9','50e-9','100e-9','250e-9','500e-9','1000e-9'] 
NFINparam = '1'
###################################
mu, sigma = 0.017e-6, 0.020e-6*0.1/2 # mean and standard deviation
Lg = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.005e-6, 0.0076e-6*0.1/2 # mean and standard deviation
Wt = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.005e-6/2, 0.0076e-6*0.1/2 # mean and standard deviation
WbR = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.005e-6/2, 0.0076e-6*0.1/2 # mean and standard deviation
WbL = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.042e-6, 0.042e-6*0.1/2 # mean and standard deviation
Hfin = np.random.normal(mu, sigma, 2000)

mu, sigma = 0.00077e-6, 0.0008e-6*0.1/2 # mean and standard deviation
tox = np.random.normal(mu, sigma, 2000)

mu, sigma = 7e24, 7e24*0.1/2 # mean and standard deviation
Doping_FIN = np.random.normal(mu, sigma, 2000)


#
#prepare hspice
os.system('hspice -C')

index=0
for L in Lg:
  preparehspice.preparehspiceidvgGEO1(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,str(L),str(Hfin[index]),str(Wt[index]),str(WbR[index]+WbL[index]+Wt[index]),str(tox[index]),str(Doping_FIN[index]),NFINparam)

  #run hspise
  filenameoutput = 'idvgaux'
  os.system('hspice -C ' + wheretosimpath +'idvgaux.sp -o ' + wheretosimpath+filenameoutput)

  #parse results
  parsehspice.parsehspicev2(wheretosimpath,filenameoutput,vds,str(L),str(Wt[index]),str(WbR[index]),str(WbL[index]),str(Hfin[index]),str(tox[index]),str(Doping_FIN[index]))
  index +=1
  print index
 
os.system('hspice -C -K')
      
