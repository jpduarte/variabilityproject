#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parsehspice

###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/mchspicereduced/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplateGEO4.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo4template_extreme_MCreduced.nmos"#modelcard-109template.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)
NFINparam = '1'
###################################
"""Weff
sigma = 4.2412981124e-09
mean = 9.19717434862e-08

Ach-Weff relations
5.89760127488e-09 -6.36648180252e-17 1.46795688713e-17

Cins-Weff relations
0.0416526092526 1.44347020205e-10 1.96760768758e-10


"""
mu, sigma = 0.020e-6, 0.020e-6*0.1/2 # mean and standard deviation
Lg = np.random.normal(mu, sigma, 2000)

mu, sigma = 9.19717434862e-08, 4.2412981124e-09# mean and standard deviation
Weff = np.random.normal(mu, sigma, 2000)

mu, sigma = 0, 1.96760768758e-10 # mean and standard deviation
Cinsadded = np.random.normal(mu, sigma, 2000)
Cins = 0.0416526092526*Weff + 1.44347020205e-10+Cinsadded

mu, sigma = 0, 1.46795688713e-17 # mean and standard deviation
Achadded = np.random.normal(mu, sigma, 2000)
Ach = 5.89760127488e-09*Weff + -6.36648180252e-17 +Achadded

mu, sigma = 6e24, 6e24*0.1/2 # mean and standard deviation
Doping_FIN = np.random.normal(mu, sigma, 2000)


#
#prepare hspice
os.system('hspice -C')

index=0
for L in Lg:
  preparehspice.preparehspiceidvgGEO4(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,str(L),str(Ach[index]),str(Cins[index]),str(Weff[index]),str(Doping_FIN[index]),NFINparam)
  #Lparam,Ach_UFCMparam,Cins_UFCMparam,W_UFCMparam,NBODYparam,NFINparam

  #run hspise
  filenameoutput = 'idvgaux'
  os.system('hspice -C ' + wheretosimpath +'idvgaux.sp -o ' + wheretosimpath+filenameoutput)

  #parse results
  parsehspice.parsehspicev3(wheretosimpath,filenameoutput,vds,str(L),str(Weff[index]),str(Cins[index]),str(Ach[index]),str(Doping_FIN[index]))#Lg, Weff, Cins, Ach, Doping_FIN 
  index +=1
  print index
 
os.system('hspice -C -K')
      
