#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parsehspice

###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/idvgscaling14nm/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplateGEO114nm.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo4template_mc.nmos"#modelcard-109template.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)
Lparam_array = ['20e-9','50e-9','100e-9','250e-9','500e-9','1000e-9'] 
NFINparam = '1'
DEVTYPEparam = '1'
###################################
#prepare hspice


Lg_array    = 0.018e-6*np.array([1])#
Wt_array    = 0.0076e-6*np.array([1])
WbR_array   = 0.0076e-6/2*np.array([1])
WbL_array   = 0.0076e-6/2*np.array([1])
Hfin_array  = 0.042e-6*np.array([1])
tox_array   = 0.0008e-6*np.array([1])
Doping_FIN_array = 6e24*np.array([1])
PHIG = 4.20*np.array([1])#
RSHSparam = 300*np.array([1])#
RSHDparam = 300*np.array([1])#

Lg, Wt, WbR,WbL, Hfin, tox, Doping_FIN = sf.meshgrid2(Lg_array,Wt_array,WbR_array,WbL_array, Hfin_array, tox_array,Doping_FIN_array )
#
#prepare hspice
index=0
for L in Lg:
  preparehspice.preparehspiceidvgGEO1v2(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,str(L),str(Hfin[index]),str(Wt[index]),str(WbR[index]+WbL[index]+Wt[index]),str(tox[index]),str(Doping_FIN[index]),NFINparam,str(PHIG[index]),str(RSHSparam[index]),str(RSHDparam[index]))

  #run hspise
  filenameoutput = 'idvgaux'
  os.system('hspice ' + wheretosimpath +'idvgaux.sp -o ' + wheretosimpath+filenameoutput)

  #parse results
  parsehspice.parsehspicev1(wheretosimpath,filenameoutput,vds,str(L))


      
