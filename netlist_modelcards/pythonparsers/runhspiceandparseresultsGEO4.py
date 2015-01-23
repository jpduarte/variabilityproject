#this file parse the data from hspice .out results
import re
import os
import numpy as np
import supportfunctions as sf
import shutil
import preparehspice
import parsehspice

###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/idvgscalinggeo4/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplateGEO4.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo4template_extreme_MCreduced.nmos"#modelcard-109template.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)
Lparam_array = ['20e-9','50e-9','100e-9','250e-9','500e-9','1000e-9'] 
NFINparam = '1'
DEVTYPEparam = '1'
Ach = '4.78747853612e-16'
Cins = '3.97521011392e-11'
Weff = '9.19717434862e-08'
Doping_FIN = '6e24'
"""
Cins
2.6443138247e-12
3.97521011392e-11
Weff
4.2412981124e-09
9.19717434862e-08
Ach
2.90028305832e-17
4.78747853612e-16
"""
###################################
#prepare hspice
for Lparam in Lparam_array:
  #preparehspice.preparehspiceidvg(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,NFINparam,DEVTYPEparam)
  preparehspice.preparehspiceidvgGEO4(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,Ach,Cins,Weff,Doping_FIN,NFINparam)
  #run hspise
  filenameoutput = 'idvgaux'
  os.system('hspice ' + wheretosimpath +'idvgaux.sp -o ' + wheretosimpath+filenameoutput)

  #parse results
  parsehspice.parsehspicev1(wheretosimpath,filenameoutput,vds,Lparam)


      
