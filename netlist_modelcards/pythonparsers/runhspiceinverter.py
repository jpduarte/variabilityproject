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
###################################
wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/SRAM/'
"""templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplate.sp"
#modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
#modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109-geo1template.nmos"#modelcard-109template.nmos"
vgssample = 50
vdssample = 2
vgs = np.linspace(0,0.86,vgssample)
vds = np.linspace(0.05,0.86,vdssample)
Lparam_array = ['20e-9','50e-9','100e-9','250e-9','500e-9','1000e-9'] 
NFINparam = '1'"""
###################################
#prepare hspice
#for Lparam in Lparam_array:
#preparehspice.preparehspiceidvg(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,NFINparam)

#run hspise
filetorun = 'inverter1'
os.system('hspice ' + wheretosimpath +filetorun+'.sp -o ' + wheretosimpath+filetorun)

#parse results
finalnamedata = 'VTC'
outputnamefiles = parseinverterout.parseinvv1(wheretosimpath,filetorun,finalnamedata,"vout")#return list with names of files 

figurenumber=1
for namefile in outputnamefiles:
  plotgeneral.plotinverter(wheretosimpath,namefile,'Vin','Vout',figurenumber,0)
  figurenumber+=1
  
pylab.show()  
  
  


      
