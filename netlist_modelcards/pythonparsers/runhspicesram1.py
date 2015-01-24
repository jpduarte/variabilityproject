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
###################################
#prepare hspice
#for Lparam in Lparam_array:
#preparehspice.preparehspiceidvg(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,NFINparam)

#run hspise
filetorun = 'sram1'
os.system('hspice ' + wheretosimpath +filetorun+'.sp -o ' + wheretosimpath+filetorun)

#parse results
finalnamedata = 'sram'
outputnamefiles = parseinverterout.parseinvv1(wheretosimpath,filetorun,finalnamedata,"vl2")#return list with names of files 
print outputnamefiles
figurenumber=1
for namefile in outputnamefiles:
  plotgeneral.plotinverter(wheretosimpath,namefile,'Vin','Vout',figurenumber,0)
  figurenumber+=1
  
figurenumber=1
for namefile in outputnamefiles:
  plotgeneral.plotinverter(wheretosimpath,namefile,'Vout','Vin',figurenumber,0)
  figurenumber+=1  
  
pylab.show()  
  
  


      
