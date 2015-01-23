#extract DIBL

import os
import numpy as np
import pylab
######################
filenameaux = 'montecarloresultsallfilevdsatmchspicereduced'
vdssatfile = open(filenameaux, 'r') 

filenameaux = 'montecarloresultsallfilevdlinmchspicereduced'
vdslinfile = open(filenameaux, 'r') 


montecarloresultsall = open('summarydevicemchspicereduced', 'w')
stringtoprint = 'Ioffsat'+' '+'Ionsat'+' '+'Lgsat'+' '+'Weffsat'+' '+'Cinssat'+' '+'Achsat'+' '+'Nfinsat'+' '+'vdsat'+' '+'Vthsat'+' gmaxsat'+' SSsat'+' Iofflin'+' '+'Ionlin'+' '+'Lglin'+' '+'Wefflin'+' '+'Cinslin'+' '+'Achlin'+' '+'Nfinlin'+' '+'vdlin'+' '+'Vthlin'+' gmaxlin'+' SSlin'+' DIBL'+'\n'
montecarloresultsall.write(stringtoprint)

count=0

for line in vdssatfile:
  #if count>10:
  #  break
  if count>0: 
    header = str.split(line)
    findWt = line.find("Weff") 
    stringtoprint = line[0:findWt] + ' '
    montecarloresultsall.write(stringtoprint)
    for line2 in vdslinfile:
      findfilename = line2.find(header[-1])    
      if findfilename > -1:
        header2 = str.split(line2)
        DIBL = -(float(header[8])-float(header2[8]))/0.81
        findWt = line2.find("Weff") 
        stringtoprint = line2[0:findWt] + ' '+str(DIBL)+'\n'
        montecarloresultsall.write(stringtoprint)           
    vdslinfile.seek(0)
  count+=1
  print count     
vdssatfile.close()
vdslinfile.close()
montecarloresultsall.close()        
