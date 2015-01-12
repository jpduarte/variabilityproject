#extract DIBL for file sin v9 folder

import os
import numpy as np
import pylab
######################
filenameaux = 'montecarloresultsallfilevdsat2'
vdssatfile = open(filenameaux, 'r') 

filenameaux = 'montecarloresultsallfilevdlin2'
vdslinfile = open(filenameaux, 'r') 


montecarloresultsall = open('summarydeviceMC2', 'w')
stringtoprint = 'Ioffsat'+' '+'Ionsat'+' '+'Wtsat'+' '+'Lgsat'+' '+'Hfinsat'+' '+'toxsat'+' '+'WbRsat'+' '+'WbLsat' +' '+'Nfinsat'+' '+'WFsat'+' '+'vdsat'+' '+'Vthsat'+' gmaxsat'+' SSsat'+' Iofflin'+' '+'Ionlin'+' '+'Wtlin'+' '+'Lglin'+' '+'Hfinlin'+' '+'toxlin'+' '+'WbRlin'+' '+'WbLlin' +' '+'Nfinlin'+' '+'WFlin'+' '+'vdlin'+' '+'Vthlin'+' gmaxlin'+' SSlin'+' DIBL'+'\n'
montecarloresultsall.write(stringtoprint)

count=0

for line in vdssatfile:
  #if count>3:
  #  break
  if count>0: 
    header = str.split(line)
    findWt = line.find("Wt") 
    stringtoprint = line[0:findWt] + ' '
    montecarloresultsall.write(stringtoprint)
    for line2 in vdslinfile:
      findfilename = line2.find(header[-1])    
      if findfilename > -1:
        header2 = str.split(line2)
        DIBL = -(float(header[11])-float(header2[11]))/0.86
        #print DIBL,header[11],header2[11]
        findWt = line2.find("Wt") 
        stringtoprint = line2[0:findWt] + ' '+str(DIBL)+'\n'
        montecarloresultsall.write(stringtoprint)           
    vdslinfile.seek(0)
  count+=1
  print count     
vdssatfile.close()
vdslinfile.close()
montecarloresultsall.close()        
