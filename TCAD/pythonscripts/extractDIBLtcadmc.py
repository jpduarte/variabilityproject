#extract DIBL for file sin v9 folder

import os
import numpy as np
import pylab
######################
filenameaux = 'montecarloresultsallfilevdsattcadmc'
vdssatfile = open(filenameaux, 'r') 

filenameaux = 'montecarloresultsallfilevdlintcadmc'
vdslinfile = open(filenameaux, 'r') 


montecarloresultsall = open('summarydevicetcadmc', 'w')
stringtoprint = 'Ioffsat'+' '+'Ionsat'+' '+'Wtsat'+' '+'Lgsat'+' '+'Hfinsat'+' '+'toxsat'+' '+'WbRsat'+' '+'WbLsat' +' '+'Nfinsat'+' '+'vdsat'+' '+'Vthsat'+' gmaxsat'+' SSsat'+' Iofflin'+' '+'Ionlin'+' '+'Wtlin'+' '+'Lglin'+' '+'Hfinlin'+' '+'toxlin'+' '+'WbRlin'+' '+'WbLlin' +' '+'Nfinlin'+' '+'vdlin'+' '+'Vthlin'+' gmaxlin'+' SSlin'+' DIBL'+'\n'
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
        DIBL = -(float(header[10])-float(header2[10]))/0.81
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
