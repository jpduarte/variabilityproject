#extracresults script from I-V results
#Juan Duarte

import os
import numpy as np
import pylab
######################
def findIds(vg_array,ids_array,vgref):
  index = 0 
  indexids = -1
  for vg in vg_array:
    if abs(vg-vgref)<1e-3:
      indexids = index
    index+=1
  return ids_array[indexids]

def findvth(vg_array,ids_array,idsref):
#find vg for a given current level
  index = 0 
  diff = 1000
  for ids in ids_array:
    if abs(ids-idsref)<diff:
      indexids = index
      diff = abs(ids-idsref)
    index+=1
  return vg_array[indexids]

def findgmax(vg_array,ids_array,dervnumber):
#find gmax, or max derivative, dervnumber=1 is for gmax
  dI = np.diff(ids_array,n=dervnumber) 
  dV = np.diff(vg_array,n=dervnumber) 
  gmax = 0
  index=0
  for dVi in dV:
    if abs(dVi)<1e-4:
      gmaux = 0
    else:
      gmaux = dI[index]/dVi
    if gmaux>gmax:
      gmax = gmaux
    index+=1
  return gmaux

def findSS(vg_array,ids_array,vgi,level):
#find SS with the current at given vgi wrt current*level
  index = 0 
  for vg in vg_array:
    if abs(vg-vgi)<1e-3:
      indexids1 = index
    index+=1
  Ioff = ids_array[indexids1]
  indexids2 = 0
  index = 0 
  diff = 1000
  for ids in ids_array:
    if (abs(ids-Ioff*level)<diff and index>indexids1):#abs(np.log(abs(ids))-np.log(abs(Ioff*1e3)))<diff 
      indexids2 = index
      diff = abs(ids-Ioff*level)
    index+=1
  if indexids2 == 0:
    SS = 0
  else:
    SS = np.log(10)*(vg_array[indexids2]-vgi)/(np.log(ids_array[indexids2]/Ioff))
  return SS
######################
root_dir = '/home/jpduarte/STDB/FINFETSRC2014/v10'
data_files = [(x[0], x[2]) for x in os.walk(root_dir)]
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number
#filenames = data_files[0][1]

factorIDS = 1e6/(42e-3*2+7.6e-3)#this factor is for normalize

montecarloresultsall = open('montecarloresultsallfile2', 'w')
stringtoprint = 'Ioff'+' '+'Ion'+' '+'Wt'+' '+'Lg'+' '+'Hfin'+' '+'tox'+' '+'WbR'+' '+'WbL' +' '+'Nfin'+' WF '+'vd'+' '+'Vth'+' gmax'+' SS'+' filename'+'\n'
montecarloresultsall.write(stringtoprint)

montecarloresultsallvdlin = open('montecarloresultsallfilevdlin2', 'w')
montecarloresultsallvdlin.write(stringtoprint)
montecarloresultsallvdmed = open('montecarloresultsallfilevdmed2', 'w')
montecarloresultsallvdmed.write(stringtoprint)
montecarloresultsallvdsat = open('montecarloresultsallfilevdsat2', 'w')
montecarloresultsallvdsat.write(stringtoprint)

count = 1
for namefile in data_files[0][1]:
  #if count>10:
	#	break
  flagfiledata =  namefile.find("DATAREADY") #if it is found "DATAREADY" it returns the position, this can be used for 
  if ((flagfiledata > 0)):
    indexWt = namefile.find("Wt") 
    indexLg = namefile.find("Lg") 
    indexHfin = namefile.find("Hfin") 
    indextox = namefile.find("tox") 
    indexWbR = namefile.find("WbR") 
    indexWbL = namefile.find("WbL") 
    indexNfin = namefile.find("Nfin") 
    indexWF = namefile.find("WF") 
    indexvd = namefile.find("vd") 

    Wt =  namefile[indexWt+2:indexLg]
    Lg =  namefile[indexLg+2:indexHfin]
    Hfin =  namefile[indexHfin+4:indextox]
    tox =  namefile[indextox+3:indexWbR]
    WbR =  namefile[indexWbR+3:indexWbL]
    WbL =  namefile[indexWbL+3:indexNfin]
    Nfin =  namefile[indexNfin+4:indexWF]
    WF =  namefile[indexWF+2:indexvd]    
    vd =  namefile[indexvd+2:flagfiledata]
      
    nametoprint =  namefile[indexWt:indexvd]


    filenameaux = root_dir +'/'+ namefile
    target = open( filenameaux, 'r')
    fistline = target.readline()
    header = str.split(fistline)
    findvgindex = fistline.find("gateOuterVoltage") 
    if (len(header)>0 and findvgindex>0):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      Ion = np.amax(datalist[:,Idsindex])
      Ioff = findIds(datalist[:,vgindex],datalist[:,Idsindex],0) 
      Vth = findvth(datalist[:,vgindex],datalist[:,Idsindex],300*1e-9*(42e-9*2+7.6e-9)/(20e-9))#300nA*W/L for current reference in Vth
      gmax = findgmax(datalist[:,vgindex],datalist[:,Idsindex],1)
      SS = findSS(datalist[:,vgindex],datalist[:,Idsindex],0.0,1e3)
      stringtoprint = str(Ioff)+' '+str(Ion)+' '+Wt+' '+Lg+' '+Hfin+' '+tox+' '+WbR+' '+WbL +' '+Nfin+' '+WF+' '+vd+' '+str(Vth)+' '+str(gmax)+' '+str(SS)+' '+ nametoprint+'\n' 
      montecarloresultsall.write(stringtoprint)

      if vd =='0.05': 
        montecarloresultsallvdlin.write(stringtoprint)
      if vd =='0.43': 
        montecarloresultsallvdmed.write(stringtoprint)
      if vd =='0.86': 
        montecarloresultsallvdsat.write(stringtoprint)
      target.close()
      count +=1
      print count

montecarloresultsall.close()
montecarloresultsallvdlin.close()
montecarloresultsallvdmed.close()
montecarloresultsallvdsat.close()

