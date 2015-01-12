#this file parse the data from hspice .out results
import re
import numpy as np
###################################
def parsehspicev1(wheretosimpath,filenameoutput,vds,Lparam):
  filepath =wheretosimpath

  filenameaux = filenameoutput+'.lis'
  outputfiletoread = open(filepath+filenameaux, 'r') 
  #state machine
  #00: searching for "volt    current"
  #01: searching for x1.d and parsing parameters, +EOT ='tox_n+tox_mc*tox_n*0.1', tox_mc,1:eot           =  8.4856E-10  2:eot           =  8.2360E-10
  #02: search y, meanhwile parsing data id-vg
  x = [':']
  delimiters = "1:", "=", "\n", "         "
  regexPattern = '|'.join(map(re.escape, delimiters))
  state=0

  countruns=0
  for line in outputfiletoread:
    #print line
    if (state==0):
      if (line.find("volt    current")>-1):
        state=1
        filenameaux =''      
    elif (state==1):
      if line.find("x1.d")>-1:
        filenameaux = filenameoutput + 'Lg'+Lparam+'vd'+str(vds[countruns]) + 'DATAREADY.out'
        countruns+=1
        print filenameaux 
        hspicefile = open(filepath+filenameaux, 'w') 
        hspicefile.write('gateOuterVoltage drainTotalCurrent\n') 
        state=2
    elif (state==2):
      if line.find("y")>-1 :
        state=0
        hspicefile.close()
      else:
        stringinline = str.split(line)
        hspicefile.write(stringinline[0]+' '+stringinline[1]+'\n')
        
def parsehspicev2(wheretosimpath,filenameoutput,vds,Lg, Wt, WbR,WbL, Hfin, tox, Doping_FIN ):
  filepath =wheretosimpath

  filenameaux = filenameoutput+'.lis'
  outputfiletoread = open(filepath+filenameaux, 'r') 
  #state machine
  #00: searching for "volt    current"
  #01: searching for x1.d and parsing parameters, +EOT ='tox_n+tox_mc*tox_n*0.1', tox_mc,1:eot           =  8.4856E-10  2:eot           =  8.2360E-10
  #02: search y, meanhwile parsing data id-vg
  x = [':']
  delimiters = "1:", "=", "\n", "         "
  regexPattern = '|'.join(map(re.escape, delimiters))
  state=0

  countruns=0
  for line in outputfiletoread:
    #print line
    if (state==0):
      if (line.find("volt    current")>-1):
        state=1
        filenameaux =''      
    elif (state==1):
      if line.find("x1.d")>-1:
        filenameaux = filenameoutput + 'Wt'+Wt+'Lg'+Lg+'Hfin'+Hfin+'tox'+tox+'WbR'+WbR+'WbL'+WbL+'Nfin'+Doping_FIN +'vd'+str(vds[countruns]) + 'DATAREADY.out'
        countruns+=1
        #print filenameaux 
        hspicefile = open(filepath+filenameaux, 'w') 
        hspicefile.write('gateOuterVoltage drainTotalCurrent\n') 
        state=2
    elif (state==2):
      if line.find("y")>-1 :
        state=0
        hspicefile.close()
      else:
        stringinline = str.split(line)
        hspicefile.write(stringinline[0]+' '+stringinline[1]+'\n')        
