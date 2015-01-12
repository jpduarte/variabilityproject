#this file parse the data from hspice .out results
import re
filepath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/montecarlohspice/devicesimulations/'

filenameaux = 'idvgnmos.out'
outputfiletoread = open(filepath+filenameaux, 'r') 
montecarloruns = 3 
#state machine
#00: searching for MONTE CARLO PARAMETER DEFINITIONS
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
    if (line.find("MONTE CARLO PARAMETER DEFINITIONS")>-1):
      state=1
      filenameaux =''
      countruns+=1
      if countruns<=montecarloruns:
        vdsvalue = 'vds0.05'
      else:
        vdsvalue = 'vds0.86'
      
  elif (state==1):
    if line.find(":")>-1 :

      stringinline = re.split(regexPattern, line)#str.split(line,[':','1'])
      variablename = stringinline[1]
      variablevalue = stringinline[3]
      filenameaux = filenameaux + variablename+variablevalue.replace(" ","")
    if line.find("x1.d")>-1:
      filenameaux = filenameaux + vdsvalue + 'DEVICE.out'
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
