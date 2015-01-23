#parse data output from inverter

###################################d
def parseinvv1(wheretosimpath,filenameoutput,finalnamedata):
  filepath =wheretosimpath

  filenameaux = filenameoutput+'.lis'
  outputfiletoread = open(filepath+filenameaux, 'r') 
  #state machine
  #00: searching for "volt    voltage"
  #01: searching for x1.d ->open file
  #02: search y, meanhwile parsing data id-vg
  state=0
  outputnamefiles = []
  countruns=0
  for line in outputfiletoread:
    #print line
    if (state==0):
      if (line.find("volt    voltage")>-1):
        state=1
        filenameaux =''      
    elif (state==1):
      if line.find("x1.d")>-1:
        filenameaux = finalnamedata +'INVERTER'+str(countruns)+'.out'
        outputnamefiles.append(filenameaux)
        countruns+=1
        #print filenameaux 
        hspicefile = open(filepath+filenameaux, 'w') 
        hspicefile.write('Vin Vout\n') 
        state=2
    elif (state==2):
      if line.find("y")>-1 :
        state=0
        hspicefile.close()
      else:
        stringinline = str.split(line)
        hspicefile.write(stringinline[0]+' '+stringinline[1]+'\n')  
  return outputnamefiles         
  

