import numpy as np
import pylab

def plotinverter(wheretosimpath,namefile,namestrx,namestry,figurenumber,flaglog):
    target = open( wheretosimpath+namefile, 'r')
    header = str.split(target.readline())
    if len(header)>0:
      xindex =  header.index(namestrx)
      yindex =  header.index(namestry)
    
      datalist = np.loadtxt( wheretosimpath+namefile,skiprows = 1)
      pylab.figure(figurenumber)
      if flaglog:
        pylab.plot( datalist[:,xindex], abs(datalist[:,yindex]), lw=2 )      
        ax = pylab.gca()
        ax.set_yscale('log')      

      else:
        pylab.plot( datalist[:,xindex], (datalist[:,yindex]), lw=2 )
      target.close()
      
          
def plotsram(wheretosimpath,namefile,namestrx,namestry,figurenumber,flaglog):
    target = open( wheretosimpath+namefile, 'r')
    header = str.split(target.readline())
    if len(header)>0:
      xindex =  header.index(namestrx)
      yindex =  header.index(namestry)
    
      datalist = np.loadtxt( wheretosimpath+namefile,skiprows = 1)
      pylab.figure(figurenumber)
      
      Vdownarray = datalist[::-1,xindex]
      Vuparray = datalist[:,yindex]
      VlarrayVuparray = datalist[:,xindex]
      VlarrayVdownarray = datalist[::-1,yindex]
      
      if flaglog:
        pylab.plot( VlarrayVuparray, Vuparray, lw=2 ) 
        pylab.plot( VlarrayVdownarray, Vdownarray, lw=2 )      
        ax = pylab.gca()
        ax.set_yscale('log')      
      else:
        #pylab.subplot(1, 2, 1)
        pylab.plot( VlarrayVuparray, Vuparray, lw=2 ) 
        pylab.plot( VlarrayVdownarray, Vdownarray, lw=2 )   
      target.close()
      
      countdown = 0
      SNMarray = []
      for Vldown in VlarrayVdownarray:
        VL = Vldown
        VR = Vdownarray[countdown]
        countup = 0
        diff1=100000
        diff2=100000
        for Vlup in VlarrayVuparray:
          Vx =   Vlup + (-VL+VR)      
          diff = abs(Vx-Vuparray[countup])
          countup+=1
          if diff < diff1:
            diff1 = diff
            SNMaux = abs(Vlup-VL)
        SNMarray.append(SNMaux) 
        countdown+=1
        #break
      #pylab.figure(figurenumber+1)
      #pylab.plot( VlarrayVdownarray, SNMarray, 'o',lw=2 )
      snrlen =int(round(len(SNMarray) /2))

      SNM1 = np.amax(SNMarray[0:snrlen])
      SNM2 = np.amax(SNMarray[snrlen:-1])#round(snrlen/2)

      return min(SNM1, SNM2)
