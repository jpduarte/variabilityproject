import numpy as np
import pylab

def plotinverter(wheretosimpath,namefile,namestrx,namestry,figurenumber):
    target = open( wheretosimpath+namefile, 'r')
    header = str.split(target.readline())
    if len(header)>0:
      xindex =  header.index(namestrx)
      yindex =  header.index(namestry)
    
      datalist = np.loadtxt( wheretosimpath+namefile,skiprows = 1)
      pylab.figure(figurenumber)
      pylab.plot( datalist[:,xindex], datalist[:,yindex], lw=2 )
      target.close()

