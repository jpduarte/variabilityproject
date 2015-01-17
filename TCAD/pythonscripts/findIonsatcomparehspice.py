#findIonsat.py: read all the files and plot Id-Vg for Vd=0.86 V, in linear and log scale, compare with with two folders
import os
import numpy as np
import pylab

################################## 
root_dir = '/users/jpduarte/research/variabilityproject/netlist_modelcards/cornerhspice'
data_files = [(x[0], x[2]) for x in os.walk(root_dir)]
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number

#filenames = data_files[0][1]

factorIDS = 1e6/(42e-3*2+7.6e-3)

count = 1
for namefile in data_files[0][1]:
  #if count>2:
	#	break
  flagfiledata =  namefile.find("DATAREADY") #if it is found "DATAREADY" it returns the position, this can be used for 
  flagfileIon =  namefile.find("vd0.86") 
  if ((flagfiledata > 0) and (flagfileIon > 0)):
    #print namefile
    filenameaux = root_dir +'/'+ namefile
    target = open( filenameaux, 'r')
    fistline = target.readline()
    header = str.split(fistline)
    findvgindex = fistline.find("gateOuterVoltage")     
    if (len(header)>0 and findvgindex>-1):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
    
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      
      pylab.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS, lw=2 )
      target.close()
      count +=1
      print count

    
pylab.xlabel("VG (V)", fontsize=20)
pylab.ylabel("ID,SAT (uA/um)", fontsize=20)
ax = pylab.gca()
pylab.xlim([0,0.86])
#pylab.ylim([0,1000])
#ax.set_yscale('log')  
#pylab.savefig('IonsatvsVglin', dpi=600, bbox_inches='tight')  
ax.set_yscale('log')  
pylab.xlim([0,0.86])
pylab.ylim([1e-2,1e4])
#pylab.savefig('IonsatvsVglog', dpi=600, bbox_inches='tight') 
pylab.show()
