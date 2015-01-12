#plot Id-Vg scaling, TCAD versus Hspice
import os
import numpy as np
import pylab

root_dir = '../data/v7'
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
      pylab.figure(1)
      pylab.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS,'o',markerfacecolor=(1, 1, 1, 1),lw=1, color='k' )      
      pylab.figure(2)
      pylab.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS,'o',markerfacecolor=(1, 1, 1, 1),lw=1, color='k' )
      target.close()
      count +=1
      print count

    
################################## 

root_dir = '../../netlist_modelcards/idvgscaling'
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
    filenameaux = root_dir +'/'+ namefile
    target = open( filenameaux, 'r')
    fistline = target.readline()
    header = str.split(fistline)
    findvgindex = fistline.find("gateOuterVoltage")   
    if (len(header)>0 and findvgindex>-1):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
    
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      pylab.figure(1)
      pylab.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS, lw=2, color='r' )
      pylab.figure(2)
      pylab.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS, lw=2, color='r' )      
      target.close()
      count +=1
      print count

pylab.figure(1)  
pylab.xlabel("VG (V)", fontsize=20)
pylab.ylabel("ID,LIN (uA/um)", fontsize=20)
ax = pylab.gca()
pylab.xlim([0,0.86])
#pylab.savefig('IonlinvsVgLGscalinglinGEO1', dpi=600, bbox_inches='tight') 


pylab.figure(2)  
pylab.xlabel("VG (V)", fontsize=20)
pylab.ylabel("ID,LIN (uA/um)", fontsize=20)
ax = pylab.gca()
ax.set_yscale('log')  
pylab.xlim([0,0.86])
pylab.ylim([1e-4,1e4])
#pylab.savefig('IonlinvsVgLGscalinglogGEO1', dpi=600, bbox_inches='tight') 



pylab.show()
