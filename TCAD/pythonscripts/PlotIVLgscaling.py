#plot Id-Vg scaling, TCAD versus Hspice
import os
import numpy as np
import matplotlib.pyplot as plt

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
      plt.subplot(1, 2, 2)
      tcad1 = plt.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS,'o',markerfacecolor=(1, 1, 1, 1),lw=1, color='k' )      
      plt.subplot(1, 2, 1)
      tcad2 = plt.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS,'o',markerfacecolor=(1, 1, 1, 1),lw=1, color='k' )
      target.close()
      count +=1
      print count

    
################################## 

root_dir = '../../netlist_modelcards/idvgscaling14nm'
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
      plt.subplot(1, 2, 2)
      model1 = plt.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS, lw=2, color='r' )
      plt.subplot(1, 2, 1)
      model2 = plt.plot( datalist[:,vgindex], datalist[:,Idsindex]*factorIDS, lw=2, color='r' )      
      target.close()
      count +=1
      print count

plt.subplot(1, 2, 2) 
plt.xlabel("VG (V)", fontsize=20)
#plt.ylabel("ID,LIN (uA/um)", fontsize=20)
ax = plt.gca()
plt.xlim([0,0.86])
plt.ylim([-100,1400])
ax.arrow(0.6, 0, -0.25, 400, head_width=0.05, head_length=25, fc='k', ec='k')
ax.text(0.01, 600, r'LG = 1000nm~', fontsize=15)
ax.text(0.2, 500, r'20nm', fontsize=15)
#plt.savefig('IonlinvsVgLGscalinglinGEO1', dpi=600, bbox_inches='tight') 
plt.yticks(fontsize = 16) 
plt.legend([ tcad1,model1],['TCAD','Model'],loc=2,prop={'size':20})
plt.subplot(1, 2, 1) 
plt.xlabel("VG (V)", fontsize=20)
plt.ylabel("Drain Current (uA/um), VDS=0.86 V", fontsize=20)
ax = plt.gca()
ax.set_yscale('log')  
plt.xlim([0,0.86])
plt.ylim([1e-4,1e4])
#plt.xticks(fontsize = 15) 
plt.yticks(fontsize = 18) 
#
ax.arrow(0.6, 10, -0.25, 400, head_width=0.05, head_length=50, fc='k', ec='k')

plt.legend([ tcad2,model2],['TCAD','Model'],loc=4,prop={'size':20})
ax.text(0.01, 3000, r'LG = 1000nm~', fontsize=15)
ax.text(0.2, 1000, r'20nm', fontsize=15)
#plt.savefig('DrainCurrentScaling', dpi=600, bbox_inches='tight') 
plt.show()
