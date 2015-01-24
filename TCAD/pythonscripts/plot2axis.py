#findIonsat.py: read all the files and plot Id-Vg for Vd=0.86 V, in linear and log scale
import os
import numpy as np
import matplotlib.pyplot as plt


root_dir = '/users/jpduarte/research/variabilityproject/TCAD/data/v6'
data_files = [(x[0], x[2]) for x in os.walk(root_dir)]
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number

#filenames = data_files[0][1]
def findIds(vg_array,ids_array,vgref):
  index = 0 
  indexids = -1
  for vg in vg_array:
    if abs(vg-vgref)<1e-3:
      indexids = index
    index+=1
  return ids_array[indexids],indexids

factorIDS = 1e6/(42e-3*2+7.6e-3)

count = 1
for namefile in data_files[0][1]:
  #if count>1:
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
    if (len(header)>0 and findvgindex>0):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
      
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      xx,indexstart = findIds(datalist[:,vgindex],datalist[:,Idsindex],0)
      tcadsim = plt.plot( datalist[indexstart:-1,vgindex], datalist[indexstart:-1,Idsindex]*factorIDS, lw=1, color='k' )
      target.close()
      count +=1
      print count


################
#plot worst SS

root_dir = '/users/jpduarte/research/variabilityproject/netlist_modelcards'
data_files = ['Wt0.00836Lg0.018Hfin0.0462tox0.00088WbR0.00418WbL0.00418Nfin5.4e+18vd0.86DATAREADY']
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number

#filenames = data_files[0][1]

factorIDS = 1e6/(42e-3*2+7.6e-3)

count = 1
for namefile in data_files:
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
    if (len(header)>0 and findvgindex>0):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
      
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      xx,indexstart = findIds(datalist[:,vgindex],datalist[:,Idsindex],0)
      worstSS = plt.plot( datalist[indexstart:-1,vgindex], datalist[indexstart:-1,Idsindex]*factorIDS, '-s',lw=1, color='g',markersize=6 )
      target.close()
      count +=1
      print count
      
################
#plot best SS

root_dir = '/users/jpduarte/research/variabilityproject/TCAD/data/v6'
data_files = ['finfet_6Wt0.00684Lg0.022Hfin0.0378tox0.00088WbR0.00418WbL0.00418Nfin6.6e+18vd0.86DATAREADYDATAREADY']
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number

#filenames = data_files[0][1]

factorIDS = 1e6/(42e-3*2+7.6e-3)

count = 1
for namefile in data_files:
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
    if (len(header)>0 and findvgindex>0):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
      
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      xx,indexstart = findIds(datalist[:,vgindex],datalist[:,Idsindex],0)
      bestSS = plt.plot( datalist[indexstart:-1,vgindex], datalist[indexstart:-1,Idsindex]*factorIDS, '-o',lw=1, color='c',markersize=7 )
      target.close()
      count +=1
      print count      
#############


root_dir = '/users/jpduarte/research/variabilityproject/TCAD/data/v7'
data_files = ['finfet_7Wt0.0076Lg0.02Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.86DATAREADY']
#namefile = data_files[0][1][3]#[0][0] give the adress then, [0][1][x] give the file name where x is the string number

#filenames = data_files[0][1]

factorIDS = 1e6/(42e-3*2+7.6e-3)

count = 1
for namefile in data_files:
  if count>20:
		break
  flagfiledata =  namefile.find("DATAREADY") #if it is found "DATAREADY" it returns the position, this can be used for 
  flagfileIon =  namefile.find("vd0.86") 
  if ((flagfiledata > 0) and (flagfileIon > 0)):
    #print namefile
    filenameaux = root_dir +'/'+ namefile
    target = open( filenameaux, 'r')
    fistline = target.readline()
    header = str.split(fistline)
    findvgindex = fistline.find("gateOuterVoltage")     
    if (len(header)>0 and findvgindex>0):
      vgindex =  header.index('gateOuterVoltage')
      Idsindex =  header.index('drainTotalCurrent')
      
      datalist = np.loadtxt(filenameaux,skiprows = 1)
      xx,indexstart = findIds(datalist[:,vgindex],datalist[:,Idsindex],0)
      nominal = plt.plot( datalist[indexstart:-1,vgindex], datalist[indexstart:-1,Idsindex]*factorIDS, '-<',lw=1, color='y',markersize=7, label='Nominal' )
      target.close()
      count +=1
      print count


    
plt.xlabel("VG (V)", fontsize=20)
plt.ylabel("Drain Current (uA/um)", fontsize=20)
ax = plt.gca()
plt.xlim([0,0.86])
plt.xticks(fontsize = 18) 
plt.yticks(fontsize = 18) 
#plt.ylim([0,1000])
#ax.set_yscale('log')  
#plt.savefig('IonsatvsVglin', dpi=600, bbox_inches='tight')  
ax.set_yscale('log')  
plt.xlim([0,0.86])
plt.ylim([1e-2,1e4])
#

##legend
#

plt.legend([tcadsim,nominal,bestSS,worstSS],['TCAD','Nominal','Lowest Ioff','Highest Ioff'],loc=5,prop={'size':20})
#plt.legend(bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure)
plt.savefig('Idsextreme2', dpi=600, bbox_inches='tight')             
plt.show()
