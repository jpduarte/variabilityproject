import numpy as np
import pylab
import os
import functionsaux as fxsupport
import shutil #this library is to copy file from original to AUX for simulation

def findindexmax(arrayinput,numbertofind, tolerance):
    i=0;
    indexfind = []
    for numberaxu in arrayinput:
        if abs(numbertofind-numberaxu)<tolerance:
           #print "array %d is in range for Vg %f"% (i,numberaxu)
            indexfind.append(i)
        i+=1
    return indexfind

foldername = './v5/' 
nameproject = 'finfet_6'                 
              

#dimensions                
Lg_array    = 0.020*np.array([1, 0.9,1.1])#
Wt_array    = 0.0076*np.array([1, 0.9,1.1])
WbR_array   = 0.0076/2*np.array([1, 0.9,1.1])
WbL_array   = 0.0076/2*np.array([1, 0.9,1.1])
Hfin_array  = 0.042*np.array([1, 0.9,1.1])
tox_array   = 0.0008*np.array([1, 0.9,1.1])
Ro_array    = 0.003*np.array([1])
Doping_FIN_array = 6e18*np.array([1])
Doping_PUNCH_array = [1e21]
Doping_SD_array = [1e21]
#bias
vd_array = ['0.05']
vb_array = ['0.0']
vg_array =['0.86']

Lg, Wt, WbR,WbL, Hfin, tox, Ro, Doping_FIN, Doping_PUNCH, Doping_SD = fxsupport.meshgrid2(Lg_array,Wt_array,WbR_array,WbL_array, Hfin_array, tox_array,Ro_array, Doping_FIN_array, Doping_PUNCH_array,Doping_SD_array )

Ion = []
Ioff = []
Lgplot = []
Vthplot = []

lgindex=0
lpindexplot=0
for L in Lg:
    if lgindex < 1500:#lgindex > 13 and 
        for vg in vg_array:
            for vd in vd_array:
	            for vb in vb_array:
	                print vd
	                filenameaux = foldername + nameproject + 'Wt'+str(Wt[lgindex])+'Lg' + str(L) + 'Hfin'+str(Hfin[lgindex])+'tox'+str(tox[lgindex])+'WbR'+str(WbR[lgindex])+'WbL'+str(WbL[lgindex])+ 'vd' + vd + 'DATAREADY'
	                print "folder name: %s "% (filenameaux) 
	                print "lgindex: %d "% (lgindex)
                    target = open(filenameaux, 'r')
                    header = str.split(target.readline())
                    vgindex =  header.index('gateOuterVoltage')
                    Idsindex =  header.index('drainTotalCurrent')
                    datalist = np.loadtxt(filenameaux,skiprows = 1)
                    dervnumber = 2
                    dIdV = np.diff(datalist[:,Idsindex],n=dervnumber)   
                    indexVth = np.argmax(dIdV)
                    #pylab.plot( datalist[:,vgindex], datalist[:,Idsindex], lw=2 )
                    #pylab.plot( datalist[:-dervnumber,vgindex], dIdV, lw=2 )
                    target.close()
                    print "Vth: %f (V) for Vdd %s (V)"% (datalist[indexVth,vgindex], vd)
                    print "Ids at Vth: %e (V) for Vdd %s (V)"% (datalist[indexVth,Idsindex], vd)
                    print "Ion: %e (A) for Vdd %s (V)"% (np.amax(datalist[:,Idsindex]), vd)
                    print "Ioff: %e (A) for Vdd %s (V) "% (np.amax(datalist[findindexmax(datalist[:,vgindex],0,1e-5),Idsindex]), vd)
                    print "Ion: %f (uA/um) for Vdd %s (V)"% (np.amax(datalist[:,Idsindex])/(0.0076+2*0.042)/1e-6, vd)
                    print "Ioff: %f (nA/um) for Vdd %s (V) "% (np.amax(datalist[findindexmax(datalist[:,vgindex],0,1e-5),Idsindex])/(0.0076+2*0.042)/1e-9, vd)
                    Ion.append(np.amax(datalist[:,Idsindex]))#/(0.0076+2*0.042)/1e-6)
                    Ioff.append(np.amax(datalist[findindexmax(datalist[:,vgindex],0,1e-5),Idsindex]))#/(0.0076+2*0.042)/1e-9)
                    Vthplot.append(datalist[indexVth,vgindex])
                    Lgplot.append(L)
                    
    lgindex +=1

#9.1599991150899303e-09


################################################################plot  
pylab.scatter(Ion,Vthplot,s=80, facecolors='none', edgecolors='k') 
pylab.title("Vd Linear")    
ax = pylab.gca()
pylab.xlabel("Ion (A)", fontsize=18)
pylab.ylabel("Vth", fontsize=18)
#pylab.xlabel("Gate Voltage (V)", fontsize=18)
#pylab.ylabel("Drain Current (A)", fontsize=18)
#ax.set_yscale('log')
#ax.set_xscale('log')
#pylab.ylim([1e-10,2e-4])
#pylab.xlim([0,0.86])
#pylab.savefig('IdVgallvdsatfinfet_6lin.png', dpi=300, bbox_inches='tight')
pylab.savefig('IonVthvlin_fig6.png', dpi=300, bbox_inches='tight')
pylab.show()




