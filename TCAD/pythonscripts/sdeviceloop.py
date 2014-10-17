#run loops for sdevice simulation
import os
import functionsaux as fxsupport
import shutil #this library is to copy file from original to AUX for simulation
import numpy as np
                
foldername = './v6/' 
nameproject = 'finfet_6'               

#dimensions                
Lg_array    = 0.020*np.array([1, 0.9,1.1])#
Wt_array    = 0.0076*np.array([1, 0.9,1.1])
WbR_array   = 0.0076/2*np.array([1, 0.9,1.1])
WbL_array   = 0.0076/2*np.array([1, 0.9,1.1])
Hfin_array  = 0.042*np.array([1, 0.9,1.1])
tox_array   = 0.0008*np.array([1, 0.9,1.1])
Ro_array    = 0.003*np.array([1])
Doping_FIN_array = 6e18*np.array([1, 0.9,1.1])
Doping_PUNCH_array = [1e21]
Doping_SD_array = [1e21]
#bias
vd_array = ['0.86','0.05','0.43']
vb_array = ['0.0']
vg_array =['0.86']

Lg, Wt, WbR,WbL, Hfin, tox, Ro, Doping_FIN, Doping_PUNCH, Doping_SD = fxsupport.meshgrid2(Lg_array,Wt_array,WbR_array,WbL_array, Hfin_array, tox_array,Ro_array, Doping_FIN_array, Doping_PUNCH_array,Doping_SD_array )
#
#
lgindex=0
for L in Lg:

    filename = foldername + 'sde_dvs.cmd'
    filenameaux = foldername + 'sde_dvsAUX.cmd'
    shutil.copyfile(filename,filenameaux)

    #example of update values in aux file
    if lgindex<5000:#lgindex>43 and lgindex<47 :
        fxsupport.inplace_change(filenameaux,'@Lg_user@',str(L))
        fxsupport.inplace_change(filenameaux,'@Wfintop_user@',str(Wt[lgindex]))
        fxsupport.inplace_change(filenameaux,'@WsiBL_user@',str(WbL[lgindex]))
        fxsupport.inplace_change(filenameaux,'@WsiBR_user@',str(WbR[lgindex]))
        fxsupport.inplace_change(filenameaux,'@Hfin_user@',str(Hfin[lgindex]))
        fxsupport.inplace_change(filenameaux,'@tox_user@',str(tox[lgindex]))
        fxsupport.inplace_change(filenameaux,'@Ro_user@',str(Ro[lgindex]))
        fxsupport.inplace_change(filenameaux,'@Doping_FIN@',str(Doping_FIN[lgindex]))
        fxsupport.inplace_change(filenameaux,'@Doping_PUNCH@',str(Doping_PUNCH[lgindex]))
        fxsupport.inplace_change(filenameaux,'@DOPING_SD@',str(Doping_SD[lgindex]))

        stringmeshname = foldername + nameproject +'_'+str(lgindex)
        fxsupport.inplace_change(filenameaux,'n@node@_msh',stringmeshname)#TODO: no need ot msh in the name, it automatically will add msh
        print "Generating structure for simulation"
        if lgindex<1500:#lgindex>44 and lgindex<46 :
            print "lgindex: ", lgindex
            os.system('sde -e -l ' + foldername +'sde_dvsAUX.cmd')
            
            for vg in vg_array:
                for vd in vd_array:
	                for vb in vb_array:
	                    filename = foldername +'sdevice_des.cmd'
	                    filenameaux = foldername +'sdevice_desAUX.cmd'
	                    shutil.copyfile(filename,filenameaux)
	                    
	                    fxsupport.inplace_change(filenameaux,'@parameter@',foldername+'models.par')
	                    fxsupport.inplace_change(filenameaux,'@tdr@',stringmeshname+'_msh.tdr')
	                    stringplot = foldername + nameproject + 'Wt'+str(Wt[lgindex])+'Lg' + str(L) + 'Hfin'+str(Hfin[lgindex])+'tox'+str(tox[lgindex])+'WbR'+str(WbR[lgindex])+'WbL'+str(WbL[lgindex])+ 'vd' + vd 
	                    fxsupport.inplace_change(filenameaux,'@plot@',stringplot)
	                    fxsupport.inplace_change(filenameaux,'@WF_gate@','4.235')
	                    fxsupport.inplace_change(filenameaux,'@Vd_user@',vd)
	                    fxsupport.inplace_change(filenameaux,'@Vb_user@',vb)
	                    fxsupport.inplace_change(filenameaux,'@Vginitial_user@','0.00')
	                    fxsupport.inplace_change(filenameaux,'@Vgfinal_user@',vg)
	                    print "Running device simulation"
	                    os.system('sdevice ' + filenameaux)
	                    os.system('python ./pythonscripts/readcsdeviceplt.py createdataout '+ stringplot +' '+stringplot+'DATAREADY')
    lgindex +=1
    
    
#print Lg
