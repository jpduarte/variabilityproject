import numpy
time,Vs,Vsi,Sqf,Sdi,Ise,Ish,Is,Sq,Vg,Vgi,Gqf, Gdi,Ige,Igh,Ig,Gq, Vb, Vbi, Bqf, Bdi,Ibe,Ibh,Ib,Bq, Vd, Vdi, Dqf, Ddi,Ide,Idh,Id,Dq = numpy.loadtxt("/users/sourabh/sourabh/FinFET_Variability/variabilityproject/netlist_modelcards/Wt0.00836Lg0.018Hfin0.0462tox0.00088WbR0.00418WbL0.00418Nfin5.4e+18v0.86",skiprows=118,unpack=True)

f_handle = file('idvgsat_worst.dat','a')

f_handle.write('! VERSION = 6.00 \n')
f_handle.write('BEGIN_HEADER \n')
f_handle.write('ICCAP_INPUTS \n')
f_handle.write('vs         V  S GROUND DEFAULT 0 CON        0 \n')
f_handle.write('vd         V  D GROUND DEFAULT 0 CON        0.05 \n')
f_handle.write('vb         V  B GROUND DEFAULT 0 CON        0 \n')
f_handle.write('vg         V  G GROUND V 0 LIN        1    0          0.86       100  0.00868687 \n')
f_handle.write('ICCAP_OUTPUTS \n')
f_handle.write('id         I  D GROUND DEFAULT B \n')
f_handle.write('END_HEADER \n')

f_handle.write('\n')

f_handle.write('BEGIN_DB \n')
f_handle.write('     ICCAP_VAR vs         0               \n')
f_handle.write('     ICCAP_VAR vd         0.05 \n')
f_handle.write('     ICCAP_VAR vb         0 \n')

f_handle.write('\n')

f_handle.write(' #vg              id     \n')
numpy.savetxt(f_handle,numpy.column_stack((Vg,Id)),fmt=('%4.4f', '%4.3e'))
f_handle.write('END_DB \n')

f_handle.close()
