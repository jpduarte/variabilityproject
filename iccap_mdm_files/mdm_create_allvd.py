import numpy
time,Vs,Vsi,Sqf,Sdi,Ise,Ish,Is,Sq,Vg,Vgi,Gqf, Gdi,Ige,Igh,Ig,Gq, Vb, Vbi, Bqf, Bdi,Ibe,Ibh,Ib,Bq, Vd, Vdi, Dqf, Ddi,Ide,Idh,Id,Dq = numpy.loadtxt("/users/sourabh/sourabh/Variability/iccap_mdm_files/finfet_7LgWt0.0076Lg0.02Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.05DATAREADY",skiprows=118,unpack=True)
time2,Vs2,Vsi2,Sqf2,Sdi2,Ise2,Ish2,Is2,Sq2,Vg2,Vgi2,Gqf2, Gdi2,Ige2,Igh2,Ig2,Gq2, Vb2, Vbi2, Bqf2, Bdi2,Ibe2,Ibh2,Ib2,Bq2, Vd2, Vdi2, Dqf2, Ddi2,Ide2,Idh2,Id2,Dq2 = numpy.loadtxt("/users/sourabh/sourabh/Variability/iccap_mdm_files/finfet_7LgWt0.0076Lg0.02Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.43DATAREADY",skiprows=118,unpack=True)
time3,Vs3,Vsi3,Sqf3,Sdi3,Ise3,Ish3,Is3,Sq3,Vg3,Vgi3,Gqf3, Gdi3,Ige3,Igh3,Ig3,Gq3, Vb3, Vbi3, Bqf3, Bdi3,Ibe3,Ibh3,Ib3,Bq3, Vd3, Vdi3, Dqf3, Ddi3,Ide3,Idh3,Id3,Dq3 = numpy.loadtxt("/users/sourabh/sourabh/Variability/iccap_mdm_files/finfet_7LgWt0.0076Lg0.02Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.86DATAREADY",skiprows=118,unpack=True)



f_handle = file('IDVG_LG0.02.mdm','a')

f_handle.write('! VERSION = 6.00 \n')
f_handle.write('BEGIN_HEADER \n')
f_handle.write('ICCAP_INPUTS \n')
f_handle.write('vs         V  S GROUND DEFAULT 0 CON        0 \n')
f_handle.write('vd         V  D GROUND DEFAULT 0 LIST       2 3 0.05 0.43 0.86 \n')
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
f_handle.write('\n')

f_handle.write('BEGIN_DB \n')
f_handle.write('     ICCAP_VAR vs         0               \n')
f_handle.write('     ICCAP_VAR vd         0.43 \n')
f_handle.write('     ICCAP_VAR vb         0 \n')
f_handle.write('\n')
f_handle.write(' #vg              id     \n')
numpy.savetxt(f_handle,numpy.column_stack((Vg2,Id2)),fmt=('%4.4f', '%4.3e'))
f_handle.write('END_DB \n')
f_handle.write('\n')

f_handle.write('BEGIN_DB \n')
f_handle.write('     ICCAP_VAR vs         0               \n')
f_handle.write('     ICCAP_VAR vd         0.86 \n')
f_handle.write('     ICCAP_VAR vb         0 \n')
f_handle.write('\n')
f_handle.write(' #vg              id     \n')
numpy.savetxt(f_handle,numpy.column_stack((Vg3,Id3)),fmt=('%4.4f', '%4.3e'))
f_handle.write('END_DB \n')
f_handle.write('\n')


f_handle.close()
