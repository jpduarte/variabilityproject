import numpy as np
import pylab
from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots
from numpy import real

def unified_charge_model(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch) :
  rc  = (2*Cins/(Weff**2*ech/Ach))
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2-vt*np.log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+np.log(Cins*vt/(q*ni**2*2*Ach/Nch))
  vth_fixed_factor_Sub = np.log((qdep*rc)**2/(np.exp(qdep*rc)-qdep*rc-1))+vth_fixed_factor_SI
  vth_N_Sub = -qdep+vth_fixed_factor_Sub
  vth_N_SI  = -qdep+vth_fixed_factor_SI
  Vg_local_N  = (Vg-Vch)*vt**(-1)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5
  if (Vov>60):
    qm      = -Vov*2
    qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
    x0      = qtrc*(np.exp(qtrc)-qtrc-1)**(-1)
    x1      = qtrc*x0
    f0      = F-qm+np.log(-qm)+np.log(x1)
    f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
    f2      = -(qm**2)**(-1)
    qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
  else:
    qm      = np.exp((Vg_local_N-vth_N_Sub)*0.5)
    if(qm>1e-7):
      qm      = 2*(1-np.sqrt(1+(np.log(1+qm))**2))
      qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
      x0      = qtrc*(np.exp(qtrc)-qtrc-1)**(-1)
      x1      = qtrc*x0
      f0      = F-qm+np.log(-qm)+np.log(x1)
      f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
      f2      = -(qm**2)**(-1)
      qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
    else:
      qm      = -qm**2
  return qm

#Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch
Vg=1
Vch=0

q =		1.602190e-19 
k =		1.380650e-23 
T =		300 
eo =		8.854000e-14 
eins =		3.453060e-13 
ech =		1.035918e-12 
Eg 	=	1.169640 
Nc =		2.540000e+19 
Nv 	=	3.140000e+19 
vt = 0.0259
ni = 1e10

phi_substrate = 4.05 
phi_gate 	= 4.3
Hfin = 10e-9
Tfin = 5e-9
tins = 1e-9

alpha_MI = 1.6
Cins = (2*Hfin)*eins/tins
Ach=Hfin*Tfin
Weff = 2*Tfin
Nch = 1e15
################################# qm versus Vg
Vg_array = np.linspace(0,1,50)
qm = []
for Vg in Vg_array:
  qm.append(-unified_charge_model(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch))

pylab.figure(1)
pylab.plot( Vg_array, qm, lw=2 )
#pylab.show()
################################# Ids versus Vg
def fixed_quad(func,a,b,n,args=()) :
  [x,w] = p_roots(n)
  x = real(x)
  y = (b-a)*(x+1)/2.0 + a
  sum_integral = 0
  i=0
  for xi in x:
    sum_integral = sum_integral+w[i]*func((b-a)*(xi+1)/2.0 + a,*args)
    i=i+1
  return (b-a)/2.0*sum_integral


def current_integral(qm, a, b, c,eps):
     return -qm*(1-1/qm)/(1-a*qm+b*qm**2-c*(1/(eps+qm)))

def unified_current_model(qs,qd,a,b,c,eps):
  Is = qs**2/2-qs
  Id = qd**2/2-qd
  qave=(qs+qd)/2
  return (Is-Id)/(1-a*qave+b*qave**2-c*(1/(eps+qave)))
  
Vd_array = np.linspace(1,1,1)#np.array([0.05,0.5, 1])
Vs = 0
Ids = []
Idsnum = []
Idsgauss = []
error = []
for Vg in Vg_array:
  for Vd in Vd_array:
    qs = unified_charge_model(Vg,Vs,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch)
    qd = unified_charge_model(Vg,Vd,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch)
    Ids.append(unified_current_model(qs,qd,0.1,0.001,1,-0.1))
    Idsaux = quad(current_integral, qs, qd, args=(0.1,0.001,1,-0.1))
    Idsnum.append(Idsaux[0])
    Idsaux2 = fixed_quad(current_integral, qs, qd, 3,[0.1,0.001,1,-0.1]) #quadrature(current_integral, qs, qd, args=(Vg,Vd), maxiter=3)
    Idsgauss.append(Idsaux2)
    error.append(100*np.abs((Idsaux[0]-Idsaux2)/Idsaux[0]))

  
pylab.figure(2)
pylab.plot( Vg_array, np.reshape(Ids, (len(Vg_array),len(Ids)/len(Vg_array))), lw=2 )
pylab.plot( Vg_array, np.reshape(Idsgauss, (len(Vg_array),len(Idsgauss)/len(Vg_array))),'bo' )
pylab.xlabel("Gate Voltage (V)", fontsize=18)
pylab.ylabel("Drain Current (A)", fontsize=18)
#pylab.savefig('Idsmobility2linoldvsnew.png', dpi=300, bbox_inches='tight')

pylab.figure(3)
pylab.plot( Vg_array, np.reshape(Ids, (len(Vg_array),len(Ids)/len(Vg_array))), lw=2 )
pylab.plot( Vg_array, np.reshape(Idsgauss, (len(Vg_array),len(Idsgauss)/len(Vg_array))),'bo' )
pylab.xlabel("Gate Voltage (V)", fontsize=18)
pylab.ylabel("Drain Current (A)", fontsize=18)
#pylab.savefig('Idsmobility2logoldvsnew.png', dpi=300, bbox_inches='tight')
ax = pylab.gca()
ax.set_yscale('log')



pylab.figure(4)
pylab.plot( Vg_array, np.reshape(Idsnum, (len(Vg_array),len(Idsnum)/len(Vg_array))) )
pylab.plot( Vg_array, np.reshape(Idsgauss, (len(Vg_array),len(Idsgauss)/len(Vg_array))),'bo' )
ax = pylab.gca()
pylab.xlabel("Gate Voltage (V)", fontsize=18)
pylab.ylabel("Drain Current (A)", fontsize=18)
#pylab.savefig('Idsmobility2lin.png', dpi=300, bbox_inches='tight')
#ax.set_yscale('log')

pylab.figure(5)
pylab.plot( Vg_array, np.reshape(Idsnum, (len(Vg_array),len(Idsnum)/len(Vg_array))) )
pylab.plot( Vg_array, np.reshape(Idsgauss, (len(Vg_array),len(Idsgauss)/len(Vg_array))),'bo' )
ax = pylab.gca()
pylab.xlabel("Gate Voltage (V)", fontsize=18)
pylab.ylabel("Drain Current (A)", fontsize=18)
ax.set_yscale('log')
#pylab.savefig('Idsmobility2log.png', dpi=300, bbox_inches='tight')

pylab.figure(6)
pylab.plot( Vg_array, np.reshape(error, (len(Vg_array),len(error)/len(Vg_array))) )
pylab.xlabel("Gate Voltage (V)", fontsize=18)
pylab.ylabel("Error %", fontsize=18)
#pylab.savefig('ErrorIdsmobility2.png', dpi=300, bbox_inches='tight')

dIdV = np.diff(Ids,n=1)
dIdVgauss = np.diff(Idsgauss,n=1)

pylab.figure(7)
pylab.plot( Vg_array[:-1], dIdV,lw=2  )
pylab.plot( Vg_array[:-1], dIdVgauss,'bo'  )


pylab.show()
