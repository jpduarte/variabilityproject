#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
#this version solve for normalized surface potential 1, not charge

from numpy.matlib import repmat
from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace,minimum,pi
import pylab
import unifiedIMG
##############################################################################################################
def coth(x):
  return 1/tanh(x)

def cot(x):
  return 1/tan(x)

def llog(x):
  return log(abs(x))     

def lsqrt(x):
  return sqrt(abs(x))

def f(phi1,xn,vt,k1,k2,xg1,xg2,A0):
  q1 = -(phi1-xg1)
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)#this is q^2
  qcoth,logsinhqsq,dqcothdq,dqdqsqrt,dlogsinhqsqdqsqrt = logqoversin(qsqrt)
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-logsinhqsq
  equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+qcoth)*(k1*q1+k2*q2)
  return equation4,q2,qsqrt

def logqoversin(qsqrt):
#sin ~ x-x**3/3!
#cos ~ 1 -x**2/2!
  if qsqrt < 0:
    q     = sqrt(-qsqrt) 
    dqdqsqrt = -0.5/q 
    #q = minimum(q,pi/2*0.9999)  
    if(q<1e-4):
      qcoth = (1-qsqrt/8)/(1/2-qsqrt/24) #q*cot(q/2)
      qsqoversinhqsq = 1/((1/2-qsqrt/24))**2
      logsinhqsq = llog(qsqoversinhqsq)
      dqcothdq = (-12/(-12 + qsqrt)**2)*(-q/2.0) #dqcothdqsqrt 
      dlogsinhqsqdqsqrt = 1/(12 *(1/2 - qsqrt/24))
    else:
      sinhqsq = -(sin(q/2))**2
      logsinhqsq = llog(qsqrt/sinhqsq)
      qcoth = q*cot(q/2)
      dqcothdq = (q-sin(q))/(cos(q)-1)
      dlogsinhqsqdqsqrt = 1/qsqrt + cot(q/2)/(q)#TODO: check sign of derivative
  else:
    q     = sqrt(qsqrt)
    dqdqsqrt = 0.5/q     
    if(q<1e-4):
      
      qcoth = (1+qsqrt/8)/(1/2+qsqrt/24) #q*cot(q/2)
      qsqoversinhqsq = 1/((1/2+qsqrt/24))**2
      logsinhqsq = llog(qsqoversinhqsq)
      dqcothdq = (12/(12 + qsqrt)**2)*(q/2.0)
      dlogsinhqsqdqsqrt = 1/(12 *(1/2 + qsqrt/24))      
    else:
      sinhqsq = (sinh(q/2))**2
      logsinhqsq = llog(qsqrt/sinhqsq) 
      qcoth = q*coth(q/2)
      dqcothdq = (-q+sinh(q))/(cosh(q)-1) 
      dlogsinhqsqdqsqrt = 1/qsqrt - cot(q/2)/q#TODO: check sign of derivative           
    
  return qcoth,logsinhqsq,dqcothdq,dqdqsqrt,dlogsinhqsqdqsqrt
    
#############################################################################
  
def fprime(phi1,xn,vt,k1,k2,xg1,xg2,A0):
  q1 = -(phi1-xg1)
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  qcoth,logsinhqsq,dqcothdq,dqdqsqrt,dlogsinhqsqdqsqrt = logqoversin(qsqrt)
  if qsqrt < 0: 
      
    T1 = -qsqrt + k1**2*q1**2#A0*exp(-xn)*exp(-q1 + xg1)
    T0 = -T1/2 - k1**2*q1#-2(d(qsqrt)/dq1) = -A0*exp(-xn)*exp(-q1 + xg1)-k1**2*q1
    q = lsqrt(-qsqrt)
    #q = minimum(q,pi/2*0.9999) 
    
    sinaux = sin(q/2)
    q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-logsinhqsq 
       
    term1 = T1 + (k1 + k2*(1 + ((T0)*qcoth - (T1 + 2*k1**2*q1))/(qsqrt) + 2*(k1 + (T0)*(-cot(q/2)**2 - 1)/2 + (T0)*cot(q/2)/q)/(k1*q1 + qcoth)))*(k1*q1 + qcoth) + (k1*q1 + k2*q2)*(k1 + (T0)*(-cot(q/2)**2 - 1)/2 + (T0)*cot(q/2)/q)
  else:
    T1 = -qsqrt + k1**2*q1**2#A0*exp(-xn)*exp(-q1 + xg1)
    T0 = T1/2 + k1**2*q1#-2(d(qsqrt)/dq1) = A0*exp(-xn)*exp(-q1 + xg1)/2 +  k1**2*q1
    q = sqrt(qsqrt)
    q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-logsinhqsq
    
    term1 = T1 + (k1 + k2*(1 + ((T0)*qcoth - (2*T0))/(qsqrt) + 2*(k1 - (T0)/(2*sinh(q/2)**2) + (T0)*coth(q/2)/q)/(k1*q1 + q*coth(q/2))))*(k1*q1 + qcoth) + (k1*q1 + k2*q2)*(k1 - (T0)/(2*sinh(q/2)**2) + (T0)*coth(q/2)/q)  
  return -term1

 
vt  = 0.026
xn  = 0
tox = 1e-9 #oxide thickness
tbox = 1e-9 #box thickness
tsi = 10e-9 #silicon film thickness
es  = 11.9*8.8542e-12 #silicon permittivity 
q   = 1.6e-19
ni  = 1.45e16
eox = 3.9*8.8542e-12
cox = eox/tox
cbox = eox/tbox
csi = es/tsi
k1  = cox/csi
k2  = cbox/csi
A0  = (2*q*ni*tsi*tsi)/(es*vt)

Eg 	    = 1.169640 
Nc      = 2.540000e+19 
Nv 	    = 3.140000e+19 

Vgb_array = linspace(0, 1, 100)
Vgf_array = linspace(-1, 1, 6) #front gate fixed in this case

xg1 = Vgb_array[0] / vt
xg2 = Vgf_array[0] / vt

#####################################Newton Method####################################
xg1 = Vgb_array[0] / vt
xg2 = -6 / vt

iter=1
for Vgf in Vgf_array:
  xg1 = Vgb_array[0] / vt
  xg2 = Vgf / vt
  phinew = -9.0#0.17/vt 
  q1  = []
  q2  = []
  qtotal = []
  phi1 = []
  qsqrt = []
  itertotal = []
  fprimeall = []
  #print "solving for Vgf: %f, with initial guess %f" % (Vgf,q1guess)
  for Vgb in Vgb_array:
    #print "solving: " + str(iter)
    iter+=1
    xg1 = Vgb / vt
    xg2 = Vgf / vt
    data = (xn,vt,k1,k2,xg1,xg2,A0)
    
    iternnew = 1
    phiold = 1000
    while iternnew<100:
      deltaphis = 1e-4
      fval,q2aux,qsqrtaux = f(phinew,xn,vt,k1,k2,xg1,xg2,A0)
      #fval2 = f(phinew+deltaphis,xn,vt,k1,k2,xg1,xg2,A0)
      if abs(fval) < 1e-12:

        break
      fprimeval = fprime(phinew,xn,vt,k1,k2,xg1,xg2,A0)
      phinew = phinew - fval/fprimeval

      iternnew+=1
      
    #qfguessaux,qbguessaux=  unifiedIMG.unified_IMG_model(Vgf,Vgb,0,q,k,T,eox/3.9,eox,es,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,phi_gateb,Tsi,tins,tinsbox,Hfin)  
    fprimeall.append(fprimeval)
    phi1.append(phinew)
    q1aux = (xg1-phinew)
    #q2aux, qsqrtaux = q2value(q1aux,xn,vt,k1,k2,xg1,xg2,A0)
    qsqrt.append(qsqrtaux)
    q1.append(q1aux*vt*cox)
    q2.append(q2aux*vt*cbox)
    qtotal.append((q1aux*vt*cox+q2aux*vt*cbox) )
    itertotal.append(iternnew)
  pylab.figure(2)
  pylab.plot( Vgb_array,fprimeall,'o')
  pylab.savefig('fprimefig.png', dpi=300, bbox_inches='tight') 
  pylab.figure(3)
  pylab.plot( Vgb_array,q1,'o')
  pylab.savefig('q1fig.png', dpi=300, bbox_inches='tight')
  pylab.figure(4)
  pylab.plot( Vgb_array,q2,'o')
  pylab.savefig('q2fig.png', dpi=300, bbox_inches='tight')
  pylab.figure(5)
  pylab.plot( Vgb_array,qtotal,'o')
  pylab.savefig('qtlinfig.png', dpi=300, bbox_inches='tight')
  pylab.figure(6)
  pylab.plot( Vgb_array,qtotal,'o')
  ax = pylab.gca() 
  ax.set_yscale('log')
  pylab.savefig('qtlogfig.png', dpi=300, bbox_inches='tight') 
  pylab.figure(7)
  pylab.plot( Vgb_array,phi1,'o')
  pylab.savefig('phi1fig.png', dpi=300, bbox_inches='tight')  
  pylab.figure(8)
  pylab.plot( Vgb_array,qsqrt,'o') 
  pylab.savefig('qsqrtfig.png', dpi=300, bbox_inches='tight')   
  pylab.figure(9)
  pylab.plot( Vgb_array,itertotal,'o')   
  pylab.savefig('iterationtotalfig.png', dpi=300, bbox_inches='tight')    
pylab.show()



