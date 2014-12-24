#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
#this version solve for normalized surface potential 1, not charge
#import numpy as np
from numpy.matlib import repmat
from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace
import pylab
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
  """q     = sqrt(abs(qsqrt))
  qcoth = q*cot(q/2)
  sinhqsq = (sin(q/2))**2"""
  """if qsqrt < 0:
    q     = sqrt(-qsqrt)
    qcoth = q*cot(q/2)
    sinhqsq = -(sin(q/2))**2
  else:
    q     = sqrt(qsqrt)
    qcoth = q*coth(q/2) 
    sinhqsq = (sinh(q/2))**2"""
    
  qcoth,logsinhqsq = logqoversin(qsqrt)
     
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-logsinhqsq
  equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+qcoth)*(k1*q1+k2*q2)
  return equation4

def logqoversin(qsqrt):
  if qsqrt < 0:
    q     = sqrt(-qsqrt)
    qcoth = q*cot(q/2)
    if(q<1e-4):
      qsqoversinhqsq = 1/((1/2-qsqrt/6))**2
      logsinhqsq = llog(qsqoversinhqsq)
    else:
      sinhqsq = -(sin(q/2))**2
      logsinhqsq = llog(qsqrt/sinhqsq)
  else:
    q     = sqrt(qsqrt)
    qcoth = q*coth(q/2)
    if(q<1e-4):
      qsqoversinhqsq = 1/((1/2+qsqrt/6))**2
      logsinhqsq = llog(qsqoversinhqsq)
    else:
      sinhqsq = (sinh(q/2))**2
      logsinhqsq = llog(qsqrt/sinhqsq) 
    
  return qcoth,logsinhqsq
    
def g(phi1,xn,vt,k1,k2,xg1,xg2,A0):
  q1 = -(phi1-xg1)
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)#this is q^2
  
  if qsqrt < 0:
    q     = sqrt(-qsqrt)
    qcoth = q*cot(q/2)
    sinhqsq = -(sin(q/2))**2
  else:
    q     = sqrt(qsqrt)
    qcoth = q*coth(q/2) 
    sinhqsq = (sinh(q/2))**2
     
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-llog(qsqrt/sinhqsq)
  equation4 = (-A0*exp(-xn+xg1-q1)/(k1*q1+k2*q2)-qcoth)/k1+xg1
  return equation4  

def q2value(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  if qsqrt < 0:
    q     = sqrt(-qsqrt)
    qcoth = q*cot(q/2)
    sinhqsq = -(sin(q/2))**2
  else:
    q     = sqrt(qsqrt)
    qcoth = q*coth(q/2) 
    sinhqsq = (sinh(q/2))**2
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-llog(qsqrt/sinhqsq)
  return q2,qsqrt  
  
def fprime(phi1,xn,vt,k1,k2,xg1,xg2,A0):
  q1 = -(phi1-xg1)
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  if qsqrt < 0:
    term1 = A0*exp(-q1 + xg1 - xn) + (k1 + k2*(1 + ((-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*(-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*cos(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)/(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*sin(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**3) - (A0*exp(-xn)*exp(-q1 + xg1) + 2*k1**2*q1)/sin(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2)*sin(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2/(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2) + 2*(k1 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*(-cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2 - 1)/2 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)/lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2))/(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2))))*(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)) + (k1*q1 + k2*(q1 - xg1 + xg2 - llog(-(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/sin(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2) + 2*llog(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2))))*(k1 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*(-cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2 - 1)/2 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)/lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2))
  else:
    term1 = A0*exp(-q1 + xg1 - xn) + (k1 + k2*(1 - (-lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*(A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)*cosh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)/sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**3 + (A0*exp(-xn)*exp(-q1 + xg1) + 2*k1**2*q1)/sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2)*sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2/(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2) + 2*(k1 - (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)/(2*sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2) + (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)/lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2))/(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2))))*(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)) + (k1*q1 + k2*(q1 - xg1 + xg2 - llog((-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2) + 2*llog(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2))))*(k1 - (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)/(2*sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2) + (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)/lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2))  
  return -term1

 
vt  = 0.026
xn  = 0
tox = 1e-9 #oxide thickness
tsi = 7e-9 #silicon film thickness
tbox = 25e-9 #box thickness
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

Vgb_array = linspace(1, 1, 100)
Vgf_array = linspace(6, 6, 6) #front gate fixed in this case

xg1 = Vgb_array[0] / vt
xg2 = Vgf_array[0] / vt

#####################################Newton Method####################################
xg1 = Vgb_array[0] / vt
xg2 = -6 / vt
phi1guess = linspace(20, 25, 1000)

pylab.figure(1)
fval = []
for phi in phi1guess:

  fval = f(phi,xn,vt,k1,k2,xg1,xg2,A0)
  fprimeval = fprime(phi,xn,vt,k1,k2,xg1,xg2,A0)
  if fval>0:
    pylab.figure(1)
    pylab.plot( phi, fval,'ob')
  else:
    pylab.figure(1)
    pylab.plot( phi, fval,'or')
    
  if fprimeval>0:
    pylab.figure(2)
    pylab.plot( phi, fprimeval,'ob')
  else:
    pylab.figure(2)
    pylab.plot( phi, fprimeval,'or')
"""
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
  #print "solving for Vgf: %f, with initial guess %f" % (Vgf,q1guess)
  for Vgb in Vgb_array:
    #print "solving: " + str(iter)
    iter+=1
    xg1 = Vgb / vt
    xg2 = Vgf / vt
    data = (xn,vt,k1,k2,xg1,xg2,A0)
    
    iternnew = 1
    phiold = 1000
    while iternnew<50:
      deltaphis = 1e-4
      fval = f(phinew,xn,vt,k1,k2,xg1,xg2,A0)
      #fval2 = f(phi1guess+deltaphis,xn,vt,k1,k2,xg1,xg2,A0)
      if abs(fval) < 1e-12:
        #print "solved: %d with %d iterations" % (iter,iternnew)
        break
      fprimeval = fprime(phinew,xn,vt,k1,k2,xg1,xg2,A0)
      phinew = phinew - fval/fprimeval
      #phi1guess = phi1guess - fval/((fval2-fval)/deltaphis)
      #phiold = phi1guess
      #gvalue = g(phi1guess,xn,vt,k1,k2,xg1,xg2,A0)
      #phi1guess = gvalue
       
      #print iternnew,phi1guess,gvalue,abs(phi1guess-phiold),fprimeval,fval
      iternnew+=1
    phi1.append(phinew)
    q1aux = (xg1-phinew)
    q2aux, qsqrtaux = q2value(q1aux,xn,vt,k1,k2,xg1,xg2,A0)
    qsqrt.append(qsqrtaux)
    q1.append(q1aux*vt*cox)
    q2.append(q2aux*vt*cox)
    qtotal.append((q1aux+q2aux)*vt*cox )
    pylab.figure(3)
      pylab.plot( iternnew,fval,'o')
      ax = pylab.gca()      
      ax.set_yscale('log') 
      pylab.figure(4)
      pylab.plot( iternnew,fprimeval,'o')
      pylab.figure(5)
      pylab.plot( iternnew,phinew,'o')
  pylab.figure(3)
  pylab.plot( Vgb_array,q1,'o')
  pylab.figure(4)
  pylab.plot( Vgb_array,q2,'o')
  pylab.figure(5)
  pylab.plot( Vgb_array,qtotal,'o')
  pylab.figure(6)
  pylab.plot( Vgb_array,phi1,'o')  
  pylab.figure(7)
  pylab.plot( Vgb_array,qsqrt,'o')   """   
pylab.show()



