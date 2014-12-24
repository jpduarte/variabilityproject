#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
#this version solve for charge
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
  return log(abs(x)+1e-10)     

def lsqrt(x):
  return sqrt(abs(x))

def f(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)#this is q^2
  q     = lsqrt(qsqrt)
  if qsqrt < 0:
    qcoth = q*cot(q/2)
  else:
    qcoth = q*coth(q/2)  
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-llog(qsqrt/(sinh(qsqrt/4)))
  equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+qcoth)*(k1*q1+k2*q2)
  print qsqrt,q,qcoth,q2,equation4
  return equation4

def q2value(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  q     = lsqrt(qsqrt)
  if qsqrt < 0:
    qcoth = q*cot(q/2)
  else:
    qcoth = q*coth(q/2)  
  q2    = xg2-xg1+q1+2*llog(k1*q1+qcoth)-llog(qsqrt/(sinh(qsqrt/4)))
  return q2,qsqrt  
  
def fprime(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  if qsqrt < 0:
    term1 = A0*exp(-q1 + xg1 - xn) + (k1 + k2*(1 - (-(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*(A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1/2)*cosh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)**2 + (A0*exp(-xn)*exp(-q1 + xg1) + 2*k1**2*q1)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4))*sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)/(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2) + 2*(k1 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*(-cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2 - 1)/2 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)/lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2))/(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2))))*(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)) + (k1*q1 + k2*(q1 - xg1 + xg2 - llog((-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)) + 2*llog(k1*q1 + lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2))))*(k1 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*(-cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)**2 - 1)/2 + (-A0*exp(-xn)*exp(-q1 + xg1)/2 - k1**2*q1)*cot(lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2)/2)/lsqrt(A0*exp(-xn)*exp(-q1 + xg1) - k1**2*q1**2))
  else:
    term1 = A0*exp(-q1 + xg1 - xn) + (k1 + k2*(1 - (-(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*(A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1/2)*cosh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)**2 + (A0*exp(-xn)*exp(-q1 + xg1) + 2*k1**2*q1)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4))*sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)/(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2) + 2*(k1 - (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)/(2*sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2) + (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)/lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2))/(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2))))*(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)) + (k1*q1 + k2*(q1 - xg1 + xg2 - llog((-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/sinh(-A0*exp(-xn)*exp(-q1 + xg1)/4 + k1**2*q1**2/4)) + 2*llog(k1*q1 + lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2))))*(k1 - (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)/(2*sinh(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)**2) + (A0*exp(-xn)*exp(-q1 + xg1)/2 + k1**2*q1)*coth(lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2)/2)/lsqrt(-A0*exp(-xn)*exp(-q1 + xg1) + k1**2*q1**2))    
  return term1

 
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

Vgb_array = linspace(2, 0, 2)
Vgf_array = linspace(-6, -6, 1) #front gate fixed in this case

xg1 = Vgb_array[0] / vt
xg2 = Vgf_array[0] / vt

#this shows that f is not smooth as a function of q1, check figure 6?
q1try = linspace(50, 54, 100)
for q1tryaux in q1try:
  pylab.figure(6)
  pylab.plot( q1tryaux,f(q1tryaux,xn,vt,k1,k2,xg1,xg2,A0),'o')  
pylab.show()
#####################################Newton Method####################################
#this does not converge
iter=1
for Vgf in Vgf_array:
  xg1 = Vgb_array[0] / vt
  xg2 = Vgf / vt
  q1guess = 51 
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
    while iternnew<10:
      fval = f(q1guess,xn,vt,k1,k2,xg1,xg2,A0)
      if abs(fval) < 1e-12:
        #print "solved: %d with %d iterations" % (iter,iternnew)
        break
      q1guess = q1guess - fval/fprime(q1guess,xn,vt,k1,k2,xg1,xg2,A0)
      print q1guess,f(q1guess,xn,vt,k1,k2,xg1,xg2,A0),iternnew
      iternnew+=1
    #q1guess = fsolve(f, q1guess, args=data,xtol=1.49012e-14)
    #print "solved for Vgb: %f Vbf: %f with %d iterations, q1: %f" % (Vgb,Vgf,iternnew,q1guess*vt)
    
    phi1.append(-q1guess+xg1)
    q1aux = q1guess #(xg1-q1guess)
    q2aux, qsqrtaux = q2value(q1aux,xn,vt,k1,k2,xg1,xg2,A0)
    qsqrt.append(qsqrtaux)
    q1.append(q1aux*vt*cox)
    q2.append(q2aux*vt*cox)
    qtotal.append((q1aux+q2aux)*vt*cox )
  pylab.figure(1)
  pylab.plot( Vgb_array,q1,'o')
  pylab.figure(2)
  pylab.plot( Vgb_array,q2,'o')
  pylab.figure(3)
  pylab.plot( Vgb_array,qtotal,'o')
  pylab.figure(4)
  pylab.plot( Vgb_array,phi1,'o')  
  pylab.figure(5)
  pylab.plot( Vgb_array,qsqrt,'o')  
pylab.show()



