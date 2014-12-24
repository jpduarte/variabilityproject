#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
#this version solve for surface potential 1, not charge
#import numpy as np
from numpy.matlib import repmat
from numpy import exp,tan,tanh,log,cosh,sqrt,sinh,sin,cos,linspace
import pylab
#from scipy.optimize import fsolve
#from sympy import exp,tan,tanh,log,cosh,sqrt,sinh,sin,cos

##############################################################################################################
def coth(x):
  return 1/tanh(x)

def f(q1,xn,vt,k1,k2,xg1,xg2,A0):
  q1aux = -(q1-xg1)
  qsqrt = k1**2*q1aux**2 - A0*exp(-q1aux+xg1- xn)
  q = sqrt(qsqrt)
  q2    = q1aux-xg1+xg2+2.0*log(coth(q/2.0)*q+k1*q1aux)-log((qsqrt)/sinh(qsqrt/4.0))
  term1 = (k1*q1aux+coth(q/2.0)*q)*(k1*q1aux+k2*q2)-A0*exp(-q1aux+xg1-xn)
  
  q1 = -(q1-xg1)
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  q     = sqrt(qsqrt)
  q2    = xg2-xg1+q1+2*log(k1*q1+q*coth(q/2))-log(qsqrt/(sinh(qsqrt/4)))
  equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+q*coth(q/2))*(k1*q1+k2*q2)#
  return term1

def q2value(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  q     = sqrt(qsqrt)
  q2    = xg2-xg1+q1+2*log(k1*q1+q*coth(q/2))-log(qsqrt/(sinh(qsqrt/4)))
  return q2  
  
 
"""
def f(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  if qsqrt >= -1e20:
    q = sqrt(qsqrt)
    q2  = xg2-xg1+q1+2*log(k1*q1+q*1/tanh(q/2))-log(qsqrt/(sinh(q/2))**2)
    equation4 = -A0*exp(-xn)*exp(xg1-q1)+(k1*q1+q*1/tanh(q/2))*(k1*q1+k2*q2)#
  else:
    q = sqrt(-qsqrt)   
    q2  = xg2-xg1+q1+2*log(k1*q1+q/tan(-q/2))-log(-qsqrt/(sin(-q/2))**2)
    equation4 = -A0*exp(-xn)*exp(xg1-q1)+(k1*q1-q*1/tan(-q/2))*(k1*q1+k2*q2)
  return equation4"""

 
def fprime(q1,xn,vt,k1,k2,xg1,xg2,A0):
  #x = complex(q1,0)
  x = q1
  term1 = -A0*exp(x - xn) + (k1 + k2*(1 - 2.0*(-k1 + (-A0*exp(x - xn) + k1**2*(x - xg1)**2)**(-0.5)*(-0.5*A0*exp(x - xn) + 0.5*k1**2*(2*x - 2*xg1))/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5) - 0.5*(-0.5*A0*exp(x - xn) + 0.5*k1**2*(2*x - 2*xg1))*(-tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)**2 + 1)/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)**2)/(-k1*(x - xg1) + (-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)) + (-(-A0*exp(x - xn) + k1**2*(x - xg1)**2)*(-0.25*A0*exp(x - xn) + 0.25*k1**2*(2*x - 2*xg1))*cosh(-0.25*A0*exp(x - xn) + 0.25*k1**2*(x - xg1)**2)/sinh(-0.25*A0*exp(x - xn) + 0.25*k1**2*(x - xg1)**2)**2 + (-A0*exp(x - xn) + k1**2*(2*x - 2*xg1))/sinh(-0.25*A0*exp(x - xn) + 0.25*k1**2*(x - xg1)**2))*sinh(-0.25*A0*exp(x - xn) + 0.25*k1**2*(x - xg1)**2)/(-A0*exp(x - xn) + k1**2*(x - xg1)**2)))*(k1*(x - xg1) - (-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)) + (k1*(x - xg1) + k2*(x - xg2 + log((-A0*exp(x - xn) + k1**2*(x - xg1)**2)/sinh(-0.25*A0*exp(x - xn) + 0.25*k1**2*(x - xg1)**2)) - 2.0*log(-k1*(x - xg1) + (-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5))))*(k1 - (-A0*exp(x - xn) + k1**2*(x - xg1)**2)**(-0.5)*(-0.5*A0*exp(x - xn) + 0.5*k1**2*(2*x - 2*xg1))/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5) + 0.5*(-0.5*A0*exp(x - xn) + 0.5*k1**2*(2*x - 2*xg1))*(-tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)**2 + 1)/tanh(0.5*(-A0*exp(x - xn) + k1**2*(x - xg1)**2)**0.5)**2)
  return term1

def q1guessf(xn,vt,k1,k2,xg1,xg2,A0):
  A1  = k1**2+k1*k2+k2*k1**2;
  B1  = 2*k1+2*k2+2*k1*k2+k1*k2*(xg2-xg1)
  C1  = 2*k2*(xg2-xg1)
  D1  = B1/A1
  E1  = C1/A1
  F1  = D1+2*xg1
  G1  = xg1*D1+E1+xg1*xg1
  h1  = (F1/2)**2-G1
  xsub1 = (F1/2)-sqrt(h1) 
  neta1 = 0.5*(xsub1+xn+13-sqrt(64+(xsub1-xn-13)**2))
  return 2*neta1 #initial guess  
  
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

Vgb_array = linspace(0, 2, 100)
Vgf_array = linspace(-6, -6, 1) #front gate fixed in this case

xg1 = Vgb_array[0] / vt
xg2 = Vgf_array[0] / vt

#####################################Newton Method####################################

iter=1
for Vgf in Vgf_array:
  xg1 = Vgb_array[0] / vt
  xg2 = Vgf / vt
  q1guess = -0.17/vt #q1guessf(xn,vt,k1,k2,xg1,xg2,A0)
  q1  = []
  q2  = []
  qtotal = []
  phi1 = []
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
      #print q1guess,f(q1guess,xn,vt,k1,k2,xg1,xg2,A0),iternnew
      iternnew+=1
    #q1guess = fsolve(f, q1guess, args=data,xtol=1.49012e-14)
    #print "solved for Vgb: %f Vbf: %f with %d iterations, q1: %f" % (Vgb,Vgf,iternnew,q1guess*vt)
    
    phi1.append(q1guess)
    q1aux = (xg1-q1guess)
    q2aux = q2value(q1aux,xn,vt,k1,k2,xg1,xg2,A0)
    q1.append(q1aux*vt*cox)
    q2.append(q2aux*vt*cox)
    qtotal.append((q1aux-q2aux)*vt*cox )
  pylab.figure(1)
  pylab.plot( Vgb_array,q1,'o')
  pylab.figure(2)
  pylab.plot( Vgb_array,q2,'o')
  pylab.figure(3)
  pylab.plot( Vgb_array,qtotal,'o')
  pylab.figure(4)
  pylab.plot( Vgb_array,phi1,'o')  

pylab.show()



