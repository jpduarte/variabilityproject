#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
import numpy as np
from numpy.matlib import repmat
#from numpy import exp,tan,tanh,log,cosh,sqrt,sinh,sin,cos
import pylab
from scipy.optimize import fsolve
from sympy import exp,tan,tanh,log,cosh,sqrt,sinh,sin,cos

##############################################################################################################
def coth(x):
  return 1/tanh(x)

def f(q1,xn,vt,k1,k2,xg1,xg2,A0):
  #x = complex(q1,0)
  x = q1
  term1 = (k1*(x - xg1) - coth((k1**2*(x - xg1)**2 - A0*exp(x - xn))**0.5/2.0)*(k1**2*(x - xg1)**2 - A0*exp(x - xn))**0.5)*(k1*(x - xg1) + k2*(x - xg2 + log((k1**2*(x - xg1)**2 - A0*exp(x - xn))/sinh((k1**2*(x - xg1)**2)/4.0 - (A0*exp(x - xn))/4.0)) - 2.0*log(coth((k1**2*(x - xg1)**2 - A0*exp(x - xn))**0.5/2.0)*(k1**2*(x - xg1)**2 - A0*exp(x - xn))**(1./2) - k1*(x - xg1)))) - A0*exp(x - xn)
  return term1

"""(k1.*(x - xg1) - coth((k1.^2.*(x - xg1).^2 - A0.*exp(x - xn)).^(1./2)./2).*(k1.^2.*(x - xg1).^2 - A0.*exp(x - xn)).^(1./2)).*(k1.*(x - xg1) + k2.*(x - xg2 + log((k1.^2.*(x - xg1).^2 - A0.*exp(x - xn))./sinh((k1.^2.*(x - xg1).^2)./4 - (A0.*exp(x - xn))./4)) - 2.*log(coth((k1.^2.*(x - xg1).^2 - A0.*exp(x - xn)).^(1./2)./2).*(k1.^2.*(x - xg1).^2 - A0.*exp(x - xn)).^(1./2) - k1.*(x - xg1)))) - A0.*exp(x - xn)"""
  
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
  xsub1 = (F1/2)-np.sqrt(h1) 
  neta1 = 0.5*(xsub1+xn+13-np.sqrt(64+(xsub1-xn-13)**2))
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

Vgb_array = np.linspace(0, 2, 50)
Vgf_array = np.linspace(-6, -6, 1) #front gate fixed in this case

xg1 = -2 / vt
xg2 = -6 / vt
q1guess = q1guessf(xn,vt,k1,k2,xg1,xg2,A0)
print q1guess,f(q1guess,xn,vt,k1,k2,xg1,xg2,A0)

"""
iter=1
for Vgf in Vgf_array:
  xg1 = Vgb_array[0] / vt
  xg2 = Vgf / vt
  q1guess = q1guessf(xn,vt,k1,k2,xg1,xg2,A0)
  q1  = []
  print "solving for Vgf: %f, with initial guess %f" % (Vgf,q1guess)
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
    print "solved for Vgb: %f Vbf: %f with %d iterations, q1: %f" % (Vgb,Vgf,iternnew,q1guess*vt)
    
    q1.append(q1guess*vt)
  pylab.figure(1)
  pylab.plot( Vgb_array,q1)

pylab.show()"""
"""
xg1 = 1 / vt
xg2 = 1 / vt
     
x0_array = np.linspace(-5, 5, 100)
fvalue = []
for x0 in x0_array:
  fvalue.append(f(x0,xn,vt,k1,k2,xg1,xg2,A0))

data = (xn,vt,k1,k2,xg1,xg2,A0)
solutionq1 = fsolve(f, -1.1, args=data,xtol=1.49012e-14)   
print  solutionq1, f(solutionq1, xn,vt,k1,k2,xg1,xg2,A0)
pylab.figure(2)
pylab.plot(x0_array ,fvalue) 
pylab.show()


#data = (xn,vt,k1,k2,xg1,xg2,A0)
#q1sol = fsolve(f, 6.2, args=data,xtol=1.49012e-14)
#print q1sol, f(q1sol,xn,vt,k1,k2,xg1,xg2,A0)
   """




"""tolerance = 10^(-10);         %7 digit accuracy is desired
epsilon = 10^(-5) ;         %Don't want to divide by a number smaller than this
 
maxIterations = 20 ;         %Don't allow the iterations to continue indefinitely
haveWeFoundSolution = false ;%Were we able to find the solution to the desired tolerance? not yet.
 
for i = 1 : maxIterations 
 
    y = f(x0);
    yprime =fprime(x0);

    if(abs(yprime) < epsilon)                         %Don't want to divide by too small of a number
        fprintf('WARNING: denominator is too small\n')
        break;                                        %Leave the loop
    end
 
    x1 = x0 - y/yprime   ;                             %Do Newton's computation
 
    if(abs(x1 - x0)/abs(x1) < tolerance)              %If the result is within the desired tolerance
        haveWeFoundSolution = true;
        break;                                        %Done, so leave the loop
    end
 
    x0 = x1 ;                                    %Update x0 to start the process again                  
 
end


plot(vg1,x1*vt)  % correct trend but sharp rise in x2
hold on
grid
end
xlim([0,1.4])"""

"""term1 =(k1*(x - xg1) - 1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1.0/2.0)/2.0)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1.0/2.0))*(k1*(x - xg1) + k2*(x - xg2 + np.log((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))/np.sinh((k1**2*(x - xg1)**2)/4.0 - (A0*np.exp(x - xn))/4.0)) - 2.0*np.log(1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1.0/2.0)/2.0)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2) - k1*(x - xg1)))) - A0*np.exp(x - xn)
def fprime(x,xn,vt,k1,k2,xg1,xg2,A0):
  return (k1*(x - xg1) + k2*(x - xg2 + np.log((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))/np.sinh((k1**2*(x - xg1)**2)/4 - (A0*np.exp(x - xn))/4)) - 2*np.log(1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2) - k1*(x - xg1))))*(k1 + ((1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)**2 - 1)*(k1**2*(2*x - 2*xg1) - A0*np.exp(x - xn)))/4 - (1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)*(k1**2*(2*x - 2*xg1) - A0*np.exp(x - xn)))/(2*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2))) - A0*np.exp(x - xn) + (k1 + k2*((np.sinh((k1**2*(x - xg1)**2)/4 - (A0*np.exp(x - xn))/4)*((k1**2*(2*x - 2*xg1) - A0*np.exp(x - xn))/np.sinh((k1**2*(x - xg1)**2)/4 - (A0*np.exp(x - xn))/4) - (np.cosh((k1**2*(x - xg1)**2)/4 - (A0*np.exp(x - xn))/4)*((k1**2*(2*x - 2*xg1))/4 - (A0*np.exp(x - xn))/4)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn)))/np.sinh((k1**2*(x - xg1)**2)/4 - (A0*np.exp(x - xn))/4)**2))/(k1**2*(x - xg1)**2 - A0*np.exp(x - xn)) - (2*(k1 + ((1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)**2 - 1)*(k1**2*(2*x - 2*xg1) - A0*np.exp(x - xn)))/4 - (1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)*(k1**2*(2*x - 2*xg1) - A0*np.exp(x - xn)))/(2*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2))))/(k1*(x - xg1) - 1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)) + 1))*(k1*(x - xg1) - 1/np.tanh((k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2)/2)*(k1**2*(x - xg1)**2 - A0*np.exp(x - xn))**(1./2))"""


