#Python Implementation of paper: UTSOI2: A Complete Physical Compact Model for UTBB and Independent Double Gate MOSFETs 
#Juan Duarte (Based on Pragya Matlab implementation)
import numpy as np
from numpy import exp,tan,tanh,log,cosh,sqrt,sinh,sin,cos
import pylab
from scipy.optimize import fsolve

def f(q1,xn,vt,k1,k2,xg1,xg2,A0):
  qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
  if qsqrt >= 0:
    q = sqrt(qsqrt)
    q2  = xg2-xg1+q1+2*log(k1*q1+q*1/tanh(q/2))-log(qsqrt/(sinh(q/2))**2)
    #q2 = (-qsqrt-k1*q1*q/tanh(q/2))/(k2*q/tanh(q/2)+k2*k1*q1)
    equation4 = -A0*exp(-xn)*exp(xg1-q1)+(k1*q1+q*1/tanh(q/2))*(k1*q1+k2*q2)#
  else:
    #tanh(z) = -i tan(iz) 
    #sinh(z) = -i sin(iz)
    #q = np.complex(0,sqrt(-qsqrt))
    #q2  = xg2-xg1+q1+2*log(k1*q1+q*1/tanh(q/2))-log(qsqrt/(sinh(q/2))**2)
    #equation4 = -A0*exp(-xn)*exp(xg1-q1) +(k1*q1+q/tanh(q/2))*(k1*q1+k2*q2)#
    q = sqrt(-qsqrt)   
    q2  = xg2-xg1+q1+2*log(k1*q1+q/tan(-q/2))-log(-qsqrt/(sin(-q/2))**2)
    #q2 = (-qsqrt-k1*q1*q/tan(-q/2))/(k2*q/tan(-q/2)+k2*k1*q1)
    equation4 = -A0*exp(-xn)*exp(xg1-q1)+(k1*q1-q*1/tan(-q/2))*(k1*q1+k2*q2)
    #print k1*q1+q*1/tanh(-q/2)
  #print qsqrt, q, q2, equation4
  return equation4.real

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
#initial guess of front and back surface potential


#print f(x0,xn,vt,k1,k2,xg1,xg2,A0)
#print fprime(x0,xn,vt,k1,k2,xg1,xg2,A0)

# vg2=-6:2:6;
# vg1=-2:0.1:2;
Vg_array = np.linspace(-2, 2, 1000)
Vgf = -6 #front gate fixed in this case
q1  = []
q1guess = -1.1
xg1 = -2 / vt
xg2 = -6 / vt
q1guess = q1guessf(xn,vt,k1,k2,xg1,xg2,A0)
iter=1
for Vgb in Vg_array:
    print "solving: " + str(iter)
    iter+=1
    xg1 = Vgb / vt
    xg2 = Vgf / vt
    data = (xn,vt,k1,k2,xg1,xg2,A0)
    
    q1guess = fsolve(f, q1guess, args=data,xtol=1.49012e-14)
    q1.append(q1guess)

pylab.figure(1)
pylab.plot(Vg_array ,q1) 
#pylab.show()

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


