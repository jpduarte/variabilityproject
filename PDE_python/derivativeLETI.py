from sympy import *
import numpy as np

x = Symbol('x')
k1 = Symbol('k1')
k2 = Symbol('k2')
xg1 = Symbol('xg1')
xg2 = Symbol('xg2')
A0 = Symbol('A0')
xn = Symbol('xn')
q1 = Symbol('q1')
vt = Symbol('vt')


qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
q     = sqrt(qsqrt)
qcoth = q*coth(q/2)
q2    = xg2-xg1+q1+2*log(k1*q1+qcoth)-log(qsqrt/(sinh(q/2))**2)
equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+qcoth)*(k1*q1+k2*q2)

yprime = equation4.diff(q1)
print "derivative case 1:"
print yprime


###########################################################################


qsqrt = k1**2*q1**2-A0*exp(-xn)*exp(xg1-q1)
q     = sqrt(-qsqrt)
qcoth = q*cot(q/2)
q2    = xg2-xg1+q1+2*log(k1*q1+qcoth)-log(qsqrt/(-(sin(q/2))**2))  
equation4 = -A0*exp(-xn+xg1-q1)+(k1*q1+qcoth)*(k1*q1+k2*q2)
yprime = equation4.diff(q1)

print "derivative case 2:"
print yprime
