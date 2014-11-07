#11/5/2014 first version
#Juan Pablo Duarte
"""This code solve the Poisson Equation of a independent double gate finfet
in 1-D, assuming Boltzman charge distribution.
It is solved using finite difference method in python 2.6.6 version
To run it: python IMG.py
numpy, scipy and pylab libraries are needed
TODO: it can be implemented in a more efficient way, for example, K generated only once, etc."""

import numpy as np
from numpy.matlib import repmat
from scipy.misc import factorial
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
import pylab
from scipy import integrate
from scipy.integrate import quad

def mkfdstencil(x,xbar,k):
  maxorder            = len(x)
  h_matrix            = repmat(np.transpose(x)-xbar,maxorder,1)
  powerfactor_matrix  = np.transpose(repmat(np.arange(0,maxorder),maxorder,1))
  factorialindex      = np.transpose(repmat(factorial(np.arange(0,maxorder)),maxorder,1))
  taylormatrix        = h_matrix ** powerfactor_matrix /factorialindex
  derivativeindex     = np.zeros(maxorder)
  derivativeindex[k]  = 1
  u = np.linalg.solve(taylormatrix,derivativeindex)
  return u

def K_generator(x):
#this return matrix of Poisson Equation in Silicon Fin, with Neuman BC at both ends
  N=len(x);
  K = lil_matrix((N, N))
  K[0,:5]=mkfdstencil(x[0:5],x[0],1)
  i=1
  for xbar in x[1:-1]:
    K[i,i-1:i+2]=mkfdstencil(x[i-1:i+2],xbar,2)
    i+=1
  K[i,i-6:i+1]=mkfdstencil(x[i-6:i+1],x[i],1)  
  return K.tocsr()
  
def Efield_matrix(x):
#this return matrix of Poisson Equation in Silicon Fin, with Neuman BC at both ends
  N=len(x);
  K = lil_matrix((N, N))
  K[0,:5]=mkfdstencil(x[0:5],x[0],1)
  i=1
  for xbar in x[1:-1]:
    K[i,i-2:i+3]=mkfdstencil(x[i-1:i+2],xbar,1)
    i+=1
  K[i,i-6:i+1]=mkfdstencil(x[i-6:i+1],x[i],1)  
  return K.tocsr()  
  
def rho_phi(phi,q,vt,ni,Nch,ech):
  return -(-q*(ni**2/Nch)*np.exp(phi/vt)-q*Nch)/ech

def drho_dphi(phi,q,vt,ni,Nch,ech):
  return q*(ni**2/Nch)*np.exp(phi/vt)/(vt*ech)
  
def solvePoisson(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N,phi,x,K):
    Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    #x       = Tsi/2*np.linspace(-1, 1, N)
    rhs     = np.zeros(N)
    jrhs    = np.zeros(N)
    J = lil_matrix((N, N))
    ################################
    #K = K_generator(x)
    g   = 1000
    iter=1
    while True:
        iter+=1
        if (norm(g)/N**2) <1e-16 or iter> 50:
            print "Solve in iterations:" + str(iter)
            print "Solve with norm:" + str(norm(g)/N**2)
            break  
        rhs[1:N]    = rho_phi(phi[1:N],q,vt,ni,Nch,ech)
        rhs[0]      = -eins/(ech*tins)*(Vgf-Vfbf-phi[0])
        rhs[-1]     = eins/(ech*tins)*(Vgb-Vfbb-phi[-1])

        jrhs[1:N]   = drho_dphi(phi[1:N],q,vt,ni,Nch,ech)
        jrhs[0]     = -eins/(ech*tins)*(-1)
        jrhs[-1]    = eins/(ech*tins)*(-1)
       
        J.setdiag(jrhs)
         
        g   = K*phi-rhs
        dg  = (K-J.tocsr())
        phi = phi -  solve(dg.todense(),g)
    return phi
    
def gamma_IMG(phix,Vg,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gate,phi_substrate,Eg,N):
    Vfb     = phi_gate - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    alpha   = vt*2*q*ni**2/(ech*Nch)
    beta    = 2*q*Nch/ech
    F       = (eins/tins)*(Vg-Vfb-phix)/ech
    #print "F: %e" % (F)
    return  F**2-alpha*np.exp(phix/vt)-beta*phix
   
def charge_channel(phi,x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N):
    Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    y = q*ni**2/(Nch)*np.exp(phi/vt)
    total_charge_integration = integrate.cumtrapz(y, x)
    total_charge_gauss =  eins/(tins)*(Vgf-Vfbf-phi[0])+eins/(tins)*(Vgb-Vfbb-phi[-1])
    #print total_charge
    return  total_charge_integration[-1],total_charge_gauss 
    
def dqdphi(phi,phif,phib,x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N):
    Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    alpha   = vt*2*q*ni**2/(ech*Nch)
    beta    = 2*q*Nch/ech
    gamma    = gamma_IMG(phib,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gateb,phi_substrate,Eg,N)
    
    #print total_charge
    return  -q*ni**2/(Nch)*np.exp(phi/vt)/np.sqrt(gamma+alpha*np.exp(phi/vt)+beta*phi)  

def Efield_analytical(phi,phif,phib,x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N):
    Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*np.log(Nch/ni)
    alpha   = vt*2*q*ni**2/(ech*Nch)
    beta    = 2*q*Nch/ech
    gamma    = gamma_IMG(phib,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gateb,phi_substrate,Eg,N)
    
    #print total_charge
    return  1/np.sqrt(gamma+alpha*np.exp(phi/vt)+beta*phi)     
    
##########################PARAMETERS##############################  
#physical constants
q       = 1.6e-19
vt      = 0.0259
ni      = 1e10
k       = 1.380650e-23 
T       = 300 
eo      = 8.854000e-14 
eins    = 3.453060e-13 
ech     = 1.035918e-12 
Eg 	    = 1.169640 
Nc      = 2.540000e+19 
Nv 	    = 3.140000e+19 

#device dimensions and parameters
Tsi             = 10e-9
tins            = 1e-9
Nch             = 1e18
phi_substrate   = 4.05 
phi_gatef 	    = 4.5
phi_gateb 	    = 4.5

#############################SIMULATION##############################
Vg_array = np.linspace(0, 2, 5)
Vgf = 1.0 #front gate fixed in this case

N = 3
"""
phiguess = np.zeros(N)
for Vgb in Vg_array:
    print "Solve for Vgf %5.3f and Vg %5.3f, using %d nodes" % (Vgf,Vgb,N)
    phiguess, x = solvePoisson(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N,phiguess)
    pylab.plot(np.linspace(-Tsi/2, Tsi/2, N) ,phiguess)"""
    
N = 500
phiguess = np.zeros(N)
x       = Tsi/2*np.linspace(-1, 1, N)
K = K_generator(x)
Efield_mtx = Efield_matrix(x)
for Vgb in Vg_array:
    print "Solve for Vgf %5.3f and Vg %5.3f, using %d nodes" % (Vgf,Vgb,N)
    phiguess = solvePoisson(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N,phiguess,x,K)
    pylab.figure(1)
    pylab.plot(x ,phiguess,'o') 
    pylab.figure(2)
    pylab.plot(x ,q*ni**2/(Nch)*np.exp(phiguess/vt),'o') 
    Efield = Efield_mtx*phiguess
    pylab.figure(3)
    pylab.plot(x ,Efield,'o') 
    pylab.figure(4)
    pylab.plot(x ,Efield_analytical(phiguess,phiguess[0],phiguess[-1],x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N),'o') 
    
    print "gamma front: %e" % (gamma_IMG(phiguess[0],Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_substrate,Eg,N) )
    print "gamma back: %e" % (gamma_IMG(phiguess[-1],Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gateb,phi_substrate,Eg,N) )
    print "Total Charge: %e, %e" % (charge_channel(phiguess,x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N))
    charge = quad(dqdphi, phiguess[0], phiguess[-1], args=(phiguess[0], phiguess[-1],x,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,phi_gatef,phi_gateb,phi_substrate,Eg,N))
    print "Total Charge Integration: %e"% ( charge[0])
pylab.show()    
