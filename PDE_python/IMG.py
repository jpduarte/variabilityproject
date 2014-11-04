import numpy as np
from numpy.matlib import repmat
from scipy.misc import factorial
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand

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

"""MATLAB: function u=bvp1(x,f,sigma,beta)
N=length(x);

K(1,1:3)=mkfdstencil(x(1:3),x(1),1);
for i=2:N-1
  K(i,i-1:i+1)=mkfdstencil(x(i-1:i+1),x(i),2);
end
rhs=[sigma; f(x(2:N-1))'];
rhs(end)=rhs(end)-beta*K(N-1,N);
u=[K(1:N-1,1:N-1)\rhs;beta]; """
  
def K_generator(x):
#this return matrix of Poisson Equation in Silicon Fin, with Neuman BC at both ends
  N=len(x);
  K = lil_matrix((N, N))
  K[0,:2]=mkfdstencil(x[0:2],x[0],1)
  i=1
  for xbar in x[1:-1]:
    K[i,i-1:i+2]=mkfdstencil(x[i-1:i+2],xbar,2)
    i+=1
  K[i,i-1:i+1]=mkfdstencil(x[i-1:i+1],x[i],1)  
  return K.tocsr()
  
def rho_phi(phi,vt,ni):
  return ni*np.exp(-phi/vt)

#################################################################  
x     = np.array([1, 2, 3, 4])
xbar  = 2
k     = 2
Tsi   = 1e-9
print mkfdstencil(x,xbar,k)
print K_generator(x)

