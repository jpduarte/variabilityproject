import numpy as np
import pylab

q = np.linspace(np.complex(0,-(2*np.pi-0.01)),np.complex(0,2*np.pi-0.01),100)

a = q/np.tanh(q/2)
print q,a

pylab.figure(1)
pylab.plot(q.imag ,np.abs(a.real),'o')
ax = pylab.gca()
ax.set_yscale('log')

q2 = np.linspace(np.complex(-(2*np.pi-0.01),0),np.complex(2*np.pi-0.01,0),100)
a2 = q2/np.tanh(q2/2)
pylab.figure(1)
pylab.plot(q2 ,np.abs(a2.real))
ax = pylab.gca()
ax.set_yscale('log')
print q2,a2

      
pylab.show()  

