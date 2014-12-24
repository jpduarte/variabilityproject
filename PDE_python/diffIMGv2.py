#derivative for charge model
from sympy import *
import numpy as np

Vgb = Symbol('Vgb')
Vtb = Symbol('Vtb')
Vgf = Symbol('Vgf')
Vtf = Symbol('Vtf')
vt = Symbol('vt')
qb = Symbol('qb')
qf = Symbol('qf')
rf = Symbol('rf')
rb = Symbol('rb')


qr = qb/qf
cosplusb = exp((Vgb-Vtb-vt*log(qb**2)-rf*qb)/(-2*vt))
oneovercosaminusbsqr = (qr**2*(1-cosplusb**2)/(cosplusb**2)+1)
f = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/vt) -(oneovercosaminusbsqr) 

dfdqf = f.diff(qf)
print "df/dqf"
print dfdqf
"""2*qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**3 + (-rf - 2*vt/qf)*exp((Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)/vt
"""


dfdqb = f.diff(qb)
print "df/dqb"
print dfdqb 
"""-qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*(-rf - 2*vt/qb)*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/(qf**2*vt) - qb**2*(-rf - 2*vt/qb)/(qf**2*vt) - 2*qb*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**2
"""

eta = Symbol('eta')
Vfbf = Symbol('Vfbf')
Vfbb = Symbol('Vfbb')

g = (qf)**2-(qb)**2-eta*(exp((Vgf-Vfbf-rf*qf)/vt)-exp((Vgb-Vfbb-rb*qb)/vt))
dgdqf = g.diff(qf)
print "dg/dqf"
print dgdqf 
#eta*rf*exp((-Vfbf + Vgf - qf*rf)/vt)/vt + 2*qf
dgdqf2 = dgdqf.diff(qf)
print "d2g/dqf2"
print dgdqf2 

dgdqb = g.diff(qb)
print "dg/dqb"
print dgdqb 
#-eta*rb*exp((-Vfbb + Vgb - qb*rb)/vt)/vt - 2*qb
dgdqb2 = dgdqb.diff(qb)
print "d2g/dqb2"
print dgdqb2 

"""
dg/dqf
eta*rf*exp((-Vfbf + Vgf - qf*rf)/vt)/vt + 2*qf
d2g/dqf2
-eta*rf**2*exp((-Vfbf + Vgf - qf*rf)/vt)/vt**2 + 2
dg/dqb
-eta*rb*exp((-Vfbb + Vgb - qb*rb)/vt)/vt - 2*qb
d2g/dqb2
eta*rb**2*exp((-Vfbb + Vgb - qb*rb)/vt)/vt**2 - 2
"""

h = (Vgf-Vgb)-(Vtf-Vtb)+vt*log(qr**2)-(rf*qf-rb*qb)-vt*log(qr**2*(1-cosplusb**2)+cosplusb**2)

dhdqf = h.diff(qf)
print "dh/dqh"
print dhdqf 
"""2*qb**2*vt*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))/(qf**3*(qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))/qf**2 + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))) - rf - 2*vt/qf
"""


dhdqb = h.diff(qb)
print "dh/dqb"
print dhdqb 

"""rb - vt*(qb**2*(-rf - 2*vt/qb)*exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/(qf**2*vt) + 2*qb*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))/qf**2 - (-rf - 2*vt/qb)*exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/vt)/(qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))/qf**2 + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)) + 2*vt/qb
"""

cosaminusb = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/(-2*vt))
oneovercosaminusbsqr = (qr**2*(1-cosplusb**2)/(cosplusb**2)+1)
k = 1/cosaminusb**2-oneovercosaminusbsqr

dkdqf = k.diff(qf)
print "dk/dqh"
print dkdqf 

#2*qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**3 + (-rf - 2*vt/qf)*exp((Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)/vt


dkdqb = k.diff(qb)
print "dk/dqb"
print dkdqb 

#-qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*(-rf - 2*vt/qb)*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/(qf**2*vt) - qb**2*(-rf - 2*vt/qb)/(qf**2*vt) - 2*qb*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**2

baux = (qf-qb)/(cosaminusb/(sqrt((cosaminusb**2-1)))+cosplusb/(sqrt((cosplusb**2-1))))
i = Vgf - Vtf -rf*qf-vt*log(baux**2/(cosaminusb**2-1))

didqf = i.diff(qf)
print "di/dqf"
print didqf 

didqb = i.diff(qb)
print "di/dqb"
print didqb 

"""di/dqf
-rf - vt*(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**2*((-2*qb + 2*qf)/((1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**2) + (-qb + qf)**2*((-rf - 2*vt/qf)*exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/(vt*sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))) + (-rf - 2*vt/qf)*exp(-3*(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/(vt*(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))**(3/2)))/((1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**3) - (-qb + qf)**2*(-rf - 2*vt/qf)*exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)/(vt*(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))**2*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**2))/(-qb + qf)**2
di/dqb
-vt*(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**2*((-qb + qf)**2*((-rf - 2*vt/qb)*exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/(vt*sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))) + (-rf - 2*vt/qb)*exp(-3*(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/(vt*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))**(3/2)))/((1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**3) + (2*qb - 2*qf)/((1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt))*(exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/(2*vt))/sqrt(1 - exp(-(Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)) + exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/(2*vt))/sqrt(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)))**2))/(-qb + qf)**2
"""

