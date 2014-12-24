#Implementation of paper: An Analytic Potential Model for Symmetric and Asymmetric DG MOSFETs
#Implementation of Taur Model for IMG, using tight front and back gates
#Juan Pablo Duarte
from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace,array,pi,matrix,dot
from numpy.linalg import norm,solve,inv
import pylab

def csch(x):
  return 1/sinh(x)
  
def coth(x):
  return 1/tanh(x)
  
def csc(x):
  return 1/sin(x)
  
def cot(x):
  return 1/tan(x)

def vcritical(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  s=1.01
  f=1
  iter=1
  while iter<100:
    if abs(f)<1e-12:
      break
    f   =log((s+1)/(s-1))+r*((s+1)/(s-1)-(s-1)/(s+1))-(Vfbb-Vfbf)/(2*vt)
    f_d1=1/(s+1)-1/(s-1)
    f_d2=r*(-4*(1+s**2))/(s**2-1)**2
    f_d =f_d1+f_d2
    s=(s-f/f_d)
    iter+=1
    print s,f,iter
  return Vfbf+C1+2*vt*log(2/(s-1))+2*vt*r*(2/(s-1))

def ab(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vcritical = vcritical(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  
  if Vgf < Vcritical:
    a,b,phisf,phisb = abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(coth(a-b)-coth(a+b)) 
  else:
    a,b,phisf,phisb = abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(cot(a-b)-cot(a+b))
      
  return a,b,phisf,phisb,Qt
  
def abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  
  #initial guess    
  a=aguess
  b=bguess 
  s_m=array([a,b])

  Y = array([1,1])
  iter=0
  while iter<100:
    a = s_m[0]
    b = s_m[1]
    if norm(Y)<1e-12:
      break  
    f=log(sinh(a+b)/sinh(a-b))+2*r*b*(coth(a-b)+coth(a+b))-(Vfbb-Vfbf)/(2*vt)
    g=Vfbf+C1+2*vt*log(2*b/(sinh(a-b)))+4*vt*r*b*coth(a-b)-Vg
    
    Y = array([f,g])
      
    f_da=-coth(a-b)+coth(a+b)-2*b*r*((csch(a-b))**2+(csch(a+b))**2)
    f_db=(1+2*r)*coth(a-b)+(1+2*r)*coth(a+b)+2*b*r*((csch(a-b))**2-(csch(a+b))**2)
    g_da=-2*vt*coth(a-b)-4*b*r*vt*(csch(a-b))**2
    g_db=(2*vt/b)*(1+b*(1+2*r)*coth(a-b)+2*b**2*r*(csch(a-b))**2)
    J = array([[f_da,f_db],[g_da,g_db]])
    s_m = s_m - solve(J,Y) 
    iter+=1
    #print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((-b*Tsi/Tsi)+a))
   
  return a,b,phisf,phisb

def abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  
  #initial guess    
  a=aguess
  b=bguess 
  s_m=array([a,b])

  Y = array([1,1])
  iter=0
  while iter<100:
    a = s_m[0]
    b = s_m[1]
    if norm(Y)<1e-12:
      break  
     
    f=log(sin(a+b)/sin(a-b))+2*r*b*(cot(a-b)+cot(a+b))-(Vfbb-Vfbf)/(2*vt)
    g=Vfbf+C1+2*vt*log(2*b/(sin(a-b)))+4*vt*r*b*cot(a-b)-Vg
    
    Y = array([f,g])
      
    f_da=-cot(a-b)+cot(a+b)-2*b*r*((csc(a-b))**2+(csc(a+b))**2)
    f_db=(1+2*r)*cot(a-b)+(1+2*r)*cot(a+b)+2*b*r*((csc(a-b))**2-(csc(a+b))**2)
    g_da=-2*vt*cot(a-b)-4*b*r*vt*(csc(a-b))**2
    g_db=(2*vt/b)*(1+b*(1+2*r)*cot(a-b)+2*b**2*r*(csc(a-b))**2)
    J = array([[f_da,f_db],[g_da,g_db]])
    s_m = s_m - solve(J,Y) 
    iter+=1
    #print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((-b*Tsi/Tsi)+a))    
  return a,b,phisf,phisb

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
tinsbox         = 1e-9
Nch             = 1e15
phi_substrate   = 4.05 
phi_gatef 	    = 4.0
phi_gateb 	    = 5.0

Vg_array = linspace(-1,2,500)
a = []
b = []
phisf = []
phisb = []
Qt = []
aaux,baux = pi/9,pi/10
for Vg in Vg_array:
  aaux,baux,phisfaux,phisbaux,Qtaux = ab(aaux,baux ,Vg,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
  a.append(aaux)
  b.append(baux)
  phisf.append(phisfaux)
  phisb.append(phisbaux) 
  Qt.append(Qtaux) 
pylab.figure(1)
pylab.plot(Vg_array ,a,'r') 
pylab.plot(Vg_array ,b,'b') 

pylab.plot(Vg_array ,[x*-1 for x in a] ,'r') 
pylab.figure(2)
pylab.plot(Vg_array ,phisf,'r') 
pylab.plot(Vg_array ,phisb,'b') 

pylab.figure(3)
pylab.plot(Vg_array ,Qt,'r') 


pylab.show()  


#print vcritical(0.2,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
#print ab(1.1*pi/9,pi/10,1.0,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
