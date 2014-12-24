#Implementation of paper: An Analytic Potential Model for Symmetric and Asymmetric DG MOSFETs
#Implementation of Taur Model for IMG, using independent front and back gates
#Juan Pablo Duarte
#TODO: check initial guess for roboust aand b newton method, heuristic is fine
#TODO: sweep Vgf by Vgb depending on which one is larger

from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace,array,pi,matrix,dot
from numpy.linalg import norm,solve,inv,det
import pylab

def csch(x):
  return 1/sinh(x)
  
def coth(x):
  return 1/tanh(x)
  
def csc(x):
  return 1/sin(x)
  
def cot(x):
  return 1/tan(x)
  
def unified_charge_model(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch) :
  rc  = (2*Cins/(Weff**2*ech/Ach))
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2*2*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2/(exp(qdep*rc)-qdep*rc-1))+vth_fixed_factor_SI
  vth_N_Sub = -qdep+vth_fixed_factor_Sub
  vth_N_SI  = -qdep+vth_fixed_factor_SI
  Vg_local_N  = (Vg-Vch)*vt**(-1)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5
  if (Vov>60):
    qm      = -Vov*2
    qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
    x0      = qtrc*(exp(qtrc)-qtrc-1)**(-1)
    x1      = qtrc*x0
    f0      = F-qm+log(-qm)+log(x1)
    f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
    f2      = -(qm**2)**(-1)
    qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5)
    if(qm>1e-7):
      qm      = 2*(1-sqrt(1+(log(1+qm))**2))
      qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
      x0      = qtrc*(exp(qtrc)-qtrc-1)**(-1)
      x1      = qtrc*x0
      f0      = F-qm+log(-qm)+log(x1)
      f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
      f2      = -(qm**2)**(-1)
      qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
    else:
      qm      = -qm**2
  return qm  

def vcritical(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)#+(Vgf-Vgb)
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
  Vcritical = Vfbf+C1+2*vt*log(2/(s-1))+2*vt*r*(2/(s-1))
  print "Vcritical: "+ str(Vcritical)
  return Vcritical,s

def ab(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vcritical,s = vcritical(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
  
  if Vgf < Vcritical:
    a,b,phisf,phisb = abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(coth(a-b)-coth(a+b)) 
    case = 1
  else:
    a,b,phisf,phisb = abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(cot(a-b)-cot(a+b))
    case = 2  
  return a,b,phisf,phisb,Qt,Vcritical,s,case
  
def abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)#+(Vgf-Vgb)
  
  #initial guess    
  a=aguess
  b=bguess 
  s_m=array([a,b])

  Y = array([1,1])
  iter=0
  print "solving abcase1 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
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
    print det(J) 
    if  abs(det(J)  )<1e-10:
      a=aguess*10
      b=bguess*10 
    else:
      s_m = s_m - solve(J,Y) 
    iter+=1
    print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((-b*Tsi/Tsi)+a))
   
  return a,b,phisf,phisb

def abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)#+(Vgf-Vgb)
  
  #initial guess    
  a=aguess
  b=bguess 
  s_m=array([a,b])

  Y = array([1,1])
  iter=0
  print "solving abcase2 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
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
    print det(J)
    if  abs(det(J)  )<1e-10:
      a=aguess*10
      b=bguess*10 
    else:
      s_m = s_m - solve(J,Y)  
    iter+=1
    print iter,a,b,norm(Y)
    
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
Tsi             = 10e-7
Hfin            = 1e-6
tins            = 1e-7
tinsbox         = tins
Nch             = 1e10
phi_substrate   = 4.05 
phi_gatef 	    = 5.0
phi_gateb 	    = 4.0

#Vth1 calculation, Vthback
alpha_MI = 1.6
Cins = (Hfin)*eins/tins
Ach=Hfin*(Tsi/2)
Weff = Hfin

Vgb_array = linspace(-1,2,500)
Vgf = -0.5
a = []
b = []
phisf = []
phisb = []
Qt = []
Vcritical = []
s = []
case = []
aaux,baux = pi/9,pi/10
qmb = []
qmf = []
qmt = []
for Vgb in Vgb_array:
  #Vgf=Vgb
  if abs((Vgf-phi_gatef)-(Vgb-phi_gateb))<1e-10:
    phi_gateb+=1e-5
    
  if (Vgf-phi_gatef)>(Vgb-phi_gateb):
    aaux,baux,phisfaux,phisbaux,Qtaux,Vcriticalaux,saux,caseaux = ab(aaux,baux ,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
  else:
    aaux,baux,phisbaux,phisfaux,Qtaux,Vcriticalaux,saux,caseaux = ab(aaux,baux ,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg)
  a.append(aaux)
  b.append(baux)
  phisf.append(phisfaux)
  phisb.append(phisbaux) 
  Qt.append(Qtaux) 
  Vcritical.append(Vcriticalaux)
  s.append(saux)
  case.append(caseaux)
  qmbaux = -(vt*Cins/Weff)*unified_charge_model(Vgf,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gateb,alpha_MI,Cins,Ach,Weff,Nch)
  qmfaux = -(vt*Cins/Weff)*unified_charge_model(Vgf,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,alpha_MI,Cins,Ach,Weff,Nch)
  qmb.append(qmbaux)
  qmf.append(qmfaux)
  qmt.append(qmfaux+qmbaux)    
  
#delta Vth1
tins_p=(ech/eins)*tins
slope_s=(-(phi_gateb-phi_gatef)/vt)/(2*tins_p+Tsi)
v1a = (exp((tins_p+Tsi)*slope_s)-exp((tins_p)*slope_s))/slope_s
deltaVth1 = -vt*log((2/Tsi)*v1a)
print "delta Vth1: " + str(deltaVth1)






Vfbf = (phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni))
Vthfaux = Vfbf+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))+vt*(-q*Nch*Ach)/(vt*Cins)
print "Vth1 original: " + str(Vthfaux) 
print "Vth1 final: " + str(Vthfaux+deltaVth1)  

#Vth2 calculation, Vthback
alpha_MI = 1.6
Cins = (2*Hfin)*eins/tins
Ach=Hfin*(Tsi/2)
Weff = 2*(Tsi/2)

Vfbb = (phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni))+vt*(-q*Nch*Ach)/(vt*Cins)
Vthb = Vfbb+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))
print "Vth2: " + str(Vthb)  
print phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
print phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)



##########plots  
pylab.figure(1)
pylab.plot(Vgb_array ,a,'r') 
pylab.plot(Vgb_array ,b,'b') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("alpha (red) and beta (blue)", fontsize=18)

pylab.figure(2)
pylab.plot(Vgb_array ,a,'r') 
pylab.plot(Vgb_array ,b,'b') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("alpha (red) and beta (blue)", fontsize=18)
ax = pylab.gca()
ax.set_yscale('log')

#pylab.plot(Vgb_array ,[x*-1 for x in a] ,'r') 
pylab.figure(3)
pylab.plot(Vgb_array ,phisf,'r') 
pylab.plot(Vgb_array ,phisb,'b') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("Front and back potentials (V)", fontsize=18)

pylab.figure(4)
pylab.plot(Vgb_array ,Qt,'r')
pylab.plot(Vgb_array ,qmf,'b--')  
pylab.plot(Vgb_array ,qmb,'k--') 
pylab.plot(Vgb_array ,qmt,'g--') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("Qt", fontsize=18)

pylab.figure(5)
pylab.plot(Vgb_array ,Qt,'r') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("Qt", fontsize=18)
pylab.plot(Vgb_array ,qmf,'b--')  
pylab.plot(Vgb_array ,qmb,'k--') 
pylab.plot(Vgb_array ,qmt,'g--')  
ax = pylab.gca()
ax.set_yscale('log')

pylab.figure(6)
pylab.plot(Vgb_array ,Vcritical,'r') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("Vcritical", fontsize=18)

pylab.figure(7)
pylab.plot(Vgb_array ,s,'r') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("s", fontsize=18)
ax = pylab.gca()
#ax.set_yscale('log')

pylab.figure(8)
pylab.plot(Vgb_array ,case,'r') 
pylab.xlabel("Vgb (V)", fontsize=18)
pylab.ylabel("case", fontsize=18)

pylab.figure(9)
pylab.plot(a ,b,'k') 

pylab.xlabel("alpha ", fontsize=18)
pylab.ylabel("beta", fontsize=18)

pylab.show()  


#print vcritical(0.2,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
#print ab(1.1*pi/9,pi/10,1.0,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
