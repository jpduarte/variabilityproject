#Implementation of paper: An Analytic Potential Model for Symmetric and Asymmetric DG MOSFETs
#Implementation of Taur Model for IMG, using same front and back gates volatages
#Juan Pablo Duarte
#Implement a decouple model using unified model
#TODO: make front and back voltages not the same

from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace,array,pi,matrix,dot, arccos, arccosh
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
    if(qm>1e-3):#original 1e-7
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
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)+(Vgf-Vgb)
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
    a,b,phisf,phisb,iter = abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(coth(a-b)-coth(a+b)) 
    Qf = 4*ech*vt/(Tsi)*b*(coth(a-b))
    Qb = 4*ech*vt/(Tsi)*b*(-coth(a+b)) 
    case = 1
  else:
    a,b,phisf,phisb,iter = abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    Qt = 4*ech*vt/(Tsi)*b*(cot(a-b)-cot(a+b))
    Qf = 4*ech*vt/(Tsi)*b*(cot(a-b))
    Qb = 4*ech*vt/(Tsi)*b*(-cot(a+b))     
    case = 2  
  return a,b,phisf,phisb,Qt,Vcritical,s,case,Qf,Qb,iter
  
def abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)+(Vgf-Vgb)
  """
  #initial guess  
  if aguess>0.99*pi/2:
    a=0.98*pi/2
  else:  
    a=aguess
  if bguess>0.99*pi/2:
    b=0.95*pi/2
  else:  
    b=bguess"""  
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
      a=aguess*0.99
      b=bguess*0.9
    else:
      s_m = s_m - solve(J,Y)
    #s_m = s_m - solve(J,Y)
    iter+=1
    print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((-b*Tsi/Tsi)+a))
   
  return a,b,phisf,phisb,iter

def abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):  
  Vg = Vgf
  r        = ech*tins/(eins*Tsi)
  C1       = vt*log(2*ech*vt/(q*ni*Tsi**2))
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)+(Vgf-Vgb)
  
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
      a=aguess*0.99
      b=bguess*0.9
    else:
      s_m = s_m - solve(J,Y)
    #s_m = s_m - solve(J,Y)       
    iter+=1
    print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((-b*Tsi/Tsi)+a))    
  return a,b,phisf,phisb,iter

def unified_IMG_model(Vgf,Vgb,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,phi_gateb,Tsi,tins,tinsbox,Hfin):
  sweepflag = 0
  if (Vgf-phi_gatef) > (Vgb-phi_gateb):
    Vaux = Vgf
    Vgf = Vgb
    Vgb = Vaux
    Vaux = phi_gatef
    phi_gatef = phi_gateb
    phi_gateb = Vaux
    sweepflag = 1
  
  deltaVg = Vgb - Vgf
  phi_gateb = phi_gateb - deltaVg
  
  Vg = Vgf
  
  alpha_MI = 1.6#1.6
  Cins = (Hfin)*eins/tins
  Ach=Hfin*(Tsi/2)
  Weff = Hfin
  tins_p=(ech/eins)*tins
  deltaphi = phi_gateb-phi_gatef#-(sqrt((phi_gateb-phi_gatef)**2+0.01)) #this is bias dependent 
  slope_s=(deltaphi/vt)/(2*tins_p+Tsi)
  v1a = (exp((tins_p+Tsi)*slope_s)-exp((tins_p)*slope_s))/(slope_s)
  deltaVth1 = -vt*log((v1a/Tsi))
  #print "delta Vth1: " + str(deltaVth1)
  
  

  Vfbb = (phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni))
  Vthfaux = Vfbb+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))+vt*(-q*Nch*Ach)/(vt*Cins)
  Vth1 = Vthfaux+deltaVth1
  #print "Vth1 original: " + str(Vthfaux) 
  #print "Vth1 final: " + str(Vth1)  

  Vfbf = (phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni))+vt*(-q*Nch*Ach)/(vt*Cins)
  Vth2 = Vfbf+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))
  #print "Vth2: " + str(Vth2)  

  qmatvth2 = unified_charge_model(Vth2,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gateb,alpha_MI,Cins,Ach,Weff,Nch)
  qmftvth2 = unified_charge_model(Vth2,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,alpha_MI,Cins,Ach,Weff,Nch)
  qmatvth2new = unified_charge_model(Vth2-deltaVth1,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gateb,alpha_MI,Cins,Ach,Weff,Nch)

  Cins1 = Cins*(qmatvth2+0*qmftvth2)/qmatvth2new #unified_charge_model_inverse(qmatvth2,Vth2-
  Cins2 = 2*Cins-Cins1
  
  qmb = -(vt*Cins1/Weff)*unified_charge_model(Vgf-deltaVth1,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gateb,1.5,Cins,Ach,Weff,Nch)
  qmf = -(vt*Cins2/Weff)*unified_charge_model(Vgf,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,1.5,Cins,Ach,Weff,Nch)
  #qmb = -(vt*Cins/Weff)*unified_charge_model(Vgf-0*deltaVth1,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gateb,alpha_MI,Cins,Ach,Weff,Nch)
  #qmf = -(vt*Cins/Weff)*unified_charge_model(Vgf,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,alpha_MI,Cins,Ach,Weff,Nch) 
  
  qmtotal = qmb + qmf
  
  qbguess = -slope_s*vt*ech+qmb*Cins/Cins1
  qfguess = slope_s*vt*ech+(qmtotal-qmb*Cins/Cins1)
  
  #qbguess = max(qmb,-(slope_s*vt*ech))
  #qfguess = qmf+slope_s*vt*ech
      
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  xi = q*ni/(2*ech*vt)
  rf = tins*ech*2*vt*2/(eins*Tsi)
  rb = tins*ech*2*vt*2/(eins*Tsi)
  Vtf = Vfbf-2*vt*log((Tsi/2)*sqrt(xi))
  Vtb = Vfbb-2*vt*log((Tsi/2)*sqrt(xi))
  eta = 2*q*ni*vt/ech/(4*vt/Tsi)**2
  
  qf = qbguess
  qb = qfguess
  x = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/(-2*vt))
  y = exp((Vgb-Vtb-vt*log(qb**2)-rf*qf)/(-2*vt))
  if abs(x)<1:
    amb = arccosh(x)
    apb = arccosh(y)
    aguess = (amb+apb)/2
    bguess = (-amb+apb)/2
  else:
    amb = arccos(x)
    apb = arccos(y)
    aguess = (amb+apb)/2
    bguess = (-amb+apb)/2
  
  """if sweepflag == 1:
    qaux = qfguess
    qfguess = qbguess
    qbguess = qaux"""
      
  return qmb,qmf,Cins1/Cins,deltaVth1,deltaphi,qfguess,qbguess,aguess,bguess
  

def qfqb(qfguess,qbguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  #fix qb then solve for qf
  qf  = qfguess/(4*ech*vt/(Tsi))
  qb = qbguess/(4*ech*vt/(Tsi))
  
  if (Vgf-phi_gatef) > (Vgb-phi_gateb):
    Vaux = Vgf
    Vgf = Vgb
    Vgb = Vaux
    Vaux = phi_gatef
    phi_gatef = phi_gateb
    phi_gateb = Vaux
    #qb  = qfguess/(4*ech*vt/(Tsi))
    #qf = qbguess/(4*ech*vt/(Tsi))

    
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  xi = q*ni/(2*ech*vt)
  rf = tins*ech*2*vt*2/(eins*Tsi)
  rb = tins*ech*2*vt*2/(eins*Tsi)
  Vtf = Vfbf-2*vt*log((Tsi/2)*sqrt(xi))
  Vtb = Vfbb-2*vt*log((Tsi/2)*sqrt(xi))
  eta = 2*q*ni*vt/ech/(4*vt/Tsi)**2

  #Newtonmethod
 

  g=1
  iter=0
  #print "solving abcase2 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
  while iter<0:

    #if norm(g)<1e-10:
    #  break  
     
    g = (qf)**2-(qb)**2-eta*exp((Vgb-Vfbb-rb*qb)/vt)*(exp((Vgf-Vfbf-rf*qf)/vt-(Vgb-Vfbb-rb*qb)/vt)-1)
    dgqgb = -eta*rb*exp((-Vfbb + Vgb - qb*rb)/vt)/vt - 2*qb
    dgqgb2 = eta*rb**2*exp((-Vfbb + Vgb - qb*rb)/vt)/vt**2 - 2
    
    dgqgf = eta*rf*exp((-Vfbf + Vgf - qf*rf)/vt)/vt + 2*qf
    dgqgf2 = -eta*rf**2*exp((-Vfbf + Vgf - qf*rf)/vt)/vt**2 + 2
    #qf = qf - g*(2*dgqgf/(2*dgqgf**2-g*dgqgf2))
    qb = qb - g*(2*dgqgb/(2*dgqgb**2-g*dgqgb2))
    print iter,g 
    iter+=1 
    
    
  x=exp((Vgb-Vtb-vt*log(qb**2)-rb*qb)/(-2*vt))
  y=exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/(-2*vt))
  if abs(y) > 1:
    amb = arccosh(y)
  else:
    amb = arccos(y)
  if abs(x) > 1:
    apb = arccosh(x)
  else:
    apb = arccos(x)
  aguess = (amb+apb)/2  
  bguess = (-amb+apb)/2     
  """
  g=1
  iter=0
  #print "solving abcase2 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
  while iter<2:

    #if norm(g)<1e-10:
    #  break  
     
    g = (qf)**2-(qb)**2-eta*exp((Vgb-Vfbb-rb*qb)/vt)*(exp((Vgf-Vfbf-rf*qf)/vt-(Vgb-Vfbb-rb*qb)/vt)-1)
    dgqgb = -eta*rb*exp((-Vfbb + Vgb - qb*rb)/vt)/vt - 2*qb
    dgqgb2 = eta*rb**2*exp((-Vfbb + Vgb - qb*rb)/vt)/vt**2 - 2
    
    dgqgf = eta*rf*exp((-Vfbf + Vgf - qf*rf)/vt)/vt + 2*qf
    dgqgf2 = -eta*rf**2*exp((-Vfbf + Vgf - qf*rf)/vt)/vt**2 + 2
    qf = qf - g*(2*dgqgf/(2*dgqgf**2-g*dgqgf2))
    #qb = qb - g*(2*dgqgb/(2*dgqgb**2-g*dgqgb2))
    print iter,g 
    iter+=1 
  """
  return qf*(4*ech*vt/(Tsi)),qb*(4*ech*vt/(Tsi)),abs(aguess*0.9),abs(bguess*0.85) 
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
phi_gatef 	    = 4.0
phi_gateb 	    = 5.0


Vgb_array = linspace(-1,3,500)
Vgf_array = linspace(-1,-1,1)

for Vgf in Vgf_array:
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
  Cins1 = []
  deltaVth1 = []
  deltaphi = []
  Qf = []
  Qb = []
  qfguess = []
  qbguess = []
  aguess = []
  bguess = []
  qfNRM = []
  qbNRM = []  
  gvalue = []
  iterab = []
  for Vgb in Vgb_array:
    #Vgf=Vgb
    #if abs((Vgf-phi_gatef)-(Vgb-phi_gateb))<1e-10:
      #phi_gateb+=1e-5
    qmbaux,qmfaux,Cins1aux,deltaVth1aux,deltaphiaux,qfguessaux,qbguessaux,aguessaux,bguessaux =  unified_IMG_model(Vgf,Vgb,0,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gatef,phi_gateb,Tsi,tins,tinsbox,Hfin)  
    qfNRMaux,qbNRMaux,aguessaux,bguessaux = qfqb(qfguessaux,qbguessaux,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
      
    if (Vgf-phi_gatef)>(Vgb-phi_gateb):
      aaux,baux,phisfaux,phisbaux,Qtaux,Vcriticalaux,saux,caseaux,Qfaux,Qbaux,iterabaux = ab(aguessaux,bguessaux ,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    else:
      aaux,baux,phisbaux,phisfaux,Qtaux,Vcriticalaux,saux,caseaux,Qfaux,Qbaux,iterabaux = ab(aguessaux,bguessaux ,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg)
      
    
    iterab.append(iterabaux)
    a.append(aaux)
    b.append(baux)
    phisf.append(phisfaux)
    phisb.append(phisbaux) 
    Qt.append(Qtaux) 
    Vcritical.append(Vcriticalaux)
    s.append(saux)
    case.append(caseaux)
    qmb.append(qmbaux)
    qmf.append(qmfaux)
    qmt.append(qmfaux+qmbaux)
    Cins1.append(Cins1aux)    
    deltaVth1.append(deltaVth1aux)
    deltaphi.append(deltaphiaux)
    Qf.append(Qfaux)
    Qb.append(Qbaux)
    qfguess.append(qfguessaux)
    qbguess.append(qbguessaux)
    aguess.append(aguessaux)
    bguess.append(bguessaux)
    qfNRM.append(qfNRMaux)
    qbNRM.append(qbNRMaux) 
    #gvalue.append(gvalueaux)
 
  pylab.figure(8)
  pylab.plot(Vgb_array ,case,'o')  
 
  pylab.figure(7)
  pylab.plot(Vgb_array ,iterab,'k') 
  
 
  pylab.figure(6)
  pylab.plot(Vgb_array ,a,'r') 
  pylab.plot(Vgb_array ,b,'b')
  pylab.plot(Vgb_array ,aguess,'r--') 
  pylab.plot(Vgb_array ,bguess,'b--') 
  ax = pylab.gca()
  #ax.set_yscale('log')
   
  pylab.figure(3)
  pylab.plot(Vgb_array ,a,'r') 
  pylab.plot(Vgb_array ,b,'b') 
  pylab.plot(Vgb_array ,aguess,'r--') 
  pylab.plot(Vgb_array ,bguess,'b--')   
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("alpha (red) and beta (blue)", fontsize=18)
    
  pylab.figure(4)
  pylab.plot(Vgb_array ,Qt,'k-')
  pylab.plot(Vgb_array ,qmf,'b--')  
  pylab.plot(Vgb_array ,qmb,'k--') 
  pylab.plot(Vgb_array ,qmt,'k')  
  pylab.plot(Vgb_array ,Qf,'r-')
  pylab.plot(Vgb_array ,Qb,'b-')
  pylab.plot(Vgb_array ,qbguess,'y-')  
  pylab.plot(Vgb_array ,qfguess,'g-')
  pylab.plot(Vgb_array ,qfNRM,'c--') 
  pylab.plot(Vgb_array ,qbNRM,'m--')   
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("Qt", fontsize=18)
  
  pylab.figure(5)
  pylab.plot(Vgb_array ,Qt,'r') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("Qt", fontsize=18)
  pylab.plot(Vgb_array ,qmf,'b--')  
  pylab.plot(Vgb_array ,qmb,'k--') 
  pylab.plot(Vgb_array ,qmt,'k')  
  ax = pylab.gca()
  ax.set_yscale('log')
  
  pylab.figure(1)
  pylab.plot(Vgb_array ,Cins1,'r') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("Cins1", fontsize=18)

  pylab.figure(1)
  pylab.plot(Vgb_array ,deltaVth1,'r') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("deltaVth1", fontsize=18)


  pylab.figure(2)
  pylab.plot(Vgb_array ,deltaphi,'r') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("deltaVth1", fontsize=18)

pylab.show()  
