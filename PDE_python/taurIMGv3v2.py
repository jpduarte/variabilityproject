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
    #print s,f,iter
  Vcritical = Vfbf+C1+2*vt*log(2/(s-1))+2*vt*r*(2/(s-1))
  #print "Vcritical: "+ str(Vcritical)
  return Vcritical,s

def ab(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vcritical,s = vcritical(Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
  
  if Vgf < Vcritical:
    a,b,phisf,phisb = abcase1(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    cotaminusb = b*coth(a-b)
    cotaplusb = -b*coth(a+b)
    cosminus = cosh((a-b))**2
    cosplus = cosh((a+b))**2
    Qt = 4*ech*vt/(Tsi)*b*(coth(a-b)-coth(a+b)) 
    case = 1
    
  else:
    a,b,phisf,phisb = abcase2(aguess,bguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    cotaminusb = b*cot((a-b))
    cotaplusb = -b*cot(a+b)
    cosminus = cos((a-b))**2
    cosplus = cos((a+b))**2    
    Qt = 4*ech*vt/(Tsi)*b*(cot(a-b)-cot(a+b))
    case = 2  
  return a,b,phisf,phisb,Qt,Vcritical,s,case,cotaminusb,cotaplusb,cosminus,cosplus
  
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
  #print "solving abcase1 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
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
    #print det(J) 
    if  abs(det(J)  )<1e-10:
      a=aguess*10
      b=bguess*10 
    else:
      s_m = s_m - solve(J,Y) 
    iter+=1
    #print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sinh((-b*Tsi/Tsi)+a))
  #Qfguess = (Vg
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
  #print "solving abcase2 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
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
    #print det(J)
    if  abs(det(J)  )<1e-10:
      a=aguess*10
      b=bguess*10 
    else:
      s_m = s_m - solve(J,Y)  
    iter+=1
    #print iter,a,b,norm(Y)
    
  phisf = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((b*Tsi/Tsi)+a))
  phisb = -vt*2*log((Tsi/(2*b))*sqrt((q*ni)/(2*ech*vt))*sin((-b*Tsi/Tsi)+a))    
  return a,b,phisf,phisb

def Qt_guess(Vcritical,s,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Cins = (Hfin)*eins/tins
  Ach=Hfin*(Tsi/2)
  Weff = Hfin

  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)#+(Vgf-Vgb)
  Vthf     = Vfbf+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))
  Vthb     = Vfbb+vt*log(Cins*vt/(q*ni**2*2*Ach/Nch))
  #print "Vthf %e, Vthb: %e" % (Vthf,Vthb)
  
  Qt = 2*(eins/tins)*(Vgf-0.5*(Vthb+Vthf))#+(eins/tinsbox)*(Vgf-Vthb)
  x = sqrt((log(1+exp(3*(Vgf-Vcritical-2*vt)/vt))/3+1)**2-1)
  #alpha = 1/(0.7*2/pi+1/x)#(4*ech*vt/Tsi)*(pi/2)*(s/(s-1))/Qt
  alpha = sqrt(1/(2/(pi*0.9)+1/x)**2)
  alpha2 = sqrt((log(1+exp(-3*(Vgf-Vcritical+vt)/vt))/(15*2)+1)**2-1)
  return Qt,alpha+alpha2

def equ1(a,b,Vcritical,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  xi = q*ni/(2*ech*vt)
  rf = tins*ech*2*vt*2/(eins*Tsi)
  rb = tins*ech*2*vt*2/(eins*Tsi)
  Vtf = Vfbf-2*vt*log((Tsi/2)*sqrt(xi))
  Vtb = Vfbb-2*vt*log((Tsi/2)*sqrt(xi))
  eta = 2*q*ni*vt/ech/(4*vt/Tsi)**2
  if(Vgf<Vcritical):
    qf = b*coth(a-b)
    qb = -b*coth(a+b)
    cosaminusb = cosh(a-b) 
    cosplusb = cosh(a+b) 
  else:
    qf = b*cot(a-b)
    qb = -b*cot(a+b)
    cosaminusb = cos(a-b)
    cosplusb = cos(a+b) 
  
  aux1 = (exp((Vgb-Vtb-rb*qb)/(-vt))+qb**2)  
  #f = (qf**2*exp((Vgf-Vtf-rf*qf)/vt)-qf**2+1)*aux1-qb**2*(1-aux1)
  #f = Vgf-Vtf-vt*log(qf**2)-rf*qf+vt*log(cosaminusb**2)
  #f = Vgb-Vtb-vt*log(qb**2)-rf*qb+vt*log(cosplusb**2)
  qr = qb/qf
  cosaminusb = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/(-2*vt))
  cosplusb = exp((Vgb-Vtb-vt*log(qb**2)-rf*qb)/(-2*vt))
  oneovercosaminusbsqr = (qr**2*(1-cosplusb**2)/(cosplusb**2)+1)
  #f = Vgf-Vtf-vt*log(qf**2)-rf*qf-vt*log(oneovercosaminusbsqr)
  #f = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/vt) -(oneovercosaminusbsqr) #same as previous line, but applying exponential, lower errors
  #f = (qf)**2-(qb)**2-eta*exp((Vgb-Vfbb-rb*qb)/vt)*(exp((Vgf-Vfbf-rf*qf)/vt-(Vgb-Vfbb-rb*qb)/vt)-1.0)
  #f = (Vgf-Vgb)-(Vtf-Vtb)+vt*log(qr**2)-(rf*qf-rb*qb)-vt*log(qr**2*(1-cosplusb**2)+cosplusb**2)
  #f = exp(((Vgf-Vgb)-(Vtf-Vtb)-(rf*qf-rb*qb))/vt) -( (1-cosplusb**2)+cosplusb**2/qr**2) #Doestnt work well, its is same as previous equation
  #f = 1/cosaminusb**2-oneovercosaminusbsqr
  #f = (Vgf-Vgb)-(Vtf-Vtb)+vt*log(qr**2)-(rf*qf-rb*qb)+vt*log(cosaminusb**2/cosplusb**2)#this is like saying 1=1
  #baux = (qf-qb)/(cosaminusb/(sqrt(abs(1-cosaminusb**2)))+cosplusb/(sqrt(abs(1-cosplusb**2))))
  baux = (qf-qb)/(cosaminusb/(sqrt(abs(1-cosaminusb**2)))+cosplusb/(sqrt(abs(1-cosplusb**2))))
  f = Vgf - Vtf -rf*qf-vt*log(baux**2/abs(1-cosaminusb**2))
  return abs(f)
  
def qfqb(qfguess,qbguess,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg):
  Vfbf     = phi_gatef - phi_substrate -Eg/2-vt*log(Nch/ni)
  Vfbb     = phi_gateb - phi_substrate -Eg/2-vt*log(Nch/ni)
  xi = q*ni/(2*ech*vt)
  rf = tins*ech*2*vt*2/(eins*Tsi)
  rb = tins*ech*2*vt*2/(eins*Tsi)
  Vtf = Vfbf-2*vt*log((Tsi/2)*sqrt(xi))
  Vtb = Vfbb-2*vt*log((Tsi/2)*sqrt(xi))
  eta = 2*q*ni*vt/ech/(4*vt/Tsi)**2


  #Newtonmethod
  qvec=array([qfguess,qbguess])

  Y = array([1,1])
  iter=0
  #print "solving abcase2 for Vgf: %f and Vgb: %f"% (Vgf,Vgb)
  while iter<100:
    qf = qvec[0]
    qb = qvec[1]
    if norm(Y)<1e-12:
      break  
     
    qr = qb/qf
    cosplusb = exp((Vgb-Vtb-vt*log(qb**2)-rf*qb)/(-2*vt))
    oneovercosaminusbsqr = (qr**2*(1-cosplusb**2)/(cosplusb**2)+1)
    #f = Vgf-Vtf-vt*log(qf**2)-rf*qf-vt*log(oneovercosaminusbsqr)
    f = exp((Vgf-Vtf-vt*log(qf**2)-rf*qf)/vt) -(oneovercosaminusbsqr) #same as previous line, but applying exponential, lower errors
    g = (qf)**2-(qb)**2-eta*(exp((Vgf-Vfbf-rf*qf)/vt)-exp((Vgb-Vfbb-rb*qb)/vt))
    
    Y = array([f,g])
      
    f_dqf=2*qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**3 + (-rf - 2*vt/qf)*exp((Vgf - Vtf - qf*rf - vt*log(qf**2))/vt)/vt
    f_dqb=-qb**2*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*(-rf - 2*vt/qb)*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/(qf**2*vt) - qb**2*(-rf - 2*vt/qb)/(qf**2*vt) - 2*qb*(1 - exp(-(Vgb - Vtb - qb*rf - vt*log(qb**2))/vt))*exp((Vgb - Vtb - qb*rf - vt*log(qb**2))/vt)/qf**2
    
    g_dqf=eta*rf*exp((-Vfbf + Vgf - qf*rf)/vt)/vt + 2*qf
    g_dqb=-eta*rb*exp((-Vfbb + Vgb - qb*rb)/vt)/vt - 2*qb
    J = array([[f_dqf,f_dqb],[g_dqf,g_dqb]])
    #print det(J)
    if  abs(det(J)  )<1e-10:
      qf = qfguess*10
      qb = qbguess*10
    else:
      qvec = qvec - solve(J,Y)  
    iter+=1
    print iter,qf,qb,norm(Y)  
  
  
  return qf,qb 
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
tins            = 2e-7
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

Vgb_array = linspace(-0.5,2,500)
Vgf = -0.5

Tsi_array = linspace(10e-7,10e-7,1)

for Tsi in Tsi_array:
  a = []
  b = []
  phisf = []
  phisb = []
  Qt = []
  Vcritical = []
  s = []
  case = []
  aaux,baux = 1.02*0.9*pi/2,0.9*pi/2
  qfNRMaux,qbNRMaux = -20,-20
  qmb = []
  qmf = []
  qmt = []
  cotaminusb = []
  cotaplusb = []
  Qtguess = []
  aguess = []
  bguess = []
  sdep = []
  cosminus = []
  cosplus = []
  cosminqf = []
  equation1 = []
  qfNRM = []
  qbNRM = []
  for Vgb in Vgb_array:
    Vgf=Vgb
    if abs((Vgf-phi_gatef)-(Vgb-phi_gateb))<1e-10:
      phi_gateb+=1e-5
      
    if (Vgf-phi_gatef)>(Vgb-phi_gateb):
      aaux,baux,phisfaux,phisbaux,Qtaux,Vcriticalaux,saux,caseaux,cotaminusbaux,cotaplusbaux,cosminusaux,cosplusaux = ab(aaux,baux ,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
      Qtguessaux,aguessaux = Qt_guess(Vcriticalaux,saux,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
      equation1.append(equ1(aaux,baux,Vcriticalaux,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg))
      #qfNRMaux,qbNRMaux = qfqb(qfNRMaux,qbNRMaux,Vgf,Vgb,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
    else:
      aaux,baux,phisbaux,phisfaux,Qtaux,Vcriticalaux,saux,caseaux,cotaminusbaux,cotaplusbaux,cosminusaux,cosplusaux = ab(aaux,baux ,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg)
      Qtguessaux,aguessaux = Qt_guess(Vcriticalaux,saux,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg)
      equation1.append(equ1(aaux,baux,Vcriticalaux,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg))
      #qfNRMaux,qbNRMaux = qfqb(qfNRMaux,qbNRMaux,Vgb,Vgf,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gateb,phi_gatef,phi_substrate,Eg)
      
    qfNRM.append(qfNRMaux)
    qbNRM.append(qbNRMaux)
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
    cotaminusb.append(cotaminusbaux)
    cotaplusb.append( cotaplusbaux)
    Qtguess.append(Qtguessaux)
    aguess.append(aguessaux)
    bguess.append(aguessaux/saux)
    sdep.append(aaux/baux)
    cosminus.append(log(cosminusaux))
    cosplus.append(log(cosplusaux))
    cosminqf.append(cotaminusbaux-log(cosminusaux))
    
    
    
  ##########plots  
  pylab.figure(1)
  pylab.plot(Vgb_array ,a,'r') 
  pylab.plot(Vgb_array ,b,'b') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("alpha (red) and beta (blue)", fontsize=18)
  pylab.plot(Vgb_array ,b,'b') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.plot(Vgb_array ,aguess,'k--') 
  pylab.plot(Vgb_array ,bguess,'g--') 

  pylab.figure(2)
  pylab.plot(Vgb_array ,a,'r') 
  pylab.plot(Vgb_array ,b,'b') 
  pylab.xlabel("Vgb (V)", fontsize=18)
  pylab.ylabel("alpha (red) and beta (blue)", fontsize=18)
  ax = pylab.gca()
  ax.set_yscale('log')
  pylab.plot(Vgb_array ,aguess,'k--') 
  pylab.plot(Vgb_array ,bguess,'g--') 

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
  pylab.plot(Vgb_array ,Qtguess,'r--')


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
  pylab.plot(Vgb_array ,sdep,'k')

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


  pylab.figure(10)
  pylab.plot(Vgb_array ,cotaminusb,'r') 
  pylab.plot(Vgb_array ,cotaplusb,'b') 
 
  pylab.figure(10)
  pylab.plot(Vgb_array ,cosminus,'r--') 
  pylab.plot(Vgb_array ,cosplus,'b--') 
  pylab.plot(Vgb_array ,cosminqf,'k--')  
  
  pylab.figure(11)
  pylab.plot(Vgb_array ,equation1,'o') 
  ax = pylab.gca()
  ax.set_yscale('log')
  #print equation1


    


pylab.show()  


##print vcritical(0.2,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
##print ab(1.1*pi/9,pi/10,1.0,0,Tsi,q,vt,ni,Nch,ech,eins,tins,tinsbox,phi_gatef,phi_gateb,phi_substrate,Eg)
