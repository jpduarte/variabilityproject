from numpy import exp,tan,tanh,log,abs,cosh,sqrt,sinh,sin,cos,linspace,array,pi,matrix,dot, arccos, arccosh
from numpy.linalg import norm,solve,inv,det

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
      
  return qmb,qmf#,Cins1/Cins,deltaVth1,deltaphi,qfguess,qbguess,aguess,bguess
  

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
