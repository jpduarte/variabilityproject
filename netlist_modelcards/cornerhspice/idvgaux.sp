*Sample netlist for BSIM-MG
*Id-Vg Characteristics for NMOS (T = 27 C)

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 25

.hdl "/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va"
.include "modelcardaux.nmos"

* --- Voltage Sources ---
vds supply  0 dc=0.05
vgs gate  0 dc=1
vbs bulk  0 dc=0
* --- Transistor ---
X1 supply gate 0 bulk nmos1 L=2.2e-08 HFIN=4.62e-08 TFIN_TOP=8.36e-09 TFIN_BASE=1.672e-08 EOT=8.8e-10 NBODY=6.6e+24 NFIN=1 
*+HFIN=42n 
*+TFIN_TOP = 7.6n 
*+TFIN_BASE = 15.20n
*+EOT =800.0p
*+NBODY =6.000E+24
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
