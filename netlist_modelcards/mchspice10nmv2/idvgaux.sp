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
X1 supply gate 0 bulk nmos1 L=1.75812335068e-08 HFIN=4.22920099157e-08 TFIN_TOP=3.95720581472e-09 TFIN_BASE=5.63975235978e-09 EOT=7.84379397537e-10 NBODY=5.96879441524e+24 NFIN = 1 PHIG = 4.18622863831 RSHS = 372.052646314 RSHD = 334.003918883
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
