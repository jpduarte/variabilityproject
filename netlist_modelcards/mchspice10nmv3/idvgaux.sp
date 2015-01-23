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
X1 supply gate 0 bulk nmos1 L=1.76270787271e-08 HFIN=4.12576219394e-08 TFIN_TOP=6.41747784641e-09 TFIN_BASE=7.33971979645e-09 EOT=7.46007769778e-10 NBODY=7.02211475076e+24 NFIN = 1 PHIG = 4.20112557786 RSHS = 376.924308447 RSHD = 240.327961141
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
