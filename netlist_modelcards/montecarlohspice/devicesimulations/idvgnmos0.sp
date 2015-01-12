*Sample netlist for BSIM-MG
*Id-Vg Characteristics for NMOS (T = 27 C)

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 27

.param vdsvalue = 0.05

.hdl "../../../model_code/code_109beta/bsimcmg.va"
.include "modelcard-109-geo1.nmos"

* --- Voltage Sources ---
vds supply  0 dc=vdsvalue
vgs gate  0 dc=1
vbs bulk  0 dc=0
* --- Transistor ---
X1 supply gate 0 bulk nmos1 L=20n NFIN=1 
* --- DC Analysis ---
.dc vgs 0.0 1.0 0.01 monte=3
*.print dc v(X1.d)
.print dc i(X1.d)

.alter 
.param vdsvalue = 0.86



.end
