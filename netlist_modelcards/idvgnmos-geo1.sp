*Sample netlist for BSIM-MG
*Id-Vg Characteristics for NMOS (T = 27 C)

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 27

.param tox_mc=agauss(0,1,1)
.param tox_n=0.8n

.param lint_n=0
.param lint_mc=agauss(0,1,1)

.param wint_n=0
.param wint_mc=agauss(0,1,1)

.param tfin_n=7.6n
.param tfin_mc=agauss(0,1,1)

.param hfin_n=42n
.param hfin_mc=agauss(0,1,1)

.param nbody_n=6.0e24
.param nbody_mc=agauss(0,1,1)


.hdl "../model_code/code_109beta/bsimcmg.va"
.include "modelcard-109-geo1.nmos"

* --- Voltage Sources ---
vds supply  0 dc=0.05
vgs gate  0 dc=1
vbs bulk  0 dc=0
vds2 supply2  0 dc=0.86
* --- Transistor ---
X1 supply gate 0 bulk nmos1 L=20n NFIN=1 
X2 supply2 gate 0 bulk nmos1 L=20n NFIN=1 
* --- DC Analysis ---
.dc vgs -0.5 1.0 0.01 monte=50
.probe dc ids1=par'-i(vds)'
.probe dc par'-i(vbs)'
.print dc i(X1.d)


*.measure dc Vth_lin when i(X1.d)=300n
*.measure dc Vth_sat when i(X2.d)=300n

*.measure dc Id_sat i(X2.d) when v(gate)=0.86

.end
