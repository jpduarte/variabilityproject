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
X1 supply gate 0 bulk nmos1 L=1.9254690178e-08 HFIN=4.10294845459e-08 TFIN_TOP=7.58114688727e-09 TFIN_BASE=8.57029908089e-09 EOT=7.43684544493e-10 NBODY=6.20725451539e+24 NFIN = 1 PHIG = 4.22088553838 RSHS = 289.655280815 RSHD = 308.007327318
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
