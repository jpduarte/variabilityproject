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
X1 supply gate 0 bulk nmos1 L=1.95539006831e-08 HFIN=4.15865407345e-08 TFIN_TOP=7.33580765015e-09 TFIN_BASE=1.5023328429e-08 EOT=8.05734536537e-10 NBODY=5.98998082225e+24 NFIN = 1 PHIG = 4.19294359378 RSHS = 355.223509859 RSHD = 275.033500903
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
