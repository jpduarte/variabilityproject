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
X1 supply gate 0 bulk nmos1 L=1.9931040809e-08 HFIN=4.45567090069e-08 TFIN_TOP=7.72901568089e-09 TFIN_BASE=1.51846145703e-08 EOT=8.12585045709e-10 NBODY=5.80857313999e+24 NFIN = 1 PHIG = 4.19475966434 RSHS = 264.088966979 RSHD = 228.715803361
* --- DC Analysis ---
.dc vgs 0.0 0.86 0.0175510204082 vds 0.05 0.86 0.81
.print dc i(X1.d)
.end
