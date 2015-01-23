*Sample netlist for inverter

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 25

.hdl "/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va"
.include "modelcarsram.nmos"

* --- Voltage Sources ---
vdd supply  0 dc=0.86
vinput vin  0 dc=1
* --- inverter ---
X1 vout vin 0 0 nmos1 L=50e-9 NFIN=1 DEVTYPE=1
X2 vout vin supply supply nmos1 L=50e-9 NFIN=1 DEVTYPE=0
* --- DC Analysis ---
.dc vinput 0.0 0.86 0.0175510204082 
.print dc v(X1.d)
.end
