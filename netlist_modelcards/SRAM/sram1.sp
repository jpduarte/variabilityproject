*Sample netlist for inverter

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 25

.hdl "/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va"
.include "modelcarsram.nmos"
.include "modelcarsram.pmos"

* --- Voltage Sources ---
vdd supply  0 dc=0.86
vl vl1  0 dc=0.86

* --- inverter ---
X1 vl1 vl2 0 0 nmos1 L=20e-9 NFIN=1 DEVTYPE=1
X2 vl1 vl2 supply supply pmos1 L=50e-9 NFIN=1 DEVTYPE=0

X3 vl2 vl1 0 0 nmos1 L=20e-9 NFIN=1 DEVTYPE=1
X4 vl2 vl1 supply supply pmos1 L=50e-9 NFIN=1 DEVTYPE=0

X5 vl2 supply supply 0 nmos1 L=20e-9 NFIN=1 DEVTYPE=1
X6 vl1 supply supply 0 nmos1 L=20e-9 NFIN=1 DEVTYPE=1
* --- DC Analysis ---
.dc vl 0.0 0.86 0.01 
.print dc v(vl2)
.end
