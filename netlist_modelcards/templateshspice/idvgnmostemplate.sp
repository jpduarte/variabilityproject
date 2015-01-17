*Sample netlist for BSIM-MG
*Id-Vg Characteristics for NMOS (T = 27 C)

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 25

.hdl pathmodelverilog
.include pathmodelcard

* --- Voltage Sources ---
vds supply  0 dc=0.05
vgs gate  0 dc=1
vbs bulk  0 dc=0
* --- Transistor ---
X1 supply gate 0 bulk nmos1 L=Lparam NFIN=NFINparam DEVTYPE = DEVTYPEparam
* --- DC Analysis ---
.dc vgs vgsi vgsf vgsdelta vds vdsi vdsf vdsdelta
.print dc i(X1.d)
.end
