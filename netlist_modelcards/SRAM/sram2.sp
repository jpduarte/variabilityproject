*Sample netlist for inverter

.option abstol=1e-6 reltol=1e-6 post ingold
.temp 25

.hdl "/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va"
.include "modelcarsram.nmos"
.include "modelcarsram.pmos"
.include "paramsramaux.include"

* --- Voltage Sources ---
vdd supply  0 dc=0.86
vl vl1  0 dc=0.86

* --- inverter ---
X1 vl1 vl2 0 0 nmos1 DEVTYPE=1 L=Lparam1 HFIN=HFINparam1 TFIN_TOP=TFIN_TOPparam1 TFIN_BASE=TFIN_BASEparam1 EOT=EOTparam1 NBODY=NBODYparam1 NFIN = NFINparam1 PHIG = PHIGparam1 RSHS = RSHSparam1 RSHD = RSHDparam1
X2 vl1 vl2 supply supply pmos1 DEVTYPE=0 L=Lparam2 HFIN=HFINparam2 TFIN_TOP=TFIN_TOPparam2 TFIN_BASE=TFIN_BASEparam2 EOT=EOTparam2 NBODY=NBODYparam2 NFIN = NFINparam2 PHIG = PHIGparam2 RSHS = RSHSparam2 RSHD = RSHDparam2

X3 vl2 vl1 0 0 nmos1 DEVTYPE=1 L=Lparam3 HFIN=HFINparam3 TFIN_TOP=TFIN_TOPparam3 TFIN_BASE=TFIN_BASEparam3 EOT=EOTparam3 NBODY=NBODYparam3 NFIN = NFINparam3 PHIG = PHIGparam3 RSHS = RSHSparam3 RSHD = RSHDparam3
X4 vl2 vl1 supply supply pmos1 DEVTYPE=0 L=Lparam4 HFIN=HFINparam4 TFIN_TOP=TFIN_TOPparam4 TFIN_BASE=TFIN_BASEparam4 EOT=EOTparam4 NBODY=NBODYparam4 NFIN = NFINparam4 PHIG = PHIGparam4 RSHS = RSHSparam4 RSHD = RSHDparam4

X5 vl2 supply supply 0 nmos1 DEVTYPE=1 L=Lparam5 HFIN=HFINparam5 TFIN_TOP=TFIN_TOPparam5 TFIN_BASE=TFIN_BASEparam5 EOT=EOTparam5 NBODY=NBODYparam5 NFIN = NFINparam5 PHIG = PHIGparam5 RSHS = RSHSparam5 RSHD = RSHDparam5
X6 vl1 supply supply 0 nmos1 DEVTYPE=1 L=Lparam6 HFIN=HFINparam6 TFIN_TOP=TFIN_TOPparam6 TFIN_BASE=TFIN_BASEparam6 EOT=EOTparam6 NBODY=NBODYparam6 NFIN = NFINparam6 PHIG = PHIGparam6 RSHS = RSHSparam6 RSHD = RSHDparam6
* --- DC Analysis ---
.dc vl 0.0 0.86 0.01 
.print dc v(vl2)
.end
