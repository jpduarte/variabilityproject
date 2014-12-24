Readme:

letiv5v3.py:
  - implement LETI model using Newton Method
derivativeLETI.py:
  - computes the derivative for the paper LETI

letiv5v3.py:
  - analysis of LETI mode, there seems to show oscillation due cos, and sin functions
  
letiv5v5.py
  - start to write derivatives and functions for the case q~0  

letiv5v6.py
   - change the log(q^2/sin(q/2)) term to log(q*coth)-log(cos()), maybe better numerical stability
  
taurIMGv6.py:
  - solves Taur model for assymetric double gate
  - implementation of decouple model
  - implementation of initial guess for charge
  
taurIMGv3v2.py:
  - solves Taur model for assymetric double gate
  - plot assumtions of new model in function called "equ1"
  - plot first Initial Guess

taurIMGv3v3.py:  
  - funtion "qfqb" try to solve couple equation to obtain qf and qb, HOWEVER, still not able to find two equations which are independent and can be used to solve qf and qb
  
Ideas:
 1 in taurIMGv6.py, function qfqb: fix qb then solve for qf: DO NOT WORK, NM is not stable, no idea why
 2 in taurIMGv6.py, function qfqb: fix qf then solve for qb: work but not smoth
 3 assume to add zero field after equivalent gate voltage is sweep: DONOT WORK
 4 do 2 then do 1, does not work, seems that 1 is not stable
 5 use charges from decouple model to obtain an initial guess for alpha and beta, implemented in function "equ1" and file taurIMGv6.py, it is not stable, no idea why
 5 do 2 then apply 5, it is not stable, no idea why
 6 in letiv5v5, only get derivatives without taking cosh or sinh, doesnt work
 
Ideas to do:
 1 instead of 2 equations with qf, and qb, do 3 equation with alpha as paper
 2 do 1 in ideas, but with different initial guess for qf, maybe qb=qf as initial
 3 check LETI assumtions by comparing to numerical simulations
 4 change the log(q^2/sin(q/2)) term to log(q*coth)-log(cos()), maybe better numerical stability
