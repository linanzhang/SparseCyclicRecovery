# SparseCyclicRecovery

Target equation: the Lorenz 96 equation
      udot_{i} = ( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + F

Auxiliary functions:
  - Dictionary.m (The Candidate Functions)
  - DouglasRachford.m (Optimization Routine) 
  - leg2mon.m (Legendre to Monomial Transform)
  - lorenz96.m (ODE to Test)

Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang

Date: May 9, 2018
