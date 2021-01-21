# File descriptions
* Files main_directDamping.m, main_IntegralDamping.m and main_sparsePromote.m are source code files of the simulation program entries of methods of a decentralized droop control strategy with an energy function based ESS allocation scheme,CPI-based control strategy with a stability index based allocation scheme and a centralized control strategy with a sparse-promoting ESS allocation scheme respectively.
* File nonlinear_dynamic.m is our system model.
* File ReducedY.m,is source code file of functions of radial basis function of RBFNN in methods A2C, A3C and DDPG respectively.
* File sat.m is the source code file of saturation function for energy storage systems
* File bas.m is the source code file of raidual basis function
* File CalculateEnergy.m is the source code file of energy function.
* File CalculatingL.m is the source code file of calculating A matrix for sparse promote method.
* File ReducedY.m is the source code file of calculating the admittance matrix.
# Note
the running of sparse promote methods need the cvx toolbox for solving semi-definite program problems, which can be downloaded from URL:http://cvxr.com/cvx/.
