# General description
* These are the simulation codes of manuscript "Structural Scheduling of Transient Control under Energy Storage Systems by Sparse-Promoting Reinforcement Learning". The simulation requires the installation of Mat-power toolbox in MATLAB.
* The source codes simulate the transient process of power grids with energy storage systems under a contingency, the proposed method of the manuscript is impelemented in the source code.

# Simulation environment requirements
1. Matlab, version 2019a.
2. Matpower 7.0, which can be downloaded from URL:https://matpower.org/.
3. CVX semi-definite programming toolbox, which can be downloadeded from URL:http://cvxr.com/cvx/. 

# Start the simulation
All of the simulation program entry file names are started with "main" in each folder, the user can start the simulation by run the simulation program entry files in matlab.


# Folder descriptions
* The folder "Simulation without control" includes the source code files for the power system transient simulation with a contingency happen, without control.
* The folder "Simulation under proposed method" includes the source code files for the power system transient simulation with our proposed control in the paper after a contingency, and the control is activated at 0.4 s.
* The folder "Control under various gamma" includes the source code files for the power system transient simulation with our proposed control method in the various values of gamma.
* The folder "Control under various latency" includes  the source code files for the power system transient simulation with our proposed control method in the various values of communication latencies.
* The folder "Comparison of conventional works" includes the source code files for the power system transient simulation with three conventional control method for comparisions.
* The folder "Comparison of RL methods" contains the source code files of three comparison RL methods, including A2C, A3C and DDPG. Among them, A3C needs to use MATLAB parallel computing.
