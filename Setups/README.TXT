Jul 30, 2023

This folder contains the setup files to run the replica exchange method (REM) to sample minimal-curvature ring melts in duboids of different sizes.
There is one directory per system size, e.g. L4x4x3, L4x4x4 etc., containing teh required input files.

To perform the REM-based sampling it suffices to enter the relevant directory and run the executable, e.g.

cd L4x4x3
../../Replica_exchange_source_code/replica_exchange_sampling


The run will create two output data files per each replica, containing the time series of the energies (equal to the curvature or total number of corners) and of the configurations.
The latter are encoded as a one-dimensional string of 0/1 values, named "spins" in the code. These "spins" correspond to the occupied (1) or empty (0) state of a lattice edge. 

With default parameters, the lowest-temperature replica, which acts as collector for the sought solutions has temperature equal to 0.32, and the two associated output files are:

output_energy_T3.200000e-01.dat
output_spin_configs_T3.200000e-01.dat


For each system, the ground state energy (minimal curvature) is provided in the file named "Target_n_corners.dat".


