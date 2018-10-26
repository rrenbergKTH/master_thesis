# master_thesis
The code used in my master thesis "Superconductor-Insulator Quantum PhaseTransitions in a Dissipative Environment". Follow this procedure to use the files:

1. Download all files into a folder.

2. Create one subfolder named data and one named bessel_list_files.

3. Run create_bessel_list.m. This creates lists of values of log(I_1(E_J*Delta_tau)/I_0(E_J*Delta_tau)) for different values of K and choices of Delta_tau. The function takes the following arguments: create_bessel_list(K_min, K_max, number_of_K_values, alpha_min, alpha_max, number_of_alpha_values, [list of lambda values], [list of Ec_times_beta], [list of factor values]).

A working example: create_bessel_list(0.05, 10, 200, 0.05, 1, 20,[10 0.1], [0.1 0.2 5 10 15], [1 3 5 15]).

4. Compile all .cpp files. Run the program with the following arguments: file_name Lx Ec_times_beta values of K -alpha alpha_values -lambda lambda_value -Nwarmup N_equil -Nprod N_prod -factor value_of_factor -remote 0. The remote flag should most proably be set to 0 (see code for its meaning).

A working example: data/file_name.txt 2 0.1 1.5 -alpha 0.95 1.0 -lambda 10 -Nwarmup 100000 -Nprod 100000 -factor 1 -remote 0. 

5. After the simulation is completed, plot your data with plotCppData.m with e.g. the following command: plotCppData('file_name.txt',2).
