#!/bin/bash
#SBATCH -J sa_q_s2
#SBATCH -o q_sp2_"%J".out
#SBATCH -e q_sp2_"%J".out
#SBATCH -c 40

# Step 2 of SplAdder quantification, collect all individual quantification into one joint db
#run the SplAdder script 
srun ${PWD}/run_spladder_quant_step2.sh

