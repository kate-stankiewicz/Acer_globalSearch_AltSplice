#!/bin/bash
#SBATCH -J sa_q_s2
#SBATCH -o q_sp2_"%J".out
#SBATCH -e q_sp2_"%J".out
#SBATCH -c 40

# Call events in SplAdder
#run the SplAdder script 
srun ${PWD}/run_spladder_EC.sh

