#!/bin/bash
#SBATCH -J test_spl
#SBATCH -o cross_temp_"%A".out
#SBATCH -e cross_temp_"%A".out
#SBATCH -c 5
#SBATCH -w baliga2

# Job script for spladder statistical testing

#run the SplAdder script 
srun ${PWD}/run_spladder_test_contrasts_cross_temp.sh
