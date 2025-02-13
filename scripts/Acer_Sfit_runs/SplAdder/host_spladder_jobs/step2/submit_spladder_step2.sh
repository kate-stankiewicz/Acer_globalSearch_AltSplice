#!/bin/bash
#SBATCH -J sa_s2
#SBATCH -o sp2_"%A".out
#SBATCH -e sp2_"%A".out
#SBATCH -c 18

# Step 2 of SplAdder build, merging splice graphs run in previous step
#run the SplAdder script 

srun ${PWD}/run_spladder_merge.sh

