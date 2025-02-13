#!/bin/bash
#SBATCH --array=1-56
#SBATCH -J test_arr
#SBATCH -o in_test_"%A"."%a".out
#SBATCH -e in_test_"%A"."%a".out
#SBATCH -c 5

# Job array for spladder statistical testing

# run each contrast separately in array
cont_sheet="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/SplAdder/sym_spladder_jobs/testing_contrasts/fail_pass_concat/contrast_files/all_contrasts.txt"
var1=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $1}'`
var2=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $2}'`

#check which contrast is running
echo ${var1}vs${var2}

#run the SplAdder script for each separately
srun ${PWD}/run_spladder_test_contrasts.sh ${var1} ${var2}
