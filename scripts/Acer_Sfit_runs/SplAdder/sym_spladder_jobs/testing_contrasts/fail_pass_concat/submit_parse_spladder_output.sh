#!/bin/bash
#SBATCH --array=1-56
#SBATCH -J parse_arr
#SBATCH -o parse_"%A"."%a".out
#SBATCH -e parse_"%A"."%a".out
#SBATCH -c 2

# Job array for parsing spladder statistical testing

# run each contrast separately in array
cont_sheet="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/SplAdder/sym_spladder_jobs/testing_contrasts/fail_pass_concat/contrast_files/all_contrasts.txt"
var1=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $1}'`
var2=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $2}'`

#check which contrast is running
echo ${var1}vs${var2}

#run the R script for each separately
Rscript ${PWD}/run_parse_spladder_output.R ${var1} ${var2}

