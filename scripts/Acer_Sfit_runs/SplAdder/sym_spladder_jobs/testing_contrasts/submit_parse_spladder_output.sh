#!/bin/bash
#SBATCH --array=1-30
#SBATCH -J parse_arr
#SBATCH -o parse_"%A"."%a".out
#SBATCH -e parse_"%A"."%a".out
#SBATCH --mem-per-cpu=1200
#SBATCH -c 1

# Job array for parsing spladder statistical testing

# run each contrast separately in array
cont_sheet="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/SplAdder/array_jobs/host_spladder_jobs/testing_contrasts/all_contrasts.txt"
var1=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $1}'`
var2=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $2}'`
const=`sed -n "$SLURM_ARRAY_TASK_ID"p $cont_sheet |  awk '{print $3}'`

#check which contrast is running
echo ${var1}vs${var2}_${const}

#run the R script for each separately
Rscript /proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/SplAdder/array_jobs/host_spladder_jobs/testing_contrasts/run_parse_spladder_output.R ${var1} ${var2} ${const}

