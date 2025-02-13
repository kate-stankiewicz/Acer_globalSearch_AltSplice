#!/bin/bash
#SBATCH --array=1-176
#SBATCH -J na_ind
#SBATCH -o na_ind."%A"."%a".out
#SBATCH -e na_ind."%A"."%a".out
#SBATCH -c 8

# Job array for sorting STAR bams
# load the STAR env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

# set variable for parsing SLURM array and sample list
samplesheet="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/full_sample_list.txt"
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`

#check which sample is running
echo ${name}

#run the samtools script for each separately
srun ${PWD}/run_nameSort_Index_host.sh ${name}
