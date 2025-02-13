#!/bin/bash
#SBATCH --array=1-125
#SBATCH -J ind
#SBATCH -o ind_host."%A"."%a".out
#SBATCH -e ind_host."%A"."%a".out
#SBATCH -c 8

# Job array for separating host from sym in STAR alignment files
# load the STAR env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

# set variable for parsing SLURM array and sample list
samplesheet="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample/sample_list.txt"
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`

#check which sample is running
echo ${name}

#run the samtools script for each separately
srun ${PWD}/run_index_host.sh ${name}
