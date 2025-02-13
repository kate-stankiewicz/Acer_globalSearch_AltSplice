#!/bin/bash
#SBATCH --array=1-125
#SBATCH -J sa_q_arr
#SBATCH -o q_sp1_"%A"."%a".out
#SBATCH -e q_sp1_"%A"."%a".out
#SBATCH -c 2

# Job array for building splice graphs for each sample separately using SplAdder

# sample bam files from STAR (must be indexed first)
samplesheet='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/SplAdder/host_spladder_jobs/host_bams.txt'
samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`

#check which sample is running
echo ${samp}

#run the SplAdder script for each separately
srun ${PWD}/run_spladder_quant_step1.sh ${samp}

