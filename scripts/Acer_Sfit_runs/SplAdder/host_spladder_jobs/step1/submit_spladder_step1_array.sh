#!/bin/bash
#SBATCH --array=1-125
#SBATCH -J sa_arr
#SBATCH -o ar_sp1_"%A"."%a".out
#SBATCH -e ar_sp1_"%A"."%a".out
#SBATCH -c 2

# Job array for building splice graphs for each sample separately using SplAdder

#make sure conda base env is activated
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

# sample bam files from STAR (must be indexed first)
dir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/SplAdder/host_spladder_res'
samplesheet='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/SplAdder/host_spladder_jobs/host_bams.txt'
samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`

#make results directories
mkdir -p ${dir}/array_spladder_out
mkdir -p ${dir}/array_spladder_out/spladder
mkdir -p ${dir}/array_spladder_out/tmp

#check which sample is running
echo ${samp}

#run the SplAdder script for each separately
srun ${PWD}/run_test_spladder_array.sh ${samp}
