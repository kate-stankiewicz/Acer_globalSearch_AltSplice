#!/bin/bash
#SBATCH --array=1-51
#SBATCH -J count_mr
#SBATCH -o count_mr_"%A"."%a".out
#SBATCH -e count_mr_"%A"."%a".out
#SBATCH -c 2

# Job array for counting total mapped reads for the split host / sym bams

# make sure conda base env is activated
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

# parse the array
samplesheet='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/Acer_Sfit_runs/RNAseq_pipeline/count_mapped_reads/add_fail_samples_sym/fail_sample_species_list.txt'
samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`
species=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'`

# check what is running
echo ${samp} ${species}

# run the samtools script for each separately
srun ${PWD}/run_samtools_count_mapped_reads.sh ${samp} ${species}
