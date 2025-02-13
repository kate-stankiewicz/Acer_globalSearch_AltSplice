#!/bin/bash
#SBATCH --array=1-176
#SBATCH -J htseq
#SBATCH -o htseq_"%A"."%a".out
#SBATCH -e htseq_"%A"."%a".out
#SBATCH -c 2

# load the STAR env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

# parse the array
list="/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/full_sample_list.txt"
file=`sed -n "$SLURM_ARRAY_TASK_ID"p $list |  awk '{print $1}'`


# run samtools sort
srun ${PWD}/run_htseq_count_no.sh ${file}
