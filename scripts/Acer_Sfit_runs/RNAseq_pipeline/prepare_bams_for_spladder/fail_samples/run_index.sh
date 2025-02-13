#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

DIRS='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/array_jobs/throttle/rerun_throttle'

#index host separated bam files
samtools index -@ 10 ${DIRS}/STAR_results/${1}/host_${1}.sorted.bam
