#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

DIRS='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample'

#index host separated bam files
samtools sort -@ 8 ${DIRS}/${1}/results_STAR_Salmon/${1}_star_10_0.3_0.66_0_Aligned.out.bam > ${DIRS}/${1}/results_STAR_Salmon/${1}R382_Aligned_sorted.out.bam
