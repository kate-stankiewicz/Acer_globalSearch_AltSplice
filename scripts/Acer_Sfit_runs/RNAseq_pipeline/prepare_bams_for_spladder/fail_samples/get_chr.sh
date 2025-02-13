#!/bin/bash

# script to get list of chrs for the host and the sym

# load the env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

# set the dir
dir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample/R382/results_STAR_Salmon'

# pull out the host chrs
samtools idxstats ${dir}/R382R382_Aligned_sorted.out.bam | cut -f 1 | grep "Acer" > host_chr.txt

# pull out the sym chrs
samtools idxstats ${dir}/R382R382_Aligned_sorted.out.bam | cut -f 1 | grep "Sfit" > sym_chr.txt
