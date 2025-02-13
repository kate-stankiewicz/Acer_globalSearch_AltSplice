#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

DIRS='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample'

#index host separated bam files
samtools sort -n -@ 8 ${DIRS}/${1}/results_STAR_Salmon/host_${1}_sorted.bam > ${DIRS}/${1}/results_STAR_Salmon/host_${1}_NameSorted.out.bam
