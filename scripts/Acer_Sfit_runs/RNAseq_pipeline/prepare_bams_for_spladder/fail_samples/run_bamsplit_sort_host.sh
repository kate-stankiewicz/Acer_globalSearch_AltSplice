#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

DIRS='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample'

#split host from bam
samtools view -@ 8 -b -h ${DIRS}/${1}/results_STAR_Salmon/${1}R382_Aligned_sorted.out.bam $(cat ${PWD}/host_chr.txt) | samtools sort -@ 8 - > ${DIRS}/${1}/results_STAR_Salmon/host_${1}_sorted.bam
