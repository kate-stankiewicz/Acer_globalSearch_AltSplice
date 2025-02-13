#!/bin/bash

# get count of mapped reads for each bam for sym and host separately

# set env dir
DIR=/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample

# run samtools
samtools view -c -F 260 ${DIR}/${1}/results_STAR_Salmon/${2}_${1}_sorted.bam > ${2}_${1}_total_mapped_reads_count.txt
