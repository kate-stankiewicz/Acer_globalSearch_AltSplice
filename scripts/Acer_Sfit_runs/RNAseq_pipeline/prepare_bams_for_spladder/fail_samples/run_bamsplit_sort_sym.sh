#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

DIRS='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/RNASeq_pipeline/results_per_sample'

#split sym from bam
samtools view -@ 8 -b -h ${DIRS}/${1}/results_STAR_Salmon/${1}R382_Aligned_sorted.out.bam -o /dev/null -L ${PWD}/host.bed -U - | samtools sort -@ 8 - > ${DIRS}/${1}/results_STAR_Salmon/sym_${1}_sorted.bam


# had to split this a different way because the number of contigs in the sym reference were so so many that it generates an "Argument list too long" error
