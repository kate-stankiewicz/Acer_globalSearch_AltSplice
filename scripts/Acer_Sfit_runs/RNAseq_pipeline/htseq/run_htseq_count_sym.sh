#!/bin/bash

# load the STAR env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate STAR

#run htseq count on all results
dir=/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample
gff=/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/genomics/Sfitti_EVM.all.ID.gff3

htseq-count -s no -t exon -i Parent -r pos ${dir}/${1}/results_STAR_Salmon/sym_${1}_sorted.bam ${gff} > ${dir}/${1}/results_STAR_Salmon/htseqcounts/no_sym_${1}_htseq_count.txt
