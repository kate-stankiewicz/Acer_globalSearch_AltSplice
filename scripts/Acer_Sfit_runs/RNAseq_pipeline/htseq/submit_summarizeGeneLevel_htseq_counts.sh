#!/bin/bash
#SBATCH -J sum_gene
#SBATCH -o sum_gene_"%A"."%a".out
#SBATCH -e sum_gene_"%A"."%a".out
#SBATCH -c 2

#run the R script on each sample counts file
for samp in $(cat /proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/htseq_results/host/sample_list_fixed.txt); do
Rscript ${PWD}/summarizeGeneLevel_htseq_counts.R ${samp}
done
