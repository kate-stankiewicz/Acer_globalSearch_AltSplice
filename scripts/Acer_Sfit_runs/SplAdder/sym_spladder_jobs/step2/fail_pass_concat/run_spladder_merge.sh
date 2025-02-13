#!/bin/bash
# step 2 for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#Merge the splice graphs for each sample 
spladder build -o ${workdir}/results/SplAdder/sym_spladder_res/array_spladder_out2 -a ${workdir}/genomics/Sfitti_EVM.all.ID.gff3 -b ${workdir}/scripts/SplAdder/sym_spladder_jobs/sym_bams_all.txt --merge-strat merge_graphs --no-extract-ase --no-quantify-graph --parallel 40 -v


