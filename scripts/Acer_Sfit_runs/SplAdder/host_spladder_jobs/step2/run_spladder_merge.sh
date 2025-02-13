#!/bin/bash
# step 2 for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#Merge the splice graphs for each sample 
spladder build -o ${workdir}/results/SplAdder/host_spladder_res/array_spladder_out -a ${workdir}/genomics/Acropora_cervicornis.ID.gff3 -b ${workdir}/scripts/SplAdder/host_spladder_jobs/host_bams.txt --merge-strat merge_graphs --no-extract-ase --no-quantify-graph --parallel 18 -v


