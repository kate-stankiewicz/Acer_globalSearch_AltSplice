#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#Step one: generate the splice graph for each sample independently 
spladder build -o ${workdir}/results/SplAdder/sym_spladder_res/array_spladder_out -a ${workdir}/genomics/Sfitti_EVM.all.ID.gff3 -b ${1} --merge-strat single --no-extract-ase --no-quantify-graph --parallel 2 -v
