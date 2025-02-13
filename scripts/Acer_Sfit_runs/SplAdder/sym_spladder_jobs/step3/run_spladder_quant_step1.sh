#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#quantify for each sample separately
spladder build -o ${workdir}/results/SplAdder/sym_spladder_res/array_spladder_out -a ${workdir}/genomics/Sfitti_EVM.all.ID.gff3 -b ${1} --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode single --parallel 2 -v

