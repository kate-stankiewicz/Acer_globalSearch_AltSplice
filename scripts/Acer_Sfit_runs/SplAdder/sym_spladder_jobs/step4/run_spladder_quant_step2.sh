#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#aggregate all quantification into one joint db
spladder build -o ${workdir}/results/SplAdder/sym_spladder_res/array_spladder_out -a ${workdir}/genomics/Sfitti_EVM.all.ID.gff3 -b ${workdir}/scripts/SplAdder/sym_spladder_jobs/sym_bams.txt  --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode collect --parallel 20 -v

