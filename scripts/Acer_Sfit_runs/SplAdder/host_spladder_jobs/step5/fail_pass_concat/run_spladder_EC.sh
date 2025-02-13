#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer'

#call events
spladder build -o ${workdir}/results/SplAdder/host_spladder_res/array_spladder_out2 -a ${workdir}/genomics/Acropora_cervicornis.ID.gff3 -b ${workdir}/scripts/SplAdder/host_spladder_jobs/host_bams_all.txt --parallel 40 -v
