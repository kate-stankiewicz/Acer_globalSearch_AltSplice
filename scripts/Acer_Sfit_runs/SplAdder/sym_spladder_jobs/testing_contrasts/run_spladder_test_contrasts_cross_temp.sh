#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
workdir='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/SplAdder/array_jobs/host_spladder_jobs'

#statistical test of all contrasts (run separately in job array)
spladder test -o ${workdir}/array_spladder_out --out-tag cross_temp --conditionA ${workdir}/testing_contrasts/contrast_files/host_AF_30.txt --conditionB ${workdir}/testing_contrasts/contrast_files/host_ICN_36.txt --labelA AF30 --labelB ICN36 --diagnose-plots -v --parallel 5
