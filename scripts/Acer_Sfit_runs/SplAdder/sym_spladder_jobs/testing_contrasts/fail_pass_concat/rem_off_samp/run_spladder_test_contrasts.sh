#!/bin/bash
# steps for running SplAdder on large cohorts

#set variables
out='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/SplAdder/sym_spladder_res/array_spladder_out2'
tests='/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/scripts/Acer_Sfit_runs/SplAdder/sym_spladder_jobs/testing_contrasts/fail_pass_concat/rem_off_samp/contrast_files'

#statistical test of all contrasts (run separately in job array)
spladder test -o ${out} --out-tag Rem_Prob --conditionA ${tests}/clean_${1}_test.txt --conditionB ${tests}/clean_${2}_test.txt --labelA ${1} --labelB ${2} --diagnose-plots -v --parallel 5
