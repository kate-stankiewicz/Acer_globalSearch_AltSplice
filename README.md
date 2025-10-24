Repository for data processing and analysis scripts for Global Search project investigating alternative splicing in the coral Acropora cervicornis.

Input RNA-seq files can be accessed at NCBI BioProject PRJNA1213837.


Steps for running data analysis:

Part 1: Run RNA-seq pipeline (Acer_globalSearch_AltSplice/scripts/Acer_Sfit_runs/RNAseq_pipeline)

1. After downloading RNA-seq files (NCBI BioProject PRJNA1213837), run the GS pipeline (https://github.com/baliga-lab/Global_Search ) using the .json files and the metadata (Acer_globalSearch_AltSplice/metadata/AC_clean_meta_full.tsv)
2. Prepare the resulting BAM files for SplAdder (see ReadMe in prepare_bams_for_spladder)
3. Next run htseq

Part 2: Run the SplAdder pipeline (separately for the host and the symbiont; Acer_globalSearch_AltSplice/scripts/Acer_Sfit_runs/SplAdder)
1. Step1 generates the splice graph for each sample individually 
2. Step2 merges the splice graphs
3. Step3 quantifies for each sample
4. Step4 aggregates all quantifications
5. Step5 calls events
6. Run the statistical testing
    1. First run the test contrasts using the run_spladder_test_* files
    2. Next parse the output into one data frame using run_parse_spladder_output.R

Part 3: run data analysis in R (Acer_globalSearch_AltSplice/scripts/Acer_Sfit_runs/analysis)
1. Run differential expression analysis: use htseq output in DGE_Acer_GS_pilot_DESeq2.R
2. Run Acer_AS_analysis.R for main data exploration and figure generation (see comments for instructions within R script).
3. Run CBASS_schematic to generate an example schematic of experimental design
4.  Run AcerGSPilot_DGE_GOEnrichment.R for GO enrichment analysis
5. Run NMD_DGE_exploration.R to examine the splicing and gene expression of NMD related genes.
6. Run QY_MotePilot.R to run physiology analysis.




