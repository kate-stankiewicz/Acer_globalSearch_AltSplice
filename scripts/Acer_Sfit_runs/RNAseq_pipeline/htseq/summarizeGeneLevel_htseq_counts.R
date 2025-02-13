#set the command line args
args<-commandArgs(TRUE)
samp<-args[1]

#load libraries
library(dplyr)
library(stringr)
library(tidyr)

#read in the htseq file
test_htseq <- read.table(paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/htseq_results/host/no_host_",samp,"_htseq_count.txt"))

#remove htseq summary info at bottom on table (but save it separately to add back in later)
htseq_no_summary <-test_htseq[grep('^FUN', test_htseq$V1),]
summary_info <- test_htseq[grep('^FUN', test_htseq$V1, invert = T),]

#split isoform ID
test_htseq_sep <- separate(htseq_no_summary,V1,into = c("gene_name","isoform"),sep = "-",remove = FALSE,extra = "merge")

# sum counts of each isoform to get gene-level locus expression
test_htseq_gene_level <- aggregate(. ~ gene_name, data=test_htseq_sep[-c(1, 3)], FUN=sum)

#rename the column names back to V1 and V2 to match the summary information and concatenate them
colnames(test_htseq_gene_level) <- c("V1", "V2")

gene_level_htseq_output_format <- rbind(test_htseq_gene_level, summary_info)

#write this out with no header information
write.table(gene_level_htseq_output_format, paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/htseq_results/host/fix_", samp, "_htseq_count.txt"), row.names = F, col.names = F, quote = F)

