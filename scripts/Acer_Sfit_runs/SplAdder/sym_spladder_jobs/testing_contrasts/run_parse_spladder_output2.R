args<-commandArgs(TRUE)
var1<-args[1]
var2 <-args[2]
const <-args[3]

# read in the event files for each contrast
setwd(paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/SplAdder/array_jobs/host_spladder_jobs/array_spladder_out/testing_",var1,const,"_vs_",var2,const,"_Rem_Prob"))

ldf <- list() 
listcsv <- dir(pattern = "*gene_unique.tsv") 
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.table(file = listcsv[k], header = T)
}

# combine into one DF
one_ldf <- do.call(rbind, ldf)

#add a column with the constant
one_ldf$const <- paste0(const)

#add a column with the contrast variables
one_ldf$vrs <- paste0(var1,"vs",var2)

# add a column for species name
one_ldf$spec <- c("host")

#write out the results as one table
write.table(one_ldf, paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/all_samples/scripts/SplAdder/array_jobs/host_sym_parsed_event_files/host_all_events",var1,const,"_vs_",var2,const,".tsv"), quote = F, row.names = F)
