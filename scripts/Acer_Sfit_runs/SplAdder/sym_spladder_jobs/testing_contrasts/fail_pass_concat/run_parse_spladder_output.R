args<-commandArgs(TRUE)
var1<-args[1]
var2 <-args[2]

# read in the event files for each contrast
setwd(paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/SplAdder/sym_spladder_res/array_spladder_out2/testing_",var1,"_vs_",var2,"_all"))

ldf <- list() 
listcsv <- grep(list.files(pattern = "*.tsv"), pattern='gene.unique|extended', invert = TRUE, value = TRUE)
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.table(file = listcsv[k], header = T)
}

# combine into one DF
one_ldf <- do.call(rbind, ldf)

#add a column with the contrast variables
one_ldf$vrs <- paste0(var1,"_vs_",var2)

# add a column for species name
one_ldf$spec <- c("sym")

#write out the results as one table
write.table(one_ldf, paste0("/proj/omics4tb2/kstankiewicz/scratch/alt_splice/GS_Pilot_Acer/results/SplAdder/sym_spladder_res/array_spladder_out2/all_events_parsed/sym_all_events_",var1,"_vs_",var2,".tsv"), quote = F, row.names = F)
