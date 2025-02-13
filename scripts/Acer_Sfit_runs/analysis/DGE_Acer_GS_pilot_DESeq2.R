# GS Pilot Acer DESeq2 DGE: overlap with AS"

# Load libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(rstatix)
library(ggVennDiagram)
library(UpSetR)
library(apeglm)
library(data.table)

# read in the metadata
metadata <- read.table("~/alt_splice/GS_Pilot_Acer/metadata/AC_clean_meta_full.tsv", header = T)

# rename the column 
metadata <- metadata %>% dplyr::rename(sampleName = Novogene_R.)

# select only columns of interest
metadata_clean <- metadata %>% dplyr::select(Genotype, Timepoint, Treatment, treatment_timepoint, sampleName)
metadata_clean$Genotype <- factor(metadata_clean$Genotype)
metadata_clean$Timepoint <- factor(metadata_clean$Timepoint)
metadata_clean$Treatment <- factor(metadata_clean$Treatment)
metadata_clean$treatment_timepoint <- factor(metadata_clean$treatment_timepoint)



# read in the counts files and convert to DESeqDataSet object

# for the host
directory_host <- "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/RNASeq_pipeline/results_per_sample/htseq_results/host"
sampleFiles_host <- data.frame(grep("test_fix_", list.files(directory_host), value = T))
colnames(sampleFiles_host) <- "sampleName"
sampleFiles_host$fileName <- sampleFiles_host$sampleName
sampleFiles_host$sampleName <- gsub("test_fix_|_htseq_count.txt", "", sampleFiles_host$sampleName)
sampleTable_host <- merge(sampleFiles_host, metadata_clean, by = "sampleName")

dds_host <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_host,
                                       directory = directory_host,
                                       design= ~ treatment_timepoint)


# pre-filter the data

# host
keep_host <- rowSums(counts(dds_host)) >= 10
dds_host <- dds_host[keep_host,]


# host
# relevel to set the ref to initial timpoint
#dds_host$treatment_timepoint <-relevel(dds_host$treatment_timepoint, ref = "Initial_T0")
#dds_host$treatment_timepoint <- droplevels(dds_host$treatment_timepoint)
dds_host$treatment_timepoint <- factor(dds_host$treatment_timepoint, levels = c("Field_TF","Initial_T0", "Control_T1", "Control_T2", "Control_T3", "Control_T4", "Control_T5", "Low_T1", "Low_T2", "Low_T3", "Low_T4", "Low_T5", "Medium_T1", "Medium_T2", "Medium_T3", "Medium_T4", "Medium_T5", "High_T1", "High_T2", "High_T3", "High_T4", "High_T5"))

# run Deseq
dds_host <- DESeq(dds_host)

# quick check of names
resultsNames(dds_host)


#set up a table of possible testing pairs to explore
comp_against_all <- metadata_clean %>% get_comparisons("treatment_timepoint")
df_comp_all <- data.frame(t(as.data.frame(comp_against_all)))
colnames(df_comp_all) <- c("ref", "cond")

# split into timepoint and condition
comp_against_all_wide <- df_comp_all %>% tidyr::separate(ref, sep = "_", into = c("var1", "timepoint1"), remove = FALSE) %>% tidyr::separate(cond, sep = "_", into = c("var2", "timepoint2"), remove = FALSE)

# remove rows where timepoints are not the same (excluding initial and field)
comp_all_filt <- comp_against_all_wide %>% dplyr::filter(var1 == "Initial" | var1 == "Field" | var2 == "Initial" | var2 == "Field" | var1 == "Control" & timepoint1 == timepoint2 | var2 == "Control" & timepoint1 == timepoint2 )

# swap so that the initial and field are kept as the ref for testing purposes
comp_all_filt_add <- comp_all_filt %>% dplyr::mutate(SwapRef = ifelse(var1 == "Field" | 
                                             var1 == "Initial" & var2 != "Field" |
                                             var1 == "Control" & !(var2 %in% c("Field", "Initial")), ref, cond)) %>%
  dplyr::mutate(SwapCond= ifelse(var1 == "Field" | 
                                             var1 == "Initial" & var2 != "Field" |
                                             var1 == "Control" & !(var2 %in% c("Field", "Initial")), cond, ref))



# loop through the df to test each contrast type and write results to a file
for (i in rownames(comp_all_filt_add)){
  ref<-comp_all_filt_add[i,"SwapRef"]
  cond<-comp_all_filt_add[i,"SwapCond"]
  results_comp<-results(dds_host, contrast = c("treatment_timepoint", cond, ref))
  results_comp$GeneID=row.names(results_comp)
  results_comp<-results_comp[,c(7,1:6)]
  write.csv(data.frame(results_comp), 
              file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_#runs/analysis/DGE_DESeq2/DeSeq2Res_", cond, "_vs_", ref, ".csv"), row.names = F, quote = F)
}



# read back in the above files
setwd("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2")

DEG_list <- list() 
listcsv <- dir(pattern = "DeSeq2Res_") 
for (k in 1:length(listcsv)){
 DEG_list[[k]] <- read.csv(file = listcsv[k])
}

names(DEG_list) <- gsub("DeSeq2Res_|.csv","",
                       dir(pattern = "DeSeq2Res_"))

# filter by LFC and padj
#sig_DEG_list <- lapply(DEG_list, function(x) filter(x, padj < 0.01 & abs(log2FoldChange) > 1))

# write out the significant results
#mapply(function(dname, sig_DEG_list) 
#   write.csv(sig_DEG_list, file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sig_", dname, ".csv"), row.names = F, quote = F), 
#   names(sig_DEG_list), sig_DEG_list)



# Now we need to look for overlap in the genes with AS and the DEGs for each contrast
# But first we need to flip the names of the contrast for the AS to match what was set for the DGE (e.g., for DGE its High_T1_vs_Control_T1 but for AS its Control_T1_vs_High_T1)

# read in the significant differentially spliced genes
sig_AS_Acer <- read.table("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/Acer_all_sig_events_combined.tsv", header = T)

# get names from DGE
#DGE_names_list <- data.frame(names(sig_DEG_list))
#colnames(DGE_names_list) <- "vrs_format"

sig_AS_Acer<- sig_AS_Acer %>% tidyr::separate(vrs, sep = "_vs_", into = c("var1", "var2"), remove = FALSE)

# create a flipped var column for AS
sig_AS_Acer$FlipVrs <- paste0(sig_AS_Acer$var2, "_vs_", sig_AS_Acer$var1)

# convert these vrs columns to one long format column
sig_AS_Acer <- sig_AS_Acer %>% tidyr::pivot_longer(cols = c(vrs,FlipVrs), names_to = "vrs_type", values_to = "vrs_format")

# sanity check: see which contrast format matches
#sig_AS_Acer$match <-sig_AS_Acer$vrs_format %in% DGE_names_list$vrs_format

# join the two dataframes to filter for matching contrast name
AS_DGE_cont_name_join <- merge(DGE_names_list, sig_AS_Acer, by = "vrs_format")



# add the contrast name to a column of each df for the DGE
sig_DEG_list_named <- mapply(cbind, sig_DEG_list, "vrs"=names(sig_DEG_list), SIMPLIFY=F)

# convert to one huge DF
DEG_one_df <- do.call(rbind, sig_DEG_list_named)

# split the transcript ID and base gene name
DEG_one_df<- DEG_one_df %>% tidyr::separate(GeneID, sep = "-", into = c("gene_id", "transcript_id"), remove = FALSE)

# get a list of genes, regardless of contrast for AS and DE
all_AS_DE_genes <- list(AS = unique(AS_DGE_cont_name_join$gene_id), DE = unique(DEG_one_df$gene_id))

# create venn diagram to check overlap
ggVenn_all <- ggVennDiagram(all_AS_DE_genes)


# next we need to check for overlap within each specific contrast

# first add an identifier to the contrast name to indicate AS or DE
AS_DGE_cont_name_join$contrast_type <- paste0("AS_", AS_DGE_cont_name_join$vrs_format)
DEG_one_df$contrast_type <- paste0("DE_", DEG_one_df$vrs)

# select only the columns of interest in each
DEG_one_df_genes_cont <- DEG_one_df %>% dplyr::select(gene_id, contrast_type)
AS_genes_cont <- AS_DGE_cont_name_join %>% dplyr::select(gene_id, contrast_type)

# bind these two together
AS_DE_genes_one_df <- rbind(AS_genes_cont, DEG_one_df_genes_cont)

# create a binary matrix of this
binary_AS_DE <- dcast(AS_DE_genes_one_df, gene_id~contrast_type, length)

# move the gene_id column to the rownames
binary_AS_DE <- binary_AS_DE %>% remove_rownames %>% column_to_rownames(var="gene_id")

# create upset plot for timepoint specific genes
upset_AS_DE_all <- UpSetR::upset(binary_AS_DE, nsets = 112, order.by = "freq", nintersects = 40)



# create a basic column to look for overlap of contrasts regardless of timepoint
AS_DE_genes_one_df$basic_contrast <- gsub("_TF|_T0|_T1|_T2|_T3|_T4|_T5", "", AS_DE_genes_one_df$contrast_type)

# create a binary matrix of this
counts_AS_DE_basic <- dcast(AS_DE_genes_one_df, gene_id~basic_contrast, length)
counts_AS_DE_basic <- counts_AS_DE_basic %>% remove_rownames %>% column_to_rownames(var="gene_id")
binary_AS_DE_basic <- data.frame(as.matrix((counts_AS_DE_basic > 0) + 0))

# create upset plot for timepoint specific genes
upset_AS_DE_all_basic <- UpSetR::upset(binary_AS_DE_basic, nsets = 24, order.by = "freq", nintersects = 100)



# convert to list of lists
AS_DE_all_lists <- AS_DE_genes_one_df %>% dplyr::select(gene_id,contrast_type) %>% dplyr::group_by(contrast_type) %>% dplyr::summarise(named_vec = list(gene_id)) %>% deframe()

AS_DE_basic_lists <- AS_DE_genes_one_df %>% dplyr::select(gene_id,basic_contrast) %>% dplyr::group_by(basic_contrast) %>% dplyr::summarise(named_vec = list(gene_id)) %>% deframe()

# make each list unique
AS_DE_all_lists_uni <- lapply(AS_DE_all_lists, unique)
AS_DE_basic_lists_uni <- lapply(AS_DE_basic_lists, unique)

# calculate the overlap for each pairing
overlap_all_pairwise_mat <- data.frame(crossprod(table(stack(AS_DE_all_lists_uni))))
overlap_basic_pairwise_mat <- data.frame(crossprod(table(stack(AS_DE_basic_lists_uni))))

# subset the matrices so that only DE and AS are being compared (not AS to AS or DE to DE)
overlap_basic_subset <- overlap_basic_pairwise_mat %>% dplyr:: select(starts_with("AS_"))
overlap_basic_subset <- overlap_basic_subset %>% filter(grepl("DE_", rownames(overlap_basic_subset)))

overlap_all_subset <- overlap_all_pairwise_mat %>% dplyr:: select(starts_with("AS_"))
overlap_all_subset <- overlap_all_subset %>% filter(grepl("DE_", rownames(overlap_all_subset)))

# turn the overlap matricesinto long format
overlap_basic_subset_add <- overlap_basic_subset %>% rownames_to_column(var="var1")
long_basic_overlap <- overlap_basic_subset_add %>% pivot_longer(!var1, names_to = "var2", values_to = "overlap")

overlap_all_subset_add <- overlap_all_subset %>% rownames_to_column(var="var1")
long_all_overlap <- overlap_all_subset_add %>% pivot_longer(!var1, names_to = "var2", values_to = "overlap")

# bind the above together
all_long <- rbind(long_basic_overlap, long_all_overlap)

# paste together the comparison terms for the overlap (which DE and which AS contrast)
all_long$comparison <- paste0(all_long$var1, "_", all_long$var2)

# get contrast base name for each
all_long$cont1_base <- sub('.*?[_]', "", all_long$var1)
all_long$cont2_base <- sub('.*?[_]', "", all_long$var2)

# get the maximum overlap with DE for each AS contrast type
all_long_max <- all_long %>% dplyr::group_by(var2) %>% dplyr::mutate(largest_cont = var1[which.max(overlap)]) %>% dplyr::ungroup()
all_long_max <- all_long_max %>% dplyr::group_by(var2) %>% dplyr::mutate(largest_val = max(overlap)) %>% dplyr::ungroup()

# add a column with check for if a comparison is equal to the max overlap value (helps identify cases of "ties" for max within a group)
all_long_max$is_max <- all_long_max$largest_val == all_long_max$overlap

# count how many matches per group (should only be one unless there are ties) 
all_long_max <- all_long_max %>% dplyr::group_by(var2) %>% dplyr::mutate(num_maxes = sum(is_max == "TRUE")) # looks like there are 2 cases of ties here

# if a comparison has a max overlap, check whether its the same contrast for DE and AS
all_long_max <- all_long_max %>% dplyr::mutate(max_is_match = ifelse(is_max == "TRUE" & cont1_base == cont2_base, "yes", "no"))

# ID all the rows where the contrasts match between DE and AS
all_long_max <- all_long_max %>% dplyr::mutate(cont_match = ifelse(cont1_base == cont2_base, "yes", "no"))



# count the total number of genes for each contrast type
counts_AS_DE_all <- AS_DE_genes_one_df %>% dplyr::group_by(contrast_type) %>% dplyr::summarise(count = n_distinct(gene_id))
counts_AS_DE_basic <- AS_DE_genes_one_df %>% dplyr::group_by(basic_contrast) %>% dplyr::summarise(count = n_distinct(gene_id))

# put the two counts tables together
counts_AS_DE_basic <- counts_AS_DE_basic %>% dplyr::rename(contrast_type = basic_contrast)
counts_all <- rbind(counts_AS_DE_all,counts_AS_DE_basic) 

# split the AS/DE designator from the contrast type
counts_all <- separate(counts_all,contrast_type,into = c("type","cont"),sep = "_",remove = FALSE,extra = "merge")

# get which type has the most genes per contrast (DE or AS)
counts_all <- counts_all %>% dplyr::group_by(cont) %>% dplyr::mutate(largest = type[which.max(count)]) %>% dplyr::ungroup()

# since DE is always bigger, we will filter for just counts of AS to get proportion of each overlap relative to total AS
AS_counts_all <- counts_all %>% dplyr::filter(type == "AS")



# get proportion of overlap with DE relative to total AS genes per contrast 
# rename the key columns in the two DFs to match and select the needed columns for the AS counts df
AS_counts_all_rename <- AS_counts_all %>% dplyr::rename(AS_cont = contrast_type) %>% dplyr::select(AS_cont, count)
all_long_max_rename <- all_long_max %>% dplyr::rename(AS_cont = var2, DE_cont = var1)

# join the counts to the overlap summary
AS_counts_overlap_join <- dplyr::left_join(all_long_max_rename, AS_counts_all_rename, by = "AS_cont")

# calculate proportion of overlap relative to all AS genes
AS_counts_overlap_join <- AS_counts_overlap_join %>% dplyr::mutate(overlap_perc_AS = overlap / count * 100)

# remove or rename some columns and write out this file to csv
AS_counts_overlap_join <- AS_counts_overlap_join %>% dplyr::select(DE_cont, AS_cont, overlap, largest_cont, largest_val, is_max, cont_match, count, overlap_perc_AS) %>% dplyr::rename(largest_overlap_contrast = largest_cont, largest_overlap_size = largest_val, contrasts_match = cont_match, count_AS = count)

#write.csv(AS_counts_overlap_join, file ="~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/AS_DE_contrasts_overlap.csv", row.names = F, quote = F)


# filter for just overlaps where DE and AS contrast matches
overlap_contrasts_match <- AS_counts_overlap_join %>% dplyr::filter(contrasts_match == "yes")

# add a column for the base contrast name
overlap_contrasts_match$comparison <- overlap_contrasts_match$AS_cont 
overlap_contrasts_match$comparison <- gsub("AS_", "", overlap_contrasts_match$comparison)

# mark each contrast as basic or complex (containing timepoint or not)
overlap_contrasts_match <- overlap_contrasts_match %>% dplyr::mutate(contrast_type = ifelse(str_detect(comparison, "T"), "complex", "basic"))

# split by basic and complex contrast type
overlap_basic_contrasts <- overlap_contrasts_match %>% dplyr::filter(contrast_type == "basic")
overlap_complex_contrasts <- overlap_contrasts_match %>% dplyr::filter(contrast_type == "complex")

# add in the count of differentially expressed genes
counts_DE <- counts_all %>% dplyr::filter(type == "DE") %>% select(cont, count)
counts_DE <- counts_DE %>% dplyr::rename(comparison = cont, count_DE = count)
overlap_basic_contrasts_all_counts <- dplyr::left_join(overlap_basic_contrasts, counts_DE, by = "comparison")
overlap_complex_contrasts_all_counts <- dplyr::left_join(overlap_complex_contrasts, counts_DE, by = "comparison")

# make the counts one column witha  DE or AS lable
overlap_counts_basic_long <- overlap_basic_contrasts_all_counts %>% tidyr::pivot_longer(cols = c(count_AS,count_DE), names_to = "count_type", values_to = "count")
overlap_counts_basic_long$count_type <- gsub("count_", "", overlap_counts_basic_long$count_type)

overlap_counts_complex_long <- overlap_complex_contrasts_all_counts %>% tidyr::pivot_longer(cols = c(count_AS,count_DE), names_to = "count_type", values_to = "count")
overlap_counts_complex_long$count_type <- gsub("count_", "", overlap_counts_complex_long$count_type)


# plot the counts with the overlap

# set the category order for the x axis
overlap_counts_basic_long$comparison <- factor(overlap_counts_basic_long$comparison, levels= c("Initial_vs_Field", "Control_vs_Field", "Control_vs_Initial", "Low_vs_Field", "Low_vs_Initial", "Low_vs_Control", "Medium_vs_Field", "Medium_vs_Initial", "Medium_vs_Control", "High_vs_Field", "High_vs_Initial", "High_vs_Control"))

overlap_plot <- ggplot(overlap_counts_basic_long) + geom_bar(aes(x = comparison, y = count, fill = count_type), stat = "identity", position = "dodge") + geom_point(aes(x = comparison, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + geom_line(aes(x = comparison, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + scale_y_continuous(sec.axis = sec_axis(~./150, name="percent overlap")) + theme_bw() + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))

# do the same for the complex contrast types

# select just the control vs the treatment comparisons
control_only_complex <- overlap_counts_complex_long %>% dplyr::filter(!str_detect(comparison, "Field|Initial"))

# order the columns
control_only_complex$treatment <- sub("\\_.*", "", control_only_complex$comparison)
control_only_complex$timepoint <-str_sub(control_only_complex$comparison, -2, -1)
control_only_complex$treatment <- factor(control_only_complex$treatment, levels = c("Low", "Medium", "High"))


overlap_plot_complex <- ggplot(control_only_complex) + geom_bar(aes(x = timepoint, y = count, fill = count_type), stat = "identity", position = "dodge") + geom_point(aes(x = timepoint, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + geom_line(aes(x = timepoint, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + scale_y_continuous(sec.axis = sec_axis(~./150, name="percent overlap")) + theme_bw() + facet_grid(.~treatment)

#re-do the above with lfc shrinkage method

# loop through the df to test each contrast type and write results to a file
for (i in rownames(comp_all_filt_add)){
  ref<-comp_all_filt_add[i,"SwapRef"]
  cond<-comp_all_filt_add[i,"SwapCond"]
  dds_host$treatment_timepoint <- relevel(dds_host$treatment_timepoint, ref = ref)
  dds_host <- DESeq(dds_host)
  results_comp_lfc<-lfcShrink(dds_host, coef = paste0("treatment_timepoint_", cond, "_vs_", ref), type = "apeglm", svalue = TRUE)
  df_res_lfc <- data.frame(GeneID = rownames(results_comp_lfc), results_comp_lfc)
  write.csv(df_res_lfc, file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sumGeneDESeq2", cond, "_vs_", ref, "_sval_lfc.csv"), row.names = F, quote = F)
}


# read back in the above files
setwd("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2")

DEG_list_lfc <- list() 
listcsv_lfc <- dir(pattern = "sumGeneDESeq2") 
for (k in 1:length(listcsv_lfc)){
 DEG_list_lfc[[k]] <- read.csv(file = listcsv_lfc[k])
}

names(DEG_list_lfc) <- gsub("sumGeneDESeq2|_sval_lfc.csv","",
                       dir(pattern = "sumGeneDESeq2"))

# filter for significance
#sig_DEG_list_lfc <- lapply(DEG_list_lfc, function(x) filter(x, padj < 0.05))
sig_DEG_list_lfc <- lapply(DEG_list_lfc, function(x) dplyr::filter(x, svalue < 0.005 & abs(log2FoldChange) > 1))

# write out the significant results
#mapply(function(dname, sig_DEG_list_lfc) 
#   write.csv(sig_DEG_list_lfc, file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sig_sumGene_", dname, "sval_lfc.csv"), row.names = F, quote = F), 
#   names(sig_DEG_list_lfc), sig_DEG_list_lfc)



# get names from DGE
lfc_DGE_names_list <- data.frame(names(sig_DEG_list_lfc))
colnames(lfc_DGE_names_list) <- "vrs_format"

# Now we need to look for overlap in the genes with AS and the DEGs for each contrast
# But first we need to flip the names of the contrast for the AS to match what was set for the DGE (e.g., for DGE its High_T1_vs_Control_T1 but for AS its Control_T1_vs_High_T1)

# read in the significant differentially spliced genes
sig_AS_Acer <- read.table("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/Acer_all_sig_events_combined.tsv", header = T)

# split the vars for the AS
sig_AS_Acer<- sig_AS_Acer %>% tidyr::separate(vrs, sep = "_vs_", into = c("var1", "var2"), remove = FALSE)

# create a flipped var column for AS
sig_AS_Acer$FlipVrs <- paste0(sig_AS_Acer$var2, "_vs_", sig_AS_Acer$var1)

# convert these vrs columns to one long format column
sig_AS_Acer <- sig_AS_Acer %>% tidyr::pivot_longer(cols = c(vrs,FlipVrs), names_to = "vrs_type", values_to = "vrs_format")

# sanity check: see which contrast format matches
sig_AS_Acer$match <-sig_AS_Acer$vrs_format %in% lfc_DGE_names_list$vrs_format

# join the two dataframes to filter for matching contrast name
AS_DGE_cont_name_join <- merge(lfc_DGE_names_list, sig_AS_Acer, by = "vrs_format")

# join the two dataframes to filter for matching contrast name
lfc_AS_DGE_cont_name_join <- merge(lfc_DGE_names_list, sig_AS_Acer, by = "vrs_format")



# add the contrast name to a column of each df for the DGE
lfc_sig_DEG_list_named <- mapply(cbind, sig_DEG_list_lfc, "vrs"=names(sig_DEG_list_lfc), SIMPLIFY=F)

# convert to one huge DF
lfc_DEG_one_df <- do.call(rbind, lfc_sig_DEG_list_named)

#rename the Gene ID column to match
lfc_DEG_one_df <- lfc_DEG_one_df %>% dplyr::rename(gene_id = GeneID)

#write out this table for future use
#write.csv(lfc_DEG_one_df, file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sumGene_all_sig_DGE_lfc_sval.csv"), row.names = F, quote = F)

# create unique ID for each gene for each contrast
lfc_DEG_one_df$gene_id_cont <- paste0(lfc_DEG_one_df$gene_id, "_", lfc_DEG_one_df$vrs)

# group by this new ID and take the most significant value (lowest svalue)
lfc_DEG_one_df_col_iso <- lfc_DEG_one_df %>% dplyr::group_by(gene_id_cont) %>% dplyr::slice_min(svalue) %>% dplyr::ungroup()

# get a list of genes, regardless of contrast for AS and DE
lfc_all_AS_DE_genes <- list(AS = unique(lfc_AS_DGE_cont_name_join$gene_id), DE = unique(lfc_DEG_one_df_col_iso$gene_id))

# create venn diagram to check overlap
lfc_ggVenn_all <- ggVennDiagram(lfc_all_AS_DE_genes)



# next we need to check for overlap within each specific contrast

# first add an identifier to the contrast name to indicate AS or DE
lfc_AS_DGE_cont_name_join$contrast_type <- paste0("AS_", lfc_AS_DGE_cont_name_join$vrs_format)
lfc_DEG_one_df_col_iso$contrast_type <- paste0("DE_", lfc_DEG_one_df_col_iso$vrs)

# select only the columns of interest in each
lfc_DEG_one_df_genes_cont <- lfc_DEG_one_df_col_iso %>% dplyr::select(gene_id, contrast_type)
lfc_AS_genes_cont <- lfc_AS_DGE_cont_name_join %>% dplyr::select(gene_id, contrast_type)

# bind these two together
lfc_AS_DE_genes_one_df <- rbind(lfc_AS_genes_cont, lfc_DEG_one_df_genes_cont)


# create a basic column to look for overlap of contrasts regardless of timepoint
lfc_AS_DE_genes_one_df$basic_contrast <- gsub("_TF|_T0|_T1|_T2|_T3|_T4|_T5", "", lfc_AS_DE_genes_one_df$contrast_type)

# convert to list of lists
lfc_AS_DE_all_lists <- lfc_AS_DE_genes_one_df %>% dplyr::select(gene_id,contrast_type) %>% dplyr::group_by(contrast_type) %>% dplyr::summarise(named_vec = list(gene_id)) %>% deframe()

lfc_AS_DE_basic_lists <- lfc_AS_DE_genes_one_df %>% dplyr::select(gene_id,basic_contrast) %>% dplyr::group_by(basic_contrast) %>% dplyr::summarise(named_vec = list(gene_id)) %>% deframe()

# make each list unique
lfc_AS_DE_all_lists_uni <- lapply(lfc_AS_DE_all_lists, unique)
lfc_AS_DE_basic_lists_uni <- lapply(lfc_AS_DE_basic_lists, unique)

# calculate the overlap for each pairing
lfc_overlap_all_pairwise_mat <- data.frame(crossprod(table(stack(lfc_AS_DE_all_lists_uni))))
lfc_overlap_basic_pairwise_mat <- data.frame(crossprod(table(stack(lfc_AS_DE_basic_lists_uni))))

# subset the matrices so that only DE and AS are being compared (not AS to AS or DE to DE)
lfc_overlap_basic_subset <- lfc_overlap_basic_pairwise_mat %>% dplyr:: select(starts_with("AS_"))
lfc_overlap_basic_subset <- lfc_overlap_basic_subset %>% filter(grepl("DE_", rownames(lfc_overlap_basic_subset)))

lfc_overlap_all_subset <- lfc_overlap_all_pairwise_mat %>% dplyr:: select(starts_with("AS_"))
lfc_overlap_all_subset <- lfc_overlap_all_subset %>% filter(grepl("DE_", rownames(lfc_overlap_all_subset)))

# turn the overlap matricesinto long format
lfc_overlap_basic_subset_add <- lfc_overlap_basic_subset %>% rownames_to_column(var="var1")
lfc_long_basic_overlap <- lfc_overlap_basic_subset_add %>% pivot_longer(!var1, names_to = "var2", values_to = "overlap")

lfc_overlap_all_subset_add <- lfc_overlap_all_subset %>% rownames_to_column(var="var1")
lfc_long_all_overlap <- lfc_overlap_all_subset_add %>% pivot_longer(!var1, names_to = "var2", values_to = "overlap")

# bind the above together
lfc_all_long <- rbind(lfc_long_basic_overlap, lfc_long_all_overlap)

# paste together the comparison terms for the overlap (which DE and which AS contrast)
lfc_all_long$comparison <- paste0(lfc_all_long$var1, "_", lfc_all_long$var2)

# get contrast base name for each
lfc_all_long$cont1_base <- sub('.*?[_]', "", lfc_all_long$var1)
lfc_all_long$cont2_base <- sub('.*?[_]', "", lfc_all_long$var2)

# get the maximum overlap with DE for each AS contrast type
lfc_all_long_max <- lfc_all_long %>% dplyr::group_by(var2) %>% dplyr::mutate(largest_cont = var1[which.max(overlap)]) %>% dplyr::ungroup()
lfc_all_long_max <- lfc_all_long_max %>% dplyr::group_by(var2) %>% dplyr::mutate(largest_val = max(overlap)) %>% dplyr::ungroup()

# add a column with check for if a comparison is equal to the max overlap value (helps identify cases of "ties" for max within a group)
lfc_all_long_max$is_max <- lfc_all_long_max$largest_val == lfc_all_long_max$overlap

# count how many matches per group (should only be one unless there are ties) 
lfc_all_long_max <- lfc_all_long_max %>% dplyr::group_by(var2) %>% dplyr::mutate(num_maxes = sum(is_max == "TRUE")) # looks like there are 2 cases of ties here

# if a comparison has a max overlap, check whether its the same contrast for DE and AS
lfc_all_long_max <- lfc_all_long_max %>% dplyr::mutate(max_is_match = ifelse(is_max == "TRUE" & cont1_base == cont2_base, "yes", "no"))

# ID all the rows where the contrasts match between DE and AS
lfc_all_long_max <- lfc_all_long_max %>% dplyr::mutate(cont_match = ifelse(cont1_base == cont2_base, "yes", "no"))



# count the total number of genes for each contrast type
lfc_counts_AS_DE_all <- lfc_AS_DE_genes_one_df %>% dplyr::group_by(contrast_type) %>% dplyr::summarise(count = n_distinct(gene_id))
lfc_counts_AS_DE_basic <- lfc_AS_DE_genes_one_df %>% dplyr::group_by(basic_contrast) %>% dplyr::summarise(count = n_distinct(gene_id))

# put the two counts tables together
lfc_counts_AS_DE_basic <- lfc_counts_AS_DE_basic %>% dplyr::rename(contrast_type = basic_contrast)
lfc_counts_all <- rbind(lfc_counts_AS_DE_all,lfc_counts_AS_DE_basic) 

# split the AS/DE designator from the contrast type
lfc_counts_all <- separate(lfc_counts_all, contrast_type,into = c("type","cont"),sep = "_",remove = FALSE,extra = "merge")

# get which type has the most genes per contrast (DE or AS)
lfc_counts_all <- lfc_counts_all %>% dplyr::group_by(cont) %>% dplyr::mutate(largest = type[which.max(count)]) %>% dplyr::ungroup()

# since DE is always bigger, we will filter for just counts of AS to get proportion of each overlap relative to total AS
lfc_AS_counts_all <- lfc_counts_all %>% dplyr::filter(type == "AS")



# get proportion of overlap with DE relative to total AS genes per contrast 
# rename the key columns in the two DFs to match and select the needed columns for the AS counts df
lfc_AS_counts_all_rename <- lfc_AS_counts_all %>% dplyr::rename(AS_cont = contrast_type) %>% dplyr::select(AS_cont, count)
lfc_all_long_max_rename <- lfc_all_long_max %>% dplyr::rename(AS_cont = var2, DE_cont = var1)

# join the counts to the overlap summary
lfc_AS_counts_overlap_join <- dplyr::left_join(lfc_all_long_max_rename, lfc_AS_counts_all_rename, by = "AS_cont")

# calculate proportion of overlap relative to all AS genes
lfc_AS_counts_overlap_join <- lfc_AS_counts_overlap_join %>% dplyr::mutate(overlap_perc_AS = overlap / count * 100)

# remove or rename some columns and write out this file to csv
lfc_AS_counts_overlap_join <- lfc_AS_counts_overlap_join %>% dplyr::select(DE_cont, AS_cont, overlap, largest_cont, largest_val, is_max, cont_match, count, overlap_perc_AS) %>% dplyr::rename(largest_overlap_contrast = largest_cont, largest_overlap_size = largest_val, contrasts_match = cont_match, count_AS = count)

#write.csv(lfc_AS_counts_overlap_join, file = paste0("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sumGene_AS_DE_contrasts_overlap_lfc_sval.csv"), row.names = F, quote = F)



# filter for just overlaps where DE and AS contrast matches
lfc_overlap_contrasts_match <- lfc_AS_counts_overlap_join %>% dplyr::filter(contrasts_match == "yes")

# add a column for the base contrast name
lfc_overlap_contrasts_match$comparison <- lfc_overlap_contrasts_match$AS_cont 
lfc_overlap_contrasts_match$comparison <- gsub("AS_", "", lfc_overlap_contrasts_match$comparison)

# mark each contrast as basic or complex (containing timepoint or not)
lfc_overlap_contrasts_match <- lfc_overlap_contrasts_match %>% dplyr::mutate(contrast_type = ifelse(str_detect(comparison, "T"), "complex", "basic"))

# split by basic and complex contrast type
lfc_overlap_basic_contrasts <- lfc_overlap_contrasts_match %>% dplyr::filter(contrast_type == "basic")
lfc_overlap_complex_contrasts <- lfc_overlap_contrasts_match %>% dplyr::filter(contrast_type == "complex")

# add in the count of differentially expressed genes
lfc_counts_DE <- lfc_counts_all %>% dplyr::filter(type == "DE") %>% select(cont, count)
lfc_counts_DE <- lfc_counts_DE %>% dplyr::rename(comparison = cont, count_DE = count)
lfc_overlap_basic_contrasts_all_counts <- dplyr::left_join(lfc_overlap_basic_contrasts, lfc_counts_DE, by = "comparison")
lfc_overlap_complex_contrasts_all_counts <- dplyr::left_join(lfc_overlap_complex_contrasts, lfc_counts_DE, by = "comparison")

# make the counts one column witha  DE or AS lable
lfc_overlap_counts_basic_long <- lfc_overlap_basic_contrasts_all_counts %>% tidyr::pivot_longer(cols = c(count_AS,count_DE), names_to = "count_type", values_to = "count")
lfc_overlap_counts_basic_long$count_type <- gsub("count_", "", lfc_overlap_counts_basic_long$count_type)

lfc_overlap_counts_complex_long <- lfc_overlap_complex_contrasts_all_counts %>% tidyr::pivot_longer(cols = c(count_AS,count_DE), names_to = "count_type", values_to = "count")
lfc_overlap_counts_complex_long$count_type <- gsub("count_", "", lfc_overlap_counts_complex_long$count_type)


# plot the counts with the overlap

# set the category order for the x axis
lfc_overlap_counts_basic_long$comparison <- factor(lfc_overlap_counts_basic_long$comparison, levels= c("Initial_vs_Field", "Control_vs_Field", "Control_vs_Initial", "Low_vs_Field", "Low_vs_Initial", "Low_vs_Control", "Medium_vs_Field", "Medium_vs_Initial", "Medium_vs_Control", "High_vs_Field", "High_vs_Initial", "High_vs_Control"))

lfc_overlap_plot <- ggplot(lfc_overlap_counts_basic_long) + geom_bar(aes(x = comparison, y = count, fill = count_type), stat = "identity", position = "dodge") + geom_point(aes(x = comparison, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + geom_line(aes(x = comparison, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + scale_y_continuous(sec.axis = sec_axis(~./150, name="percent overlap")) + theme_bw() + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))

# do the same for the complex contrast types

# select just the control vs the treatment comparisons
lfc_control_only_complex <- lfc_overlap_counts_complex_long %>% dplyr::filter(!str_detect(comparison, "Field|Initial"))

# order the columns
lfc_control_only_complex$treatment <- sub("\\_.*", "", lfc_control_only_complex$comparison)
lfc_control_only_complex$timepoint <-str_sub(lfc_control_only_complex$comparison, -2, -1)
lfc_control_only_complex$treatment <- factor(lfc_control_only_complex$treatment, levels = c("Low", "Medium", "High"))


lfc_overlap_plot_complex <- ggplot(lfc_control_only_complex) + geom_bar(aes(x = timepoint, y = count, fill = count_type), stat = "identity", position = "dodge") + geom_point(aes(x = timepoint, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + geom_line(aes(x = timepoint, y = overlap_perc_AS*150, group = 1), stat = "identity", color = "black") + scale_y_continuous(sec.axis = sec_axis(~./150, name="percent overlap")) + theme_bw() + facet_grid(.~treatment)


# Next let's take a closer look at the differentially expressed genes that are also DAS

# first filter select the columns of interest
lfc_DEG_col_int <- lfc_DEG_one_df_col_iso %>% dplyr::select(gene_id, log2FoldChange, svalue, contrast_type) %>% tidyr::separate(contrast_type,into = c("type","cont"),sep = "_",remove = FALSE,extra = "merge")

# split the AS event id column and select only cols of interest
AS_col_int <- lfc_AS_DGE_cont_name_join %>% tidyr::separate(event_id, sep = "\\.", into = c("event_type", "event_num"), remove = FALSE) %>% dplyr::select(gene_id, event_type, event_num, dPSI, contrast_type) %>% tidyr::separate(contrast_type,into = c("type","cont"),sep = "_",remove = FALSE,extra = "merge")

# join with AS
AS_DE_wide_join <- dplyr::full_join(AS_col_int, lfc_DEG_col_int, by=c("gene_id", "cont"))

# Next filter for only contrasts relating to controls
AS_DE_wide_join_controls <- AS_DE_wide_join %>% dplyr::filter(!str_detect(cont, "Field|Initial"))

#filter for NAs
AS_DE_overlap_wide <- na.omit(AS_DE_wide_join_controls)

# filter for dPSI > 0.3
AS_DE_overlap_0.3dPSI <- AS_DE_overlap_wide %>% dplyr::filter(abs(dPSI) > 0.3)



# we need to filter significant events with a lot of missing information in one (or both) of the testing groups

##read in the PSI values

setwd("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/SplAdder/host_spladder_res/array_spladder_out2")

host_ldf <- list() 
listcsv_h <- dir(pattern = "_C3.confirmed.txt") 
for (k in 1:length(listcsv_h)){
 host_ldf[[k]] <- read.table(file = listcsv_h[k], header = T)
}

# get only psi and gene name columns
host_psi <- lapply(host_ldf, function(x) x[, grepl('psi|gene|event' , names( x ) ) ])

# convert to one dataframe
host_one_df <- do.call(rbind, host_psi)


# make it long format
long_PSI <- host_one_df %>% tidyr::pivot_longer(cols = c(-event_id,  -gene_name), names_to = "sample_id", values_to = "psi")

# add a samplename column 
long_PSI$sampleName <- gsub("host_|.sorted.psi|A|B", "", long_PSI$sample_id)

#merge with metadata
long_PSI_meta <- merge(long_PSI, metadata, by = "sampleName")


# for each treatment + timepoint, count the number of NaNs
NACounts_PSI <- long_PSI_meta %>% dplyr::group_by(event_id, treatment_timepoint) %>% dplyr::summarise(count_NaN = sum(psi == "NaN", na.rm = F))

# filter for events with more NaN values for more than 4 samples
high_NaN <- NACounts_PSI %>% dplyr::filter(count_NaN > 4)

#split the sig results contrast into var1 and var2
AS_DE_overlap_0.3dPSI_vars <- AS_DE_overlap_0.3dPSI %>% tidyr::separate(cont, sep = "_vs_", into = c("var1", "var2"), remove = FALSE)

#concat the event type and number again for matching purposes
AS_DE_overlap_0.3dPSI_vars$event_id <- paste0(AS_DE_overlap_0.3dPSI_vars$event_type, ".", AS_DE_overlap_0.3dPSI_vars$event_num)

# create concatenated ID
AS_DE_overlap_0.3dPSI_vars$key1 <- paste0(AS_DE_overlap_0.3dPSI_vars$event_id, "_", AS_DE_overlap_0.3dPSI_vars$var1)

AS_DE_overlap_0.3dPSI_vars$key2 <- paste0(AS_DE_overlap_0.3dPSI_vars$event_id, "_", AS_DE_overlap_0.3dPSI_vars$var2)

high_NaN$key <- paste0(high_NaN$event_id, "_", high_NaN$treatment_timepoint)

# match for var1
AS_DE_overlap_0.3dPSI_vars_match <- AS_DE_overlap_0.3dPSI_vars %>% dplyr::mutate(match1 = ifelse(key1 %in% high_NaN$key | key2 %in% high_NaN$key, "yes", "no"))

# filter out the ones that match to high NaN proportions
AS_DE_full_filt <- AS_DE_overlap_0.3dPSI_vars_match %>% dplyr::filter(match1 == "no")

#is there more or less AS in the DE genes?
sum(AS_DE_full_filt$dPSI < 0)
sum(AS_DE_full_filt$dPSI > 0)

# make unique by gene
genes_unique_AS_DE_fullFilt <- AS_DE_full_filt %>% dplyr::distinct(gene_id, cont, .keep_all = TRUE)

sum(genes_unique_AS_DE_fullFilt$log2FoldChange < 0)
sum(genes_unique_AS_DE_fullFilt$log2FoldChange > 0)


#the sign of the dPSI is flipped relative to the logFC for gene expression, so lets fix that
AS_DE_full_filt$flipped_dPSI <- AS_DE_full_filt$dPSI * -1

# test the correlation between dPSI and logFC
cor(AS_DE_full_filt$flipped_dPSI, AS_DE_full_filt$log2FoldChange, method = "pearson")

#filter into up and down reg by gene expression
up_reg <- AS_DE_full_filt %>% dplyr::filter(log2FoldChange > 0)
down_reg <- AS_DE_full_filt %>% dplyr::filter(log2FoldChange < 0)

# is there more or less AS in up and down reg?
sum(up_reg$flipped_dPSI < 0)
sum(up_reg$flipped_dPSI > 0)

sum(down_reg$flipped_dPSI < 0)
sum(down_reg$flipped_dPSI > 0)

length(unique(up_reg$gene_id))
length(unique(down_reg$gene_id))

table(up_reg$cont)
table(down_reg$cont)

more_AS <- AS_DE_full_filt %>% dplyr::filter(flipped_dPSI > 0)
less_AS <- AS_DE_full_filt %>% dplyr::filter(flipped_dPSI < 0)

table(more_AS$cont)
table(less_AS$cont)


# count how many contrasts each event is present in
counts_by_event <- AS_DE_full_filt %>% dplyr::add_count(event_id, name = "num_event_id")

# filter events that occur <4x (we can't track them through time easily)
counts_by_event_filt <- counts_by_event %>% dplyr::filter(num_event_id >= 4)


# split var1 into treatment and timepoint 
counts_by_event_filt <- counts_by_event_filt %>% tidyr::separate(var1, sep = "_", into = c("treatment", "timepoint"), remove = FALSE)

gene_FUN_002455 <- counts_by_event_filt %>% dplyr::filter(gene_id == "FUN_002455") %>% ungroup()

#filter this by the events that don't have enough timepoints (med and low treatments, and alt5' events)
gene_FUN_002455_filt <- gene_FUN_002455 %>% dplyr::filter(event_type != "alt_5prime") %>% dplyr::filter(treatment != "Medium") %>% dplyr::filter(treatment != "Low" )

# make it long format
gene_FUN_002455_filt_long <- gene_FUN_002455_filt %>% tidyr::pivot_longer(cols = c(flipped_dPSI,log2FoldChange), names_to = "test_type", values_to = "change")


# annot files ran emapper on annot from Nick to get GO terms
# read in the annotation files
annot_host <- read.delim("~/coral_genomes/other/Acer/Acerv.emapper.annotations", sep = "\t", header = T, skip = 4)

annot_host <- annot_host[1:20482,]

# split the transcript ID and gene ID
annot_host_clean_sep <- annot_host %>% tidyr::separate(X.query, sep = "-", into = c("gene_id", "trans_id"), remove = FALSE)

# filter to just one transcript per gene
annot_host_clean <- annot_host_clean_sep %>% dplyr::group_by(gene_id) %>% filter(trans_id == min(trans_id)) %>% distinct()


# annotate the shared events between AS and DE
shared_events_annot <- left_join(AS_DE_full_filt, annot_host_clean, by="gene_id")


# check the sign of the L2FC for gene expression
AS_DE_full_filt$sign_L2FC <- ifelse(AS_DE_full_filt$log2FoldChange < 0, "neg", "pos")

# check the sign for dPSI for AS
AS_DE_full_filt$sign_dPSI <- ifelse(AS_DE_full_filt$flipped_dPSI < 0, "neg", "pos")

# count number of events per gene
AS_DE_full_filt <- AS_DE_full_filt %>% dplyr::group_by(gene_id, cont) %>% dplyr::mutate(num_events = n()) %>% ungroup()

# check if sign matches for all events for each gene
AS_DE_full_filt <- AS_DE_full_filt %>% dplyr::group_by(gene_id, cont) %>% dplyr::mutate(num_event_signs = n_distinct(sign_dPSI)) %>% ungroup()

# check if sign matches between AS and DE for each event
AS_DE_full_filt$DE_AS_signs_match = ifelse(AS_DE_full_filt$sign_L2FC == AS_DE_full_filt$sign_dPSI, "yes", "no")


summary_by_event_class <- AS_DE_full_filt %>% dplyr::count(event_type, cont, DE_AS_signs_match)

# how much up vs down regulation do we see in ALL DEGs (not just those that overlap with AS)
lfc_DEG_col_int$sign_L2FC <- ifelse(lfc_DEG_col_int$log2FoldChange < 0, "neg", "pos")

# filter for just control contrasts
lfc_DGE_controls <- lfc_DEG_col_int %>% dplyr::filter(!str_detect(cont, "Field|Initial"))

# add a count of each gene (tells us how many contrasts the gene is present in)
lfc_DGE_controls <- lfc_DGE_controls %>% dplyr::add_count(gene_id)

# summarize 
summary_DGE <- lfc_DGE_controls %>% dplyr::count(cont, sign_L2FC)


# let's look at the types of relationships between number and sign of AS events and DGE
AS_DE_full_filt$num_event_class <- ifelse(AS_DE_full_filt$num_events > 1, "multi", "single")

# collapse these
AS_DE_full_filt$broad_class <- paste0(AS_DE_full_filt$num_event_signs, "_", AS_DE_full_filt$DE_AS_signs_match, "_", AS_DE_full_filt$num_event_class)

mixed <- AS_DE_full_filt %>% dplyr::filter(num_event_signs == "2")
agree_single <- AS_DE_full_filt %>% dplyr::filter(broad_class == "1_yes_single")
agree_multi <- AS_DE_full_filt %>% dplyr::filter(broad_class == "1_yes_multi")
disagree_single <- AS_DE_full_filt %>% dplyr::filter(broad_class == "1_no_single")
disagree_multi  <- AS_DE_full_filt %>% dplyr::filter(broad_class == "1_no_multi")


# prepare files to write out so we can do enrichment test on dell laptop

# write out this file
#write.csv(lfc_DGE_controls, file = "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/Acer_sumGene_DGE_sig_res_controls.csv", row.names = F, quote = F)

# also write out the annotated AS/DE overlap
#write.csv(AS_DE_full_filt, file = "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/Acer_sumGene_AS_DE_overlap_controls.csv", row.names = F, quote = F)


# make the df long format
AS_DE_full_filt_long <- AS_DE_full_filt %>% tidyr::pivot_longer(cols = c(sign_L2FC,sign_dPSI), names_to = "regulation_type", values_to = "sign")

AS_DE_full_filt_long <- AS_DE_full_filt_long %>% group_by(regulation_type, var1) %>% mutate(pos = sum(sign == "pos"), neg = sum(sign == "neg"))

AS_DE_full_filt_long2 <- AS_DE_full_filt_long %>% tidyr::pivot_longer(cols = c(pos,neg), names_to = "up_down", values_to = "count")

split_timepoint <- AS_DE_full_filt_long2 %>% tidyr::separate(var1, sep = "_", into = c("treatment", "timepoint"), remove = FALSE)

split_timepoint$treatment <- factor(split_timepoint$treatment, levels= c("Low", "Medium", "High"))
split_timepoint$up_down <- factor(split_timepoint$up_down, levels= c("pos", "neg"))

# make a bar chart comparing up and down regulated gene counts
up_down_plot <- ggplot(split_timepoint) + geom_bar(aes(x = regulation_type, y = count, fill = up_down), stat = "identity", position = "stack") + facet_grid(treatment ~ timepoint, scales = "free") + theme_bw()


summary_DGE_AS_overlap <- AS_DE_full_filt %>% dplyr::count(cont, sign_L2FC)

summary_joins <- full_join(summary_DGE, summary_DGE_AS_overlap, by = c("cont", "sign_L2FC"))
summary_joins <- summary_joins %>% dplyr::rename(all_DEGs = n.x, AS_DGEs = n.y)
summary_joins[is.na(summary_joins)] <- 0

summary_joins <- summary_joins %>% tidyr::pivot_longer(cols = c(all_DEGs,AS_DGEs), names_to = "type", values_to = "count")

summary_joins_split_vars <- summary_joins %>% tidyr::separate(cont, sep = "_vs_", into = c("var1", "var2"), remove = FALSE)

summary_joins_split_time <- summary_joins_split_vars %>% tidyr::separate(var1, sep = "_", into = c("treatment", "timepoint"), remove = FALSE)

summary_joins_split_time$treatment <- factor(summary_joins_split_time$treatment, levels= c("Low", "Medium", "High"))
summary_joins_split_time$sign_L2FC <- factor(summary_joins_split_time$sign_L2FC, levels= c("pos", "neg"))

up_down_plot_DGEs_all_overlap <- ggplot(summary_joins_split_time) + geom_bar(aes(x = timepoint, y = count, fill = sign_L2FC), stat = "identity", position = "stack") + facet_grid(type ~ treatment, scales = "free") + theme_bw()


