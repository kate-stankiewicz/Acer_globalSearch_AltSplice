# script for data exploration and figure plotting for A.cervicornis alternative splicing analysis

# Load libraries
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(rstatix)
library(readr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(fs)
library(ggbiplot)
library(ggfortify)
library(pheatmap)
library("RColorBrewer")
library(UpSetR)
library(ggVennDiagram)
library(viridis)

# set manuscript color schemes
#Treatments
Cont_col="#FED976"
Low_col = "#FD8D3C"
Med_col = "#E31A1C"
High_col = "#800026"

# timepoints
TF_col="#005100"
In_col="#BCFFBC"
T1_col = "#BB00BB"
T2_col = "#BFB2FF"
T3_col = "#FFF1FF"
T4_col = "#CCFAFF"
T5_col = "#2166AC"

# AS classes
alt_3_col ="#999999"
alt_5_col ="#E69F00"
ES_col = "#56B4E9"
IR_col = "#009E73"
multi_ES_col = "#D55E00"
mutex_col = "#CC79A7"

# read in the spladder output files for the sym
setwd("~/SplAdder/sym_spladder_res/array_spladder_out2")

sym_ldf <- list() 
listcsv <- dir(pattern = "_C3.confirmed.txt") 
for (k in 1:length(listcsv)){
  sym_ldf[[k]] <- read.table(file = listcsv[k], header = T)
}

# get only psi and gene name columns
sym_psi <- lapply(sym_ldf, function(x) x[, grepl('psi|gene|event' , names( x ) ) ])

# convert to one dataframe
sym_one_df <- do.call(rbind, sym_psi)


## now do the same for the the spladder output files for the coral host

setwd("~/SplAdder/host_spladder_res/array_spladder_out2")

host_ldf <- list() 
listcsv_h <- dir(pattern = "_C3.confirmed.txt") 
for (k in 1:length(listcsv_h)){
  host_ldf[[k]] <- read.table(file = listcsv_h[k], header = T)
}

# get only psi and gene name columns
host_psi <- lapply(host_ldf, function(x) x[, grepl('psi|gene|event' , names( x ) ) ])

# convert to one dataframe
host_one_df <- do.call(rbind, host_psi)


# count the number of events per sample

# for host
host_subset <- subset(host_one_df, select = -c(event_id, gene_name))
test_res_host <- as.data.frame(colSums(host_subset > 0.5, na.rm = T, ))

# for sym
sym_subset <- subset(sym_one_df, select = -c(event_id, gene_name))
test_res_sym <- as.data.frame(colSums(sym_subset > 0.5, na.rm = T, ))

# for each add in metadata for sample types

# make col name meaningful
colnames(test_res_host) <- "count"
colnames(test_res_sym) <- "count"

# add col name of sample IDs from rownames
test_res_host$sample_id <- row.names(test_res_host)
test_res_sym$sample_id <- row.names(test_res_sym)

# set the basename
test_res_host$sample_base <- gsub("host_|.sorted.psi|A|B", "", test_res_host$sample_id)
test_res_sym$sample_base <- gsub("sym_|.sorted.psi|A|B", "", test_res_sym$sample_id)

# read in the metadata file
metadata <- read.table("~/metadata/AC_clean_meta_full.tsv", header = T)

# rename the column to match
metadata <- metadata %>% dplyr::rename(sample_base = Novogene_R.)

# merge the data frames
merged_host <- merge(test_res_host, metadata, by = "sample_base")
merged_sym <- merge(test_res_sym, metadata, by = "sample_base")


## make some diagnostic plots of counts
#first set the appropriate order of the factors
merged_host$Treatment <- factor(merged_host$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))
merged_sym$Treatment <- factor(merged_sym$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))

merged_host$Timepoint <- factor(merged_host$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))
merged_sym$Timepoint <- factor(merged_sym$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

# next create the plots
plot_host <- ggplot(merged_host, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col))
plot_sym <- ggplot(merged_sym, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) 

# combine host and sym into one plot
# add a column to designation species
merged_host$species <- "host"
merged_sym$species <- "sym"

# merge the dataframes
combined_sym_host <- rbind(merged_host, merged_sym)

# create the joined plot
plot_combo <- ggplot(combined_sym_host, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(species ~ Treatment, scales = "free", space = "free_x")


# now lets consider event type
# split the event ID into type and number
split_ID_host <- host_one_df %>% tidyr::separate(event_id, sep = "\\.", into = c("event_type", "event_num"), remove = FALSE)
split_ID_sym <- sym_one_df %>% tidyr::separate(event_id, sep = "\\.", into = c("event_type", "event_num"), remove = FALSE)

# convert to long format
long_host <- split_ID_host %>% tidyr::pivot_longer(cols = c(-event_id, -event_type, -event_num, -gene_name), names_to = "sample_id", values_to = "psi")
long_sym <- split_ID_sym %>% tidyr::pivot_longer(cols = c(-event_id, -event_type, -event_num, -gene_name), names_to = "sample_id", values_to = "psi")

# count the number of events by event type
host_counts_type <- long_host %>% dplyr::group_by(event_type, sample_id) %>% dplyr::summarise(count = sum(psi > 0.5, na.rm = T))
sym_counts_type <- long_sym %>% dplyr::group_by(event_type, sample_id) %>% dplyr::summarise(count = sum(psi > 0.5, na.rm = T))

# join this new dataframe with the metadata again
# set the basename
host_counts_type$sample_base <- gsub("host_|.sorted.psi|A|B", "", host_counts_type$sample_id)
sym_counts_type$sample_base <- gsub("sym_|.sorted.psi|A|B", "", sym_counts_type$sample_id)

# merge with metadata
merged_host_type <- merge(host_counts_type, metadata, by = "sample_base")
merged_sym_type <- merge(sym_counts_type, metadata, by = "sample_base")


# Now generate plots showing counts by event type
#first set the appropriate order of the factors
merged_host_type$Treatment <- factor(merged_host_type$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))
merged_sym_type$Treatment <- factor(merged_sym_type$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))

merged_host_type$Timepoint <- factor(merged_host_type$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))
merged_sym_type$Timepoint <- factor(merged_sym_type$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

# create the plots
grid_host <- ggplot(merged_host_type, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "YlOrRd") + facet_grid(event_type ~., scales = "free", space = "free_x")
grid_sym <- ggplot(merged_sym_type, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "YlOrRd") + facet_grid(event_type ~ ., scales = "free", space = "free_x")

# look at overall counts per type
type_basic_host <- ggplot(merged_host_type, aes(x = reorder(event_type, -count), y = count, fill = event_type)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(alt_3_col, IR_col, ES_col, alt_5_col, mutex_col, multi_ES_col)) + labs(x = "Event type")
type_basic_sym <- ggplot(merged_sym_type, aes(x = event_type, y = count, fill = event_type)) + geom_boxplot() + theme_bw()


## calculate the number of NaNs per sample per species

# calculate the proportion for each sample in the host
host_NaN_prop <- as.data.frame(colMeans(is.na(host_one_df[,-c(1:2)])))

# change the column name
colnames(host_NaN_prop) <- "prop_NaN"

# add a column for sample ID
host_NaN_prop$sample_id <- rownames(host_NaN_prop)

# get the base name
host_NaN_prop$sample_base <- gsub("host_|.sorted.psi|A|B", "", host_NaN_prop$sample_id)

# now do the same in the sym
sym_NaN_prop <- as.data.frame(colMeans(is.na(sym_one_df[,-c(1:2)])))
colnames(sym_NaN_prop) <- "prop_NaN"
sym_NaN_prop$sample_id <- rownames(sym_NaN_prop)
sym_NaN_prop$sample_base <- gsub("sym_|.sorted.psi|A|B", "", sym_NaN_prop$sample_id)


# does NaN proportion correlate with count of events?
host_counts_perc_NaN <- merge(merged_host, host_NaN_prop, by = c("sample_base", "sample_id"))
sym_counts_perc_NaN <- merge(merged_sym, sym_NaN_prop, by = c("sample_base", "sample_id"))

# plot correlation between count and proportion NaN
host_counts_perc_NaN$prop_NaN <- as.numeric(host_counts_perc_NaN$prop_NaN)
host_counts_perc_NaN$Treatment <- factor(host_counts_perc_NaN$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))

host_na_count <- ggplot(host_counts_perc_NaN, aes(x = prop_NaN, y = count, col = Treatment)) + geom_point(pch=21, colour = "black", size = 2, aes(fill = Treatment)) + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(Treatment ~ Timepoint)

# repeat for the sym
sym_counts_perc_NaN$prop_NaN <- as.numeric(sym_counts_perc_NaN$prop_NaN)
sym_counts_perc_NaN$Treatment <- factor(sym_counts_perc_NaN$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))

sym_na_count <- ggplot(sym_counts_perc_NaN, aes(x = prop_NaN, y = count, col = Treatment)) + geom_point(pch=21, colour = "black", size = 2, aes(fill = Treatment)) + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(Treatment ~ Timepoint, scales = "free")


# now lets check boxplots of the distribution of NaN by sample type
# plot boxplots of mapping rate by sample
box_NaN_host <- ggplot(host_counts_perc_NaN, aes(x = Timepoint, y = prop_NaN, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free")

box_NaN_sym <- ggplot(sym_counts_perc_NaN, aes(x = Timepoint, y = prop_NaN, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free")

# next lets look at the total number of mapped reads per sample
# read in the mapped reads counts files
setwd("~/RNAseq_pipeline/count_mapped_reads")

file_names <- dir(pattern = "_total_mapped_reads_count.txt") 
counts_per_species <- do.call(rbind, lapply(file_names, function(x) cbind(read.table(x), name=strsplit(x,'\\.')[[1]][1])))

counts_per_species$name <- gsub("_total_mapped_reads_count", "", counts_per_species$name)

counts_per_species <- counts_per_species %>% dplyr::rename(total_mapped_reads = V1)

counts_per_species <- counts_per_species %>% tidyr::separate(name, sep = "\\_", into = c("species", "sample_base"), remove = FALSE)

counts_per_species <- counts_per_species %>% dplyr::select(total_mapped_reads, species, sample_base)

# split the host and sym
counts_per_species_host <- counts_per_species %>% dplyr::filter(species == "host")
counts_per_species_sym <- counts_per_species %>% dplyr::filter(species == "sym")

# calculate the average num of mapped reads for each species
counts_per_species_avg <- counts_per_species %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(average_mapped_reads = mean(total_mapped_reads, na.rm = TRUE)) %>%
  ungroup()

# join the counts per species with the existing meta file
all_meta_host <- merge(host_counts_perc_NaN, counts_per_species_host, by = c("sample_base", "species"))
all_meta_sym <- merge(sym_counts_perc_NaN, counts_per_species_sym, by = c("sample_base", "species"))

# calculate the count proportional to the total mapped reads
all_meta_host$count_per_total_reads_mapped <- all_meta_host$count / all_meta_host$total_mapped_reads
all_meta_sym$count_per_total_reads_mapped <- all_meta_sym$count / all_meta_sym$total_mapped_reads

# create a boxplots of these counts
box_total_number_mapped_host <- ggplot(all_meta_host, aes(x = Timepoint, y = total_mapped_reads, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free") 
box_total_number_mapped_sym <- ggplot(all_meta_sym, aes(x = Timepoint, y = total_mapped_reads, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free")


box_adj_count_host <- ggplot(all_meta_host, aes(x = Timepoint, y = count_per_total_reads_mapped, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free", space = "free") 
box_adj_count_sym <- ggplot(all_meta_sym, aes(x = Timepoint, y = count_per_total_reads_mapped, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free", space = "free")

# create correlation plots with total reads mapped per species
total_mapped_count_host <- ggplot(all_meta_host, aes(x = total_mapped_reads, y = prop_NaN, col = Treatment)) + geom_point(pch=21, colour = "black", size = 2, aes(fill = Treatment)) + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(Treatment ~ Timepoint, scales = "free")
total_mapped_count_sym <- ggplot(all_meta_sym, aes(x = total_mapped_reads, y = prop_NaN, col = Treatment)) + geom_point(pch=21, colour = "black", size = 2, aes(fill = Treatment)) + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(Treatment ~ Timepoint, scales = "free")


# combine host and sym together to get overall read map numbers
all_meta_both <- rbind(all_meta_host, all_meta_sym)

# sum the host and sym read counts for each sample
sum_read_counts <- all_meta_both %>% dplyr::group_by(sample_base) %>% dplyr::summarise(holobiont_total_mapped_reads = sum(total_mapped_reads), .groups = 'drop')

all_meta_both <- merge(all_meta_both, sum_read_counts, by = c("sample_base"))

# plot the total reads mapped for the holobiont
box_holo_mapped_reads <- ggplot(all_meta_both, aes(x = Timepoint, y = holobiont_total_mapped_reads, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free", space = "free") 

# split by host and sym in one plot 
box_total_number_mapped_host_sym_facet <- ggplot(all_meta_both, aes(x = Timepoint, y = total_mapped_reads, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(species ~ Treatment, scales = "free", space = "free") 

# plot the proportion NA for sym and host in one plot
box_propNA_host_sym_facet <- ggplot(all_meta_both, aes(x = Timepoint, y = prop_NaN, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(species ~ Treatment, scales = "free", space = "free") 


## next lets start global analysis of splicing events
# first we need to remove low quality samples
# let's try dropping samples with high NaN proportion

# set a threshold
threshold <- 0.6

# host
#remove the event id and gene id columns
host_just_psi <- host_one_df %>% dplyr::select(-c(event_id, gene_name))
host_drop_thresh_samps <- host_just_psi %>% dplyr::select(which(colMeans(is.na(.)) < threshold))

#sym
sym_just_psi <- sym_one_df %>% dplyr::select(-c(event_id, gene_name))
sym_drop_thresh_samps <- sym_just_psi %>% dplyr::select(which(colMeans(is.na(.)) < threshold))

# get list of samples passing the threshold
names_host <- as.data.frame(colnames(host_drop_thresh_samps))
names_sym <- as.data.frame(colnames(sym_drop_thresh_samps))

# change the column name
colnames(names_host) <- "sample_id"
colnames(names_sym) <- "sample_id"

# get the base sample name
names_host$sample_base <- gsub("host_|.sorted.psi", "", names_host$sample_id)
names_sym$sample_base <- gsub("sym_|.sorted.psi", "", names_sym$sample_id)

# merge with metadata
host_passing_samps_meta <- merge(names_host, metadata, by = "sample_base")
sym_passing_samps_meta <- merge(names_sym, metadata, by = "sample_base")

# check what is remaining
table(host_passing_samps_meta$treatment_timepoint)
table(sym_passing_samps_meta$treatment_timepoint)

table(metadata$treatment_timepoint)

# next check the variance in PSI across samples
variance_all <- cbind(var = apply(host_just_psi, 1, function(x) var(na.omit(x))), host_just_psi)

# remove rows with NAs
host_drop_na_events <- na.omit(host_just_psi)
host_drop_na_events_thresh <- na.omit(host_drop_thresh_samps)

sym_drop_na_events <- na.omit(sym_just_psi)
sym_drop_na_events_thresh <- na.omit(sym_drop_thresh_samps)

# calculate variance on the NA free events
variance_na_omit <- cbind(var = apply(host_drop_na_events, 1, function(x) var(na.omit(x))), host_drop_na_events)
sym_variance_na_omit <- cbind(var = apply(sym_drop_na_events_thresh, 1, function(x) var(na.omit(x))), sym_drop_na_events_thresh)


## run a PCA on all non-NA events
# transform and convert to a matrix
host_mat_t <- t(as.matrix(host_drop_na_events_thresh))
sym_mat_t <- t(as.matrix(sym_drop_na_events_thresh))

# run pca
host_pca <- prcomp(host_mat_t)
sym_pca <- prcomp(sym_mat_t)

# get the meta information incorporated 
host_meta_pca_samples <- as.data.frame(rownames(host_mat_t))
colnames(host_meta_pca_samples) <- "sample_id"
host_meta_pca_samples$sample_base <- gsub("host_|.sorted.psi", "", host_meta_pca_samples$sample_id)

host_meta_pca <- merge(host_meta_pca_samples, metadata, by = c("sample_base"))

# plot the pca
host_plot <- ggbiplot::ggbiplot(host_pca,ellipse=TRUE, var.axes=FALSE) + geom_point(size = 4, aes(shape = as.factor(host_meta_pca$Timepoint), color = as.factor(host_meta_pca$Treatment))) + theme_bw() + scale_color_manual(breaks=c('Field', 'Initial', 'Control', 'Low', 'Medium', 'High'), values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + scale_shape_manual(breaks=c('TF', 'T0', 'T1', 'T2', 'T3', 'T4', 'T5'), values=c(4,8,15,16,17,18,12)) + labs(color = "Treatment", shape = "Timepoint")

# do the same thing for the symbiont
# get the meta information incorporated 
sym_meta_pca_samples <- as.data.frame(rownames(sym_mat_t))
colnames(sym_meta_pca_samples) <- "sample_id"
sym_meta_pca_samples$sample_base <- gsub("sym_|.sorted.psi", "", sym_meta_pca_samples$sample_id)

sym_meta_pca <- merge(sym_meta_pca_samples, metadata, by = c("sample_base"))

# plot the pca
sym_plot <- ggbiplot::ggbiplot(sym_pca,ellipse=TRUE, var.axes=FALSE) + geom_point(size = 3, aes(shape = as.factor(sym_meta_pca$Timepoint), color = as.factor(sym_meta_pca$Treatment))) + theme_bw() + scale_color_manual(breaks=c('Field', 'Initial', 'Control', 'Low', 'Medium', 'High'), values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + scale_shape_manual(breaks=c('TF', 'T0', 'T1', 'T2', 'T3', 'T4', 'T5'), values=c(4,8,15,16,17,18,12)) + labs(color = "Treatment", shape = "Timepoint")


## check if there's more or less alt usage based on temp
## create a heatmap of psi values for host

# get just samples with NAs above threshold and filter original data frame for these
host_filt_thresh <- split_ID_host %>% dplyr::select(c(event_id, event_type, event_num, gene_name, c(names_host$sample_id)))
sym_filt_thresh <- split_ID_sym %>% dplyr::select(c(event_id, event_type, event_num, gene_name, c(names_sym$sample_id)))

# remove events with NAs from the df that includes event type information
host_no_na_event_ID <- na.omit(host_filt_thresh)
sym_no_na_event_ID <- na.omit(sym_filt_thresh)

# remove rows where all PSI values across samples are exactly the same
keep <- apply(host_no_na_event_ID[5:ncol(host_no_na_event_ID)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
host_rem_psi_ident <- host_no_na_event_ID[keep, ]

keep_sym <- apply(sym_no_na_event_ID[5:ncol(sym_no_na_event_ID)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
sym_rem_psi_ident <- sym_no_na_event_ID[keep_sym, ]

# calculate variance on remaining rows to identify those with little variance between all samples
low_var_estimates <- cbind(var = apply(host_rem_psi_ident[5:ncol(host_rem_psi_ident)], 1, function(x) var(na.omit(x))), host_rem_psi_ident)
low_var_estimates_sym <- cbind(var = apply(sym_rem_psi_ident[5:ncol(sym_rem_psi_ident)], 1, function(x) var(na.omit(x))), sym_rem_psi_ident)

# remove rows with variance below a given threshold
low_var_removed <- low_var_estimates %>% dplyr::filter(var > 0.04)
low_var_removed_sym <- low_var_estimates_sym %>% dplyr::filter(var > 0.002)

# create the heatmap of the highest variance events across all samples
# filter for just the psi columns
psi_only_low_var <- low_var_removed %>% dplyr::select(!(var:gene_name))
psi_only_low_var_sym <- low_var_removed_sym %>% dplyr::select(!(var:gene_name))

# set the event ID as the rowname
rownames(psi_only_low_var) <- low_var_removed$event_id
rownames(psi_only_low_var_sym) <- low_var_removed_sym$event_id

# reformat as a matrix
mat_low_var_rem <- as.matrix(psi_only_low_var)
mat_low_var_rem_sym <- as.matrix(psi_only_low_var_sym)

# create the row and column annotations for the files
# get the sample order
sample_order_mat <- as.data.frame(colnames(mat_low_var_rem))
colnames(sample_order_mat) <- "sample_id"
sample_order_mat$sample_base <- gsub("host_|.sorted.psi", "", sample_order_mat$sample_id)

sample_order_mat_sym <- as.data.frame(colnames(mat_low_var_rem_sym))
colnames(sample_order_mat_sym) <- "sample_id"
sample_order_mat_sym$sample_base <- gsub("sym_|.sorted.psi", "", sample_order_mat_sym$sample_id)

# merge with metadata information
sample_order_meta <- merge(sample_order_mat, metadata, by = "sample_base")
sample_order_meta_sym <- merge(sample_order_mat_sym, metadata, by = "sample_base")

# generate the col annotation
col_annot <- sample_order_meta %>% dplyr::select(Timepoint, Treatment)
rownames(col_annot) <- sample_order_meta$sample_id

col_annot_sym <- sample_order_meta_sym %>% dplyr::select(Timepoint, Treatment)
rownames(col_annot_sym) <- sample_order_meta_sym$sample_id

# generate the row annotation
row_annot <- data.frame(low_var_removed$event_type)
rownames(row_annot) <- low_var_removed$event_id
colnames(row_annot) <- c("Event_Type")

row_annot_sym <- data.frame(low_var_removed_sym$event_type)
rownames(row_annot_sym) <- low_var_removed_sym$event_id
colnames(row_annot_sym) <- c("Event_Type")

# set the colors for the heatmap
annot_colors = list(
  Treatment = c(Field=TF_col, Initial=In_col, Control=Cont_col, Low = Low_col, Medium = Med_col, High = High_col),
  Timepoint = c(TF = TF_col, T0 = In_col, T1 = T1_col, T2 = T2_col, T3 = T3_col, T4 = T4_col, T5 = T5_col),
  #Genotype = c("07" = "pink", "31" = "purple", "34" = "dark green", "41" = "maroon", "48" = "black", "50" = "grey", "62" = "yellow", "CM5" = "white"),
  Event_Type  = c(alt_3prime=alt_3_col, alt_5prime=alt_5_col, exon_skip = ES_col, intron_retention = IR_col, mult_exon_skip = multi_ES_col, mutex_exons = mutex_col))

# create a basic heatmap
heat_low_var_rem <- pheatmap(mat_low_var_rem, cluster_cols = F, annotation_col = col_annot, annotation_row = row_annot, annotation_colors = annot_colors, show_rownames = F, show_colnames = F, color = viridis(10))

heat_low_var_rem_sym <- pheatmap(mat_low_var_rem_sym, cluster_cols = F, annotation_col = col_annot_sym, annotation_row = row_annot_sym, annotation_colors = annot_colors, show_rownames = F, show_colnames = F, color = viridis(10))


## now let's look at the PSI distribution for each group for these high variance events
long_host_no_na_ident <- low_var_removed %>% tidyr::gather(sample_id, psi, -c("var", "event_id", "event_type", "event_num", "gene_name"))
long_sym_no_na_ident <- low_var_removed_sym %>% tidyr::gather(sample_id, psi, -c("var", "event_id", "event_type", "event_num", "gene_name"))

# add base sample name
long_host_no_na_ident$sample_base <- gsub("host_|.sorted.psi", "", long_host_no_na_ident$sample_id)
long_sym_no_na_ident$sample_base <- gsub("sym_|.sorted.psi", "", long_sym_no_na_ident$sample_id)

# add the metadata 
long_host_no_na_ident <- dplyr::left_join(long_host_no_na_ident, metadata, by = "sample_base")
long_sym_no_na_ident <- dplyr::left_join(long_sym_no_na_ident, metadata, by = "sample_base")

# plot this distribution
long_host_no_na_ident$Treatment <- factor(long_host_no_na_ident$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))
long_host_no_na_ident$Timepoint <- factor(long_host_no_na_ident$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

long_sym_no_na_ident$Treatment <- factor(long_sym_no_na_ident$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium"))
long_sym_no_na_ident$Timepoint <- factor(long_sym_no_na_ident$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

plot_psi_box <- ggplot(long_host_no_na_ident, aes(fill=Timepoint, y=psi, x=Treatment)) + geom_boxplot(position=position_dodge2(preserve = "single")) + theme_bw() + scale_fill_manual(breaks=c('TF', 'T0', 'T1', 'T2', 'T3', 'T4', 'T5'), values=c(TF_col, In_col, T1_col, T2_col, T3_col, T4_col, T5_col))

plot_psi_dens <- ggplot(long_host_no_na_ident, aes(fill = factor(Timepoint, levels = c('T5', 'T4', 'T3', 'T2', 'T1', 'T0', 'TF')), x=psi)) + geom_density(alpha = 0.6) + theme_bw() + facet_grid(. ~ Treatment) + scale_fill_manual(values=c(T5_col, T4_col, T3_col,T2_col, T1_col, In_col, TF_col)) + guides(fill = guide_legend(reverse = TRUE)) + labs(fill = "Timepoint")

plot_psi_box_sym <- ggplot(long_sym_no_na_ident, aes(fill=Timepoint, y=psi, x=Treatment)) + geom_boxplot(position=position_dodge2(preserve = "single")) + theme_bw() + scale_fill_manual(breaks=c('TF', 'T0', 'T1', 'T2', 'T3', 'T4', 'T5'), values=c(TF_col, In_col, T1_col, T2_col, T3_col, T4_col, T5_col))

plot_psi_dens_sym <- ggplot(long_sym_no_na_ident, aes(fill = factor(Timepoint, levels = c('T5', 'T4', 'T3', 'T2', 'T1', 'T0', 'TF')), x=psi)) + geom_density(alpha = 0.5) + theme_bw() + facet_grid(. ~ Treatment) + scale_fill_manual(values=c(T5_col, T4_col, T3_col,T2_col, T1_col, In_col, TF_col)) + guides(fill = guide_legend(reverse = TRUE)) + labs(fill = "Timepoint")


# Next check the total number of each event type
table(row_annot)

# get the ones with more splicing in the hottest treatments (indexes selected after viewing heatmap)
index_up <- as.vector(heat_low_var_rem$tree_row$order[2:55])
heat_up_events <- heat_low_var_rem$tree_row$labels[index_up]
heat_up_event_type <- row_annot %>% dplyr::filter(row.names(row_annot) %in% heat_up_events)
table(heat_up_event_type$Event_Type)

#get the ones with LESS splicing in the hotest treatments
index_down <- as.vector(heat_low_var_rem$tree_row$order[101:106])
heat_down_events <- heat_low_var_rem$tree_row$labels[index_down]
heat_down_event_type <- row_annot %>% dplyr::filter(row.names(row_annot) %in% heat_down_events)
table(heat_down_event_type$Event_Type)


## count the events by type for each sample after removing samples with high NaN and then removing all events with any NaN
# first convert to long format
long_host_thresh_no_NA <- host_no_na_event_ID %>% tidyr::pivot_longer(cols = c(-event_id, -event_type, -event_num, -gene_name), names_to = "sample_id", values_to = "psi")

# count by group
long_host_thresh_no_NA_counts <- long_host_thresh_no_NA %>% dplyr::group_by(event_type, sample_id) %>% dplyr::summarise(count = sum(psi > 0.5, na.rm = T)) %>% ungroup()

# set the base name and add in metainformation 
long_host_thresh_no_NA_counts$sample_base <- gsub("host_|.sorted.psi|A|B", "", long_host_thresh_no_NA_counts$sample_id)

lh_thresh_noNA_countmeta <- merge(long_host_thresh_no_NA_counts, metadata, by = "sample_base")

# get count of number of samples for each category
lh_tnoNA_countMeta_numSamp <- lh_thresh_noNA_countmeta %>% dplyr::add_count(treatment_timepoint)

# divide count by number of samples
lh_tnoNA_countMeta_numSamp$adj_count <- lh_tnoNA_countMeta_numSamp$count / lh_tnoNA_countMeta_numSamp$n

# plot these counts by sample type
lh_tnoNA_countMeta_numSamp$Timepoint <- factor(lh_tnoNA_countMeta_numSamp$Timepoint, levels= c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))
lh_tnoNA_countMeta_numSamp$Treatment <- factor(lh_tnoNA_countMeta_numSamp$Treatment, levels= c("Field", "Initial", "Control", "Low", "Medium", "High"))

box_adj_counts_threshNoNA_host <- ggplot(lh_tnoNA_countMeta_numSamp, aes(x = Timepoint, y = count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "YlOrRd") + facet_grid(event_type ~ Treatment, scales = "free", space = "free") 


# Next we will examine the differentially spliced events identified by Spladder test mode
# the output from Spladder test mode was parsed by a different script to collect all the events in merged files

# read in testing files
setwd("~/SplAdder/host_spladder_res/array_spladder_out2/all_events_parsed")

ldf <- list() 
listcsv <- dir(pattern = ".tsv") 
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.table(file = listcsv[k], header = T)
}

# combine into one DF
one_df <- do.call(rbind, ldf)

#filter for only significant events
sig_all <- one_df %>% dplyr::filter(p_val_adj <= 0.05)

split_ID <- sig_all %>% tidyr::separate(event_id, sep = "\\.", into = c("event_type", "event_num"), remove = FALSE)

# write out the sig events to a file
write.table(sig_all, "~/results/Acer_Sfit_runs/analysis/Acer_all_sig_events_combined.tsv", quote = F, row.names = F)

# repeat for sym
# read in testing files
setwd("~/SplAdder/sym_spladder_res/array_spladder_out2/all_events_parsed")

ldf_symt <- list() 
listcsv_symt <- dir(pattern = ".tsv") 
for (k in 1:length(listcsv_symt)){
  ldf_symt[[k]] <- read.table(file = listcsv_symt[k], header = T)
}

# combine into one DF
one_df_symt <- do.call(rbind, ldf_symt)

#filter for only significant events
sig_all_sym <- one_df_symt %>% dplyr::filter(p_val_adj <= 0.05)

split_ID_symt <- sig_all_sym %>% tidyr::separate(event_id, sep = "\\.", into = c("event_type", "event_num"), remove = FALSE)

# write out the sig events to a file
write.table(sig_all_sym, "~/results/Acer_Sfit_runs/analysis/Sfit_all_sig_events_combined.tsv", quote = F, row.names = F)


# plot the significant results per contrast
event_counts <- split_ID %>% dplyr::count(event_type)
event_counts_sym <- split_ID_symt %>% dplyr::count(event_type)

# split the vrs column
vrs_split <- split_ID %>% tidyr::separate(vrs, sep = "_", into = c("var1", "timepoint1", "vs", "var2", "timepoint2"), remove = FALSE)
vrs_split_sym <- split_ID_symt %>% tidyr::separate(vrs, sep = "_", into = c("var1", "timepoint1", "vs", "var2", "timepoint2"), remove = FALSE)

# put vars together to create a vrs column independent of timepoint 
vrs_split$treatment_contrast <- paste0(vrs_split$var1,"_vs_",vrs_split$var2)
vrs_split_sym$treatment_contrast <- paste0(vrs_split_sym$var1,"_vs_",vrs_split_sym$var2)

# filter for just the contrast types with heat treatments 
no_field_control_initial <- vrs_split %>% dplyr::filter(!treatment_contrast %in% c("Field_vs_Initial","Field_vs_Control","Initial_vs_Control"))
no_field_control_initial_sym <- vrs_split_sym %>% dplyr::filter(!treatment_contrast %in% c("Field_vs_Initial","Field_vs_Control","Initial_vs_Control"))

just_treat_v_cont <- vrs_split %>% dplyr::filter(!(str_detect(treatment_contrast, "Field|Initial")))

# set the order of the variables
no_field_control_initial$var1 <- factor(no_field_control_initial$var1, levels = c("Field", "Initial", "Control"))
no_field_control_initial$var2 <- factor(no_field_control_initial$var2, levels = c("Low", "Medium", "High"))

just_treat_v_cont$var2 <- factor(just_treat_v_cont$var2, levels = c("Low", "Medium", "High"))

no_field_control_initial_sym$var1 <- factor(no_field_control_initial_sym$var1, levels = c("Field", "Initial", "Control"))
no_field_control_initial_sym$var2 <- factor(no_field_control_initial_sym$var2, levels = c("Low", "Medium", "High"))

# plot the counts
grid_t <- ggplot(no_field_control_initial, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(var2~var1, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_just_v_cont <- ggplot(just_treat_v_cont, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(.~var2, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_sym <- ggplot(no_field_control_initial_sym, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(var2~var1, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

# now pull out the CBASS impact contrasts
field_control_initial <- vrs_split %>% dplyr::filter(treatment_contrast %in% c("Field_vs_Initial","Field_vs_Control","Initial_vs_Control"))
field_control_initial_sym <- vrs_split_sym %>% dplyr::filter(treatment_contrast %in% c("Field_vs_Initial","Field_vs_Control","Initial_vs_Control"))

# plot the counts
grid_t_cbass <- ggplot(field_control_initial, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(.~treatment_contrast, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_cbass_sym <- ggplot(field_control_initial_sym, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(.~treatment_contrast, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

# filter for dPSI > 0.3
dpsi_filtered_swapped_vrs_split <- no_field_control_initial %>% dplyr::filter(abs(dPSI) > 0.3)
dpsi_filtered_field_control_initial <-field_control_initial %>% dplyr::filter(abs(dPSI) > 0.3)

dpsi_filtered_just_v_cont <- just_treat_v_cont %>% dplyr::filter(abs(dPSI) > 0.3)


dpsi_filtered_swapped_vrs_split_sym <- no_field_control_initial_sym %>% dplyr::filter(abs(dPSI) > 0.3)
dpsi_filtered_field_control_initial_sym <-field_control_initial_sym %>% dplyr::filter(abs(dPSI) > 0.3)

# create the same plots on the dpsi filtered dataframes
grid_t_filt<- ggplot(dpsi_filtered_swapped_vrs_split, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "stack") + theme_bw() + facet_grid(var2~var1, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_filt_just_v_cont <- ggplot(dpsi_filtered_just_v_cont, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge", colour = "black") + theme_bw() + facet_grid(.~var2, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_filt_just_v_cont_rel_prop <- ggplot(dpsi_filtered_just_v_cont, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "fill", colour = "black") + theme_bw() + labs(x = "Timepoint") + facet_grid(.~var2) + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col)) + 
  scale_y_continuous(labels = scales::percent) + 
  ylab("Percent (%)") 

grid_t_filt_sym <- ggplot(dpsi_filtered_swapped_vrs_split_sym, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(var2~var1, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_cbass_filt <- ggplot(dpsi_filtered_field_control_initial, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "stack", colour = "black") + theme_bw() + facet_grid(.~treatment_contrast, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))

grid_t_cbass_filt_sym <- ggplot(dpsi_filtered_field_control_initial_sym, aes(fill=event_type, x = as.factor(timepoint2))) + geom_bar(position = "dodge") + theme_bw() + facet_grid(.~treatment_contrast, scales = "free", space = "free") + labs(x = "Timepoint", y = "Count") + scale_fill_manual(values = c(alt_3_col, alt_5_col, ES_col, IR_col, multi_ES_col, mutex_col))


# join the test result file with psi file by event id
# rename the gene_name column 
split_ID_host_rename <- host_filt_thresh %>% dplyr::rename(gene_id = gene_name)
joined_treatment_cont_psi <- dplyr::left_join(dpsi_filtered_swapped_vrs_split, split_ID_host_rename, by = c('event_id', 'event_type', 'event_num', 'gene_id'))
joined_cbass_psi <- dplyr::left_join(dpsi_filtered_field_control_initial, split_ID_host_rename, by = c('event_id', 'event_type', 'event_num', 'gene_id'))

split_ID_sym_rename <- sym_filt_thresh %>% dplyr::rename(gene_id = gene_name)
joined_treatment_cont_psi_sym <- dplyr::left_join(dpsi_filtered_swapped_vrs_split_sym, split_ID_sym_rename, by = c('event_id', 'event_type', 'event_num', 'gene_id'))
joined_cbass_psi_sym <- dplyr::left_join(dpsi_filtered_field_control_initial_sym, split_ID_sym_rename, by = c('event_id', 'event_type', 'event_num', 'gene_id'))

# how many are left if we get rid of all events where ANY sample has missing information?
joined_treat_cont_psi_no_na <- na.omit(joined_treatment_cont_psi)
joined_cbass_psi_no_na <- na.omit(joined_cbass_psi)


joined_treat_cont_psi_no_na_sym <- na.omit(joined_treatment_cont_psi_sym)
joined_cbass_psi_no_na_sym <- na.omit(joined_cbass_psi_sym)


## check for overlap of events across contrasts
# filter dPSI values 
vrs_dpsi0.3 <- vrs_split %>% dplyr::filter(abs(dPSI) > 0.3)

vrs_dpsi0.3_sym <- vrs_split_sym %>% dplyr::filter(abs(dPSI) > 0.3)

# generate a binary matrix of presence absence for each contrast
event_id_contrast <- vrs_dpsi0.3 %>% dplyr::select(event_id, vrs)
wide_dpsi0.3 <- pivot_wider(event_id_contrast, names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

event_id_contrast_sym <- vrs_dpsi0.3_sym %>% dplyr::select(event_id, vrs)
wide_dpsi0.3_sym <- pivot_wider(event_id_contrast_sym, names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

# create a complex upset plot of overlaps of event ids between contrasts
upsetR_events <- as.data.frame(select(wide_dpsi0.3, -c(event_id)))
list_all <- colnames(upsetR_events)
high <- grep('High', list_all, ignore.case = F, value = T)
med <- grep('High', list_all, ignore.case = F, value = T)
events_upsetR_plot <- UpSetR::upset(upsetR_events, nsets = 56, order.by = "freq", nintersects = 40) #, query.legend = "top", queries = list(list(query = intersects, params = list(high), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list(med), color = "#E31A1C", active = T, query.name="Medium")))

upsetR_events_sym <- as.data.frame(select(wide_dpsi0.3_sym, -c(event_id)))
list_all_sym <- colnames(upsetR_events_sym)
high_sym <- grep('High', list_all_sym, ignore.case = F, value = T)
med_sym <- grep('High', list_all_sym, ignore.case = F, value = T)
events_upsetR_plot_sym <- UpSetR::upset(upsetR_events_sym, nsets = 56, order.by = "freq", nintersects = 40) #, query.legend = "top", queries = list(list(query = intersects, params = list(high), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list(med), color = "#E31A1C", active = T, query.name="Medium")))

# now look at the overlap of genes
gene_id_contrast <- unique(vrs_dpsi0.3[ , c("gene_id", "vrs")])  
wide_dpsi0.3_genes <- pivot_wider(gene_id_contrast, names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

upsetR_genes <- as.data.frame(select(wide_dpsi0.3_genes, -c(gene_id)))
genes_upsetR_plot <- UpSetR::upset(upsetR_genes, nsets = 56, order.by = "freq", nintersects = 20) #, query.legend = "top", queries = list(list(query = intersects, params = list(high), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list(med), color = "#E31A1C", active = T, query.name="Medium")))

gene_id_contrast_sym <- unique(vrs_dpsi0.3_sym[ , c("gene_id", "vrs")])  
wide_dpsi0.3_genes_sym <- pivot_wider(gene_id_contrast_sym, names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

upsetR_genes_sym <- as.data.frame(select(wide_dpsi0.3_genes_sym, -c(gene_id)))
genes_upsetR_plot_sym <- UpSetR::upset(upsetR_genes_sym, nsets = 56, order.by = "freq", nintersects = 20) #, query.legend = "top", queries = list(list(query = intersects, params = list(high), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list(med), color = "#E31A1C", active = T, query.name="Medium")))

# look at intersection regardless of timepoint
event_id_contrast_basic <- vrs_dpsi0.3 %>% dplyr::select(event_id, treatment_contrast)
wide_dpsi0.3_basic <- pivot_wider(event_id_contrast_basic, names_from = "treatment_contrast", values_from = 'treatment_contrast', values_fill = 0, values_fn = function(x) 1)

upsetR_events_basic <- as.data.frame(select(wide_dpsi0.3_basic, -c(event_id)))
events_basic_upsetR_plot <- UpSetR::upset(upsetR_events_basic, nsets = 12, order.by = "freq", nintersects = 20, query.legend = "top", queries = list(list(query = intersects, params = list('Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium'), color = "#E31A1C", active = T, query.name="Medium"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium', 'Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "blue", active = T, query.name="Both")))

event_id_contrast_basic_sym <- vrs_dpsi0.3_sym %>% dplyr::select(event_id, treatment_contrast)
wide_dpsi0.3_basic_sym <- pivot_wider(event_id_contrast_basic_sym, names_from = "treatment_contrast", values_from = 'treatment_contrast', values_fill = 0, values_fn = function(x) 1)

upsetR_events_basic_sym <- as.data.frame(select(wide_dpsi0.3_basic_sym, -c(event_id)))
events_basic_upsetR_plot_sym <- UpSetR::upset(upsetR_events_basic_sym, nsets = 12, order.by = "freq", nintersects = 30, query.legend = "top", queries = list(list(query = intersects, params = list('Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium'), color = "#E31A1C", active = T, query.name="Medium")))

# check this for genes now instead
gene_id_contrast_basic <- unique(vrs_dpsi0.3[ , c("gene_id", "treatment_contrast")])  
wide_dpsi0.3_genes_basic <- pivot_wider(gene_id_contrast_basic, names_from = "treatment_contrast", values_from = 'treatment_contrast', values_fill = 0, values_fn = function(x) 1)

upsetR_genes_basic <- as.data.frame(select(wide_dpsi0.3_genes_basic, -c(gene_id)))
genes_basic_upsetR_plot <- UpSetR::upset(upsetR_genes_basic, nsets = 12, order.by = "freq", nintersects = 20, query.legend = "top", queries = list(list(query = intersects, params = list('Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium'), color = "#E31A1C", active = T, query.name="Medium"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium', 'Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "blue", active = T, query.name="Both")))

gene_id_contrast_basic_sym <- unique(vrs_dpsi0.3_sym[ , c("gene_id", "treatment_contrast")])  
wide_dpsi0.3_genes_basic_sym <- pivot_wider(gene_id_contrast_basic_sym, names_from = "treatment_contrast", values_from = 'treatment_contrast', values_fill = 0, values_fn = function(x) 1)

upsetR_genes_basic_sym <- as.data.frame(select(wide_dpsi0.3_genes_basic_sym, -c(gene_id)))
genes_basic_upsetR_plot_sym <- UpSetR::upset(upsetR_genes_basic_sym, nsets = 12, order.by = "freq", nintersects = 20, query.legend = "top", queries = list(list(query = intersects, params = list('Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium'), color = "#E31A1C", active = T, query.name="Medium"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium', 'Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "blue", active = T, query.name="Both")))


## summary information

# filter for events where at least one sample has a psi above 0.5
host_no_na_psi_0.5 <- host_no_na_event_ID %>% dplyr::filter_at(vars(starts_with("host")), any_vars(. >= 0.5))

sym_no_na_psi_0.5 <- sym_no_na_event_ID %>% dplyr::filter_at(vars(starts_with("sym")), any_vars(. >= 0.5))

# how many genes are involved in alt splice events (psi >= 0.5 in at least one sample)?
length(unique(host_no_na_psi_0.5$gene_name))

length(unique(sym_no_na_psi_0.5$gene_name))

# check this proportional relative to total genes in genome
length(unique(host_no_na_psi_0.5$gene_name))/29806 * 100

length(unique(sym_no_na_psi_0.5$gene_name))/37621 * 100

# get a list of these genes with splice events 
gene_list_alt_spliced <- unique(host_no_na_psi_0.5$gene_name)

gene_list_alt_spliced_sym <- unique(sym_no_na_psi_0.5$gene_name)

# get this same information for the unfiltered files
# for the no NA file
length(unique(host_no_na_event_ID$gene_name))
length(unique(host_no_na_event_ID$gene_name))/29806 * 100

length(unique(host_no_na_event_ID$gene_name))
length(unique(sym_no_na_event_ID$gene_name))/37621 * 100

# and for the full file
length(unique(split_ID_host$gene_name))
length(unique(split_ID_host$gene_name))/29806 * 100

all_psi_0.5 <- split_ID_host %>% dplyr::filter_at(vars(starts_with("host")), any_vars(. >= 0.5))
length(unique(all_psi_0.5$gene_name))
length(unique(all_psi_0.5$gene_name))/29806 * 100

length(unique(split_ID_sym$gene_name))
length(unique(split_ID_sym$gene_name))/37621 * 100

all_psi_0.5_sym <- split_ID_sym %>% dplyr::filter_at(vars(starts_with("sym")), any_vars(. >= 0.5))
length(unique(all_psi_0.5_sym$gene_name))
length(unique(all_psi_0.5_sym$gene_name))/37621 * 100


# now count the number of splice events for each
# with at least one sample >= psi 0.5
## all
length(all_psi_0.5$event_id)

length(all_psi_0.5_sym$event_id)

## no NAs
length(host_no_na_psi_0.5$event_id)

length(sym_no_na_psi_0.5$event_id)

# all
## all
length(split_ID_host$event_id)

length(split_ID_sym$event_id)

## no NAs
length(host_no_na_event_ID$event_id)

length(sym_no_na_event_ID$event_id)


# look at the gene level for the top variance events
per_gene_counts_low_var_psi <- long_host_no_na_ident %>% dplyr::group_by(gene_name, sample_base) %>% dplyr::summarise(count = sum(psi > 0.5, na.rm = T))

# join back with the meta information
per_gene_counts_low_var_psi_meta <- dplyr::left_join(per_gene_counts_low_var_psi, metadata, by = "sample_base")

# convert back to wide format
wide_per_gene_low_var_psi <- per_gene_counts_low_var_psi %>% tidyr::pivot_wider(names_from = "sample_base", values_from = "count")

# convert to a matrix and add in gene names as row names
gene_counts_just_counts <- subset(wide_per_gene_low_var_psi, select = -c(gene_name))
gene_counts_just_counts <- sapply(gene_counts_just_counts, function(x) as.numeric(as.character(x)))
mat_gene_counts_low_var <- as.matrix(gene_counts_just_counts)
rownames(mat_gene_counts_low_var) <- wide_per_gene_low_var_psi$gene_name

mat_gene_counts_low_var_scale <- scale(mat_gene_counts_low_var)

# generate the column annotations
col_annot_genes <- metadata %>% dplyr::select(Timepoint, Treatment, Genotype)
rownames(col_annot_genes) <- metadata$sample_base

# set the colors for the heatmap
annot_colors_genes = list(
  Treatment = c(Field="blue", Initial="green", Control="#FED976", Low = "#FD8D3C", Medium = "#E31A1C", High = "#800026"),
  Timepoint = c(TF = "blue", T0 = "green", T1 = "#B2182B", T2 = "#D6604D", T3 = "#FDDBC7", T4 = "#92C5DE", T5 = "#2166AC"),
  Genotype = c("07" = "pink", "31" = "purple", "34" = "dark green", "41" = "maroon", "48" = "black", "50" = "grey", "62" = "yellow", "CM5" = "white"))

# basic heatmap
heat_low_var_rem_genes <- pheatmap(mat_gene_counts_low_var_scale, scale = "none", annotation_col = col_annot_genes, annotation_colors = annot_colors_genes)

# check distribution of the counts
per_gene_counts_low_var_psi_meta$Timepoint <- factor(per_gene_counts_low_var_psi_meta$Timepoint, levels= c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))
per_gene_counts_low_var_psi_meta$Treatment <- factor(per_gene_counts_low_var_psi_meta$Treatment, levels= c("Field", "Initial", "Control", "Low", "Medium", "High"))

plot_gene_counts_dens <- ggplot(per_gene_counts_low_var_psi_meta, aes(fill = Timepoint, x=count)) + geom_bar(position = "stack") + theme_bw() + facet_grid(. ~ Treatment) + scale_fill_manual(values=c("blue", "green", "#B2182B", "#D6604D", "#FDDBC7", "#92C5DE", "#2166AC")) + labs(fill = "Timepoint")



# plot the contrast counts showing number of genes and number of events side by side

# count the number of unique events for each contrast
event_counts_per_cont <- dpsi_filtered_swapped_vrs_split %>% dplyr::group_by(vrs, timepoint2, var1, var2) %>% dplyr::summarise(num_events = n_distinct(event_id))

# count the number of unique genes involved for each contrast
gene_counts_per_cont <- dpsi_filtered_swapped_vrs_split %>% dplyr::group_by(vrs) %>% dplyr::summarise(num_genes = n_distinct(gene_id))

# merge the two 
gene_event_counts_per_cont <- merge(event_counts_per_cont, gene_counts_per_cont, by = "vrs")

# make it long format
long_gene_event_counts_per_cont <- gene_event_counts_per_cont %>% tidyr::pivot_longer(cols = c(-vrs, -timepoint2, -var1, -var2), names_to = "level", values_to = "count")

# plot this
grid_gene_events_cont <- ggplot(long_gene_event_counts_per_cont, aes(fill=level, x = as.factor(timepoint2), y = count)) + geom_bar(stat = "identity", position = "dodge") + theme_bw() + facet_grid(var2~var1) + labs(x = "Timepoint", y = "Count")


## functional information

# read in the tables of significantly enriched terms for each contrast (output after running script for clusterprofiler)
setwd("~/ClusterProfiler_GOEnrichment_results")

enrich_ldf <- list() 
enrich_listcsv <- dir(pattern = ".csv") 
for (k in 1:length(enrich_listcsv)){
  enrich_ldf[[k]] <- read.delim(file = enrich_listcsv[k], header = T, sep = "\t")
}

names(enrich_ldf) <- gsub(".csv","",
                          dir(pattern = "*.csv"),
                          fixed = TRUE)
names(enrich_ldf) <- gsub("test_host_cont_","", names(enrich_ldf), fixed = TRUE)

# add a column with the contrast name
enrich_ldf <- map2(enrich_ldf, names(enrich_ldf), ~ mutate(.x, contrast = .y))

# convert to one dataframe
enrich_one_df <- do.call(rbind, enrich_ldf)


# check for overlap across contrasts without considering timepoint

#filter the enrichment tests by contrast level
all_enrich <- enrich_one_df %>% dplyr::filter(grepl('all', contrast))


list_all_enrich <- list(high = c(as.vector(all_enrich$ID[all_enrich$contrast == "high_all"])), med = c(as.vector(all_enrich$ID[all_enrich$contrast == "med_all"])), low = c(as.vector(all_enrich$ID[all_enrich$contrast == "low_all"])))

# create a venn diagram of overlap
gvenn_all <- ggVennDiagram(list_all_enrich)



# check the overlap of enrichment in contrasts while considering timepoint

# filter the enrichment tests for those that consider timepoint
timepoint_enrich <- enrich_one_df %>% dplyr::filter(grepl('t', contrast))

# pull out just the two columns of interest
timepoint_enrich_select <- timepoint_enrich %>% dplyr::select(ID, contrast)

# turn it into a binary matrix
enrich_timep_wide <- pivot_wider(timepoint_enrich_select, names_from = "contrast", values_from = 'contrast', values_fill = 0, values_fn = function(x) 1)

# create an upset plot from these
upsetR_enrich_t <- as.data.frame(select(enrich_timep_wide, -c(ID)))
upsetR_enrich_t_plot <- UpSetR::upset(upsetR_enrich_t, nsets = 10, order.by = "freq", nintersects = 20) #, query.legend = "top", queries = list(list(query = intersects, params = list('Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "#800026", active = T, query.name="High"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium'), color = "#E31A1C", active = T, query.name="Medium"), list(query = intersects, params = list('Field_vs_Medium', 'Initial_vs_Medium', 'Control_vs_Medium', 'Field_vs_High', 'Initial_vs_High', 'Control_vs_High'), color = "blue", active = T, query.name="Both")))



# filter for the top 30 terms in each contrast
top15_enriched <- timepoint_enrich %>% 
  group_by(contrast) %>% dplyr::slice(1:5)

#pull the GO term for every contrast if it is in the top 10 of ANY contrast
all_values_topGOterms <- timepoint_enrich %>% filter(ID %in% unique(top15_enriched$ID))

# create heatmap
# filter for just the cols of interest and make wide format
top15_enrich_wide <- all_values_topGOterms %>% dplyr::select(ID, Description, Count, contrast) %>% pivot_wider(names_from = "contrast", values_from = 'Count') %>% drop_na(Description)

# remove the ID and description cols
top15_enrich_wide_clean <- top15_enrich_wide %>% dplyr::select(!c(ID,Description))

#reorder the columns
top15_enrich_wide_clean <- top15_enrich_wide_clean[, c("med_t1", "med_t2", "med_t3", "med_t4", "med_t5", "high_t1", "high_t2", "high_t3", "high_t4", "high_t5")]

# set the ID as the rowname
rownames(top15_enrich_wide_clean) <- top15_enrich_wide$Description

# reformat as a matrix
top15_enrich_mat <- as.matrix(top15_enrich_wide_clean)

# create the row and column annotations for the files
# create the column annotation
contrast_order <- as.data.frame(colnames(top15_enrich_mat))
colnames(contrast_order) <- "contrast_id"
contrast_order <- contrast_order %>% tidyr::separate(contrast_id, sep = "_", into = c("Treatment", "Timepoint"), remove = FALSE)
col_annot_contrasts <- contrast_order %>% dplyr::select(Treatment, Timepoint)
rownames(col_annot_contrasts) <- contrast_order$contrast_id

# set the colors for the heatmap
annot_colors_contrasts = list(
  Treatment = c(med = Med_col, high = High_col),
  Timepoint = c(t1 = T1_col, t2 = T2_col, t3 = T3_col, t4 = T4_col, t5 = T5_col))

col_pal <- brewer.pal(9, "YlGn")

# basic heatmap
heat_contrasts_GO <- pheatmap(log(top15_enrich_mat), cluster_cols = F, cluster_rows = F, annotation_col = col_annot_contrasts, annotation_colors = annot_colors_contrasts, show_rownames = T, show_colnames = F, color = col_pal)


# checking concordance between replicates
# first investigate just the top variance events across the whole experiment
#calculate variance for each event for each treatment group
byGroup_var_host <- long_host_no_na_ident %>%
  dplyr::group_by(treatment_timepoint, event_id) %>%
  dplyr::mutate(groupVar = var(psi))

# create a density plot of the variances for each sample type
plot_groupVar_dens <- ggplot(byGroup_var_host, aes(fill = factor(Timepoint, levels = c('T5', 'T4', 'T3', 'T2', 'T1', 'T0', 'TF')), x=groupVar)) + geom_density(alpha = 0.5) + theme_bw() + facet_grid(. ~ Treatment) + scale_fill_manual(values=c("#2166AC", "#92C5DE", "#FDDBC7","#D6604D", "#B2182B", "green", "blue")) + guides(fill = guide_legend(reverse = TRUE)) + labs(fill = "Timepoint")

# filter out the zeros
byGroup_var_no_zeros <- byGroup_var_host %>% dplyr::filter(groupVar > 0)

plot_groupVar_no0_dens <- ggplot(byGroup_var_no_zeros, aes(fill = factor(Timepoint, levels = c('T5', 'T4', 'T3', 'T2', 'T1', 'T0', 'TF')), x=groupVar)) + geom_density(alpha = 0.5) + theme_bw() + facet_grid(Timepoint ~ Treatment) + scale_fill_manual(values=c("#2166AC", "#92C5DE", "#FDDBC7","#D6604D", "#B2182B", "green", "blue")) + guides(fill = guide_legend(reverse = TRUE)) + labs(fill = "Timepoint")

# make the psi binary (spliced in or out) based on > or < 0.5
byGroup_var_host$splice_binary <- ifelse(byGroup_var_host$psi >= 0.5, 1, 0)

# check agreement in binary decision by group
#get the most common decision for each event for each treatment 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

byGroup_var_host <- byGroup_var_host %>% dplyr::group_by(treatment_timepoint, event_id) %>% dplyr::mutate(mode=Mode(splice_binary))

# get the proportion that agree with the mode for each event for each treatment
byGroup_var_host$diff <- byGroup_var_host$splice_binary - byGroup_var_host$mode

byGroup_var_host <- byGroup_var_host %>% dplyr::group_by(treatment_timepoint,event_id) %>% dplyr::mutate(total_reps = n())

# get the proportion that agree
byGroup_var_host <- byGroup_var_host %>% dplyr::group_by(treatment_timepoint,event_id) %>% dplyr::mutate(num_agree = sum(diff == 0))

byGroup_var_host$prop_agree <- byGroup_var_host$num_agree / byGroup_var_host$total_reps

# plot proportion agreement in density plot
plot_prop_agree_dens <- ggplot(byGroup_var_host, aes(fill = factor(Timepoint, levels = c('T5', 'T4', 'T3', 'T2', 'T1', 'T0', 'TF')), x=prop_agree)) + geom_density(alpha = 0.5) + theme_bw() + facet_grid(. ~ Treatment) + scale_fill_manual(values=c("#2166AC", "#92C5DE", "#FDDBC7","#D6604D", "#B2182B", "green", "blue")) + guides(fill = guide_legend(reverse = TRUE)) + labs(fill = "Timepoint")

box_groupVar_host <- ggplot(byGroup_var_no_zeros, aes(x = Timepoint, y = groupVar, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free") 

box_decAgree_host <- ggplot(byGroup_var_host, aes(x = Timepoint, y = prop_agree, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free") 

# count number of events for each group
group_var_summary <- byGroup_var_host %>%
  dplyr::group_by(treatment_timepoint) %>%   
  dplyr::summarise(unique_event_count = n_distinct(event_id))


# get list of genotypes for each sample
group_var_genos <- byGroup_var_no_zeros %>%
  dplyr::group_by(treatment_timepoint) %>%  
  dplyr::summarise(
    unique_event_count = n_distinct(event_id),  
    unique_genotypes = list(unique(Genotype)))

# Downsample the dataframe to just four genotypes per category (so all have the same N) while preserving all columns
downsampled_group_var <- byGroup_var_no_zeros %>%
  dplyr::group_by(treatment_timepoint, event_id) %>%
  group_modify(~ {
    genotypes_in_group <- .x %>%
      distinct(Genotype, .keep_all = TRUE)
    n_genotypes <- nrow(genotypes_in_group)
    sampled_genotypes <- genotypes_in_group %>%
      slice_sample(n = min(4, n_genotypes)) 
    return(sampled_genotypes)
  }) %>%
  ungroup()

# now recalculate the variance based on the new downsampled df
downsampled_group_var <- downsampled_group_var %>%
  dplyr::group_by(treatment_timepoint, event_id) %>%
  dplyr::mutate(groupVar_down = var(psi))

# make a plot of the variance for the downsampled tabled
box_groupVar_downsample <- ggplot(downsampled_group_var, aes(x = Timepoint, y = groupVar_down, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(. ~ Treatment, scales = "free_x", space = "free") 


# read in the annotation files (emapper results run on annot previous annotion to improve GO terms)
annot_host <- read.delim("~/coral_genomes/other/Acer/Acerv.emapper.annotations", sep = "\t", header = T, skip = 4)

annot_host <- annot_host[1:20482,]

# split the transcript ID and gene ID
annot_host_clean_sep <- annot_host %>% tidyr::separate(X.query, sep = "-", into = c("gene_id", "trans_id"), remove = FALSE)

# filter to just one transcript per gene
annot_host_clean <- annot_host_clean_sep %>% dplyr::group_by(gene_id) %>% filter(trans_id == min(trans_id)) %>% distinct()


# annotate the events (combine annot and test file)
host_sig_events_annot <- left_join(split_ID, annot_host_clean, by="gene_id")

#annotate the low var events
low_var_removed_rename <- low_var_removed %>% dplyr::rename(gene_id = gene_name)

low_var_host_annot <- left_join(low_var_removed_rename, annot_host_clean, by="gene_id")


# filter for SF3B1 gene (FUN_010900)
sf3b1_splice_test <- host_sig_events_annot %>% dplyr::filter(gene_id == "FUN_010900")

# filter this for dPSI > 0.3
sf3b1_splice_test_dpsi0.3 <- sf3b1_splice_test %>% dplyr::filter(abs(dPSI) > 0.3)

#create binary matrix
sf3b1_splice_test_dpsi0.3_wide <- sf3b1_splice_test_dpsi0.3 %>% dplyr::select(event_id, vrs) %>% pivot_wider(names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

#create an upset plot
upsetR_sf3b1 <- as.data.frame(select(sf3b1_splice_test_dpsi0.3_wide, -c(event_id)))
upsetR_sf3b1_plot <- UpSetR::upset(upsetR_sf3b1, nsets = 25, keep.order = T, nintersects = 20)



#create a heatmap of SF3b1 
sf3b1_psi_all<- split_ID_host %>% dplyr::filter(gene_name == "FUN_010900")

# filter for list of sig different events in this gene
sf3b1_psi_dpsi0.3_filt <- sf3b1_psi_all %>% filter(event_id %in% sf3b1_splice_test_dpsi0.3$event_id)

#filter for samples with high NA proportion
sf3b1_drop_thresh<- sf3b1_psi_dpsi0.3_filt %>% dplyr::select(which(colMeans(is.na(.)) < 0.5))

#now drop events with high NA
sf3b1_drop_na_events <- sf3b1_drop_thresh[rowMeans(is.na(sf3b1_drop_thresh)) <= .4,]


# filter for just the psi columns
sf3b1_psi_only <- sf3b1_psi_dpsi0.3_filt %>% dplyr::select(!(event_id:gene_name))

# set the event ID as the rowname
rownames(sf3b1_psi_only) <- sf3b1_psi_dpsi0.3_filt$event_id

# reformat as a matrix
mat_sf3b1_psi_only <- as.matrix(sf3b1_psi_only)

# create the row and column annotations for the files
# get the sample order
sf3b1_sample_order_mat <- as.data.frame(colnames(mat_sf3b1_psi_only))
colnames(sf3b1_sample_order_mat) <- "sample_id"
sf3b1_sample_order_mat$sample_base <- gsub("host_|.sorted.psi", "", sf3b1_sample_order_mat$sample_id)

# merge with metadata information
sf3b1_sample_order_meta <- merge(sf3b1_sample_order_mat, metadata, by = "sample_base")

# generate the col annotation
col_annot_sf3b1 <- sf3b1_sample_order_meta %>% dplyr::select(Timepoint, Treatment)
rownames(col_annot_sf3b1) <- sf3b1_sample_order_meta$sample_id

# generate the row annotation
row_annot_sf3b1 <- data.frame(sf3b1_psi_dpsi0.3_filt$event_type)
rownames(row_annot_sf3b1) <- sf3b1_psi_dpsi0.3_filt$event_id
colnames(row_annot_sf3b1) <- c("Event_Type")

# set the colors for the heatmap
annot_colors_sf3b1 = list(
  Treatment = c(Field="blue", Initial="green", Control="#FED976", Low = "#FD8D3C", Medium = "#E31A1C", High = "#800026"),
  Timepoint = c(TF = "blue", T0 = "green", T1 = "#B2182B", T2 = "#D6604D", T3 = "#FDDBC7", T4 = "#92C5DE", T5 = "#2166AC"),
  #Genotype = c("07" = "pink", "31" = "purple", "34" = "dark green", "41" = "maroon", "48" = "black", "50" = "grey", "62" = "yellow", "CM5" = "white"),
  Event_Type  = c(alt_3prime="black", alt_5prime="white", exon_skip = "#999933", intron_retention = "#117733", mult_exon_skip = "#AA4499", mutex_exons = "#888888"))

# basic heatmap
heat_sf3b1<- pheatmap(mat_sf3b1_psi_only, cluster_cols = F, cluster_rows = T, annotation_col = col_annot_sf3b1, annotation_row = row_annot_sf3b1, annotation_colors = annot_colors_sf3b1, show_rownames = T, show_colnames = F)



# filter for HSF1 gene (FUN_002455)
hsf1_splice_test <- host_sig_events_annot %>% dplyr::filter(gene_id == "FUN_002455")

# filter this for dPSI > 0.3
hsf1_splice_test_dpsi0.3 <- hsf1_splice_test %>% dplyr::filter(abs(dPSI) > 0.3)

#create binary matrix
hsf1_splice_test_dpsi0.3_wide <- hsf1_splice_test_dpsi0.3 %>% dplyr::select(event_id, vrs) %>% pivot_wider(names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

#create an upset plot
upsetR_hsf1 <- as.data.frame(select(hsf1_splice_test_dpsi0.3_wide, -c(event_id)))
upsetR_hsf1_plot <- UpSetR::upset(upsetR_hsf1, nsets = 25, keep.order = T, nintersects = 20)


#create a heatmap of FUN_002455 
hsf1_psi_all<- split_ID_host %>% dplyr::filter(gene_name == "FUN_002455")

# filter for list of sig different events in this gene
hsf1_psi_dpsi0.3_filt <- hsf1_psi_all %>% filter(event_id %in% hsf1_splice_test_dpsi0.3$event_id)

#filter for samples with high NA proportion
hsf1_drop_thresh<- hsf1_psi_dpsi0.3_filt %>% dplyr::select(which(colMeans(is.na(.)) < 0.6))

#now drop events with high NA
hsf1_drop_na_events <- hsf1_drop_thresh[rowMeans(is.na(sf3b1_drop_thresh)) <= .3,]

#make it long format
long_hsf1 <- hsf1_drop_na_events %>% tidyr::pivot_longer(cols = c(-event_id, -event_type, -event_num, -gene_name), names_to = "sample_id", values_to = "psi")

#merge w metadata and set the order of variables
long_hsf1$sample_base <- gsub("host_|.sorted.psi", "", long_hsf1$sample_id)
long_hsf1_meta <- merge(long_hsf1, metadata, by = "sample_base")
long_hsf1_meta$Treatment <- factor(long_hsf1_meta$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))
long_hsf1_meta$Timepoint <- factor(long_hsf1_meta$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

# make a boxplot
hsf1_box <- ggplot(long_hsf1_meta, aes(x = Timepoint, y = psi, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "YlOrRd") + facet_grid(event_id~.)

## now make a heatmap
# filter for just the psi columns
hsf1_psi_only <- hsf1_drop_na_events %>% dplyr::select(!(event_id:gene_name))

# set the event ID as the rowname
rownames(hsf1_psi_only) <- hsf1_drop_na_events$event_id

# reformat as a matrix
mat_hsf1_psi_only <- as.matrix(hsf1_psi_only)

# create the row and column annotations for the files
# get the sample order
hsf1_sample_order_mat <- as.data.frame(colnames(mat_hsf1_psi_only))
colnames(hsf1_sample_order_mat) <- "sample_id"
hsf1_sample_order_mat$sample_base <- gsub("host_|.sorted.psi", "", hsf1_sample_order_mat$sample_id)

# merge with metadata information
hsf1_sample_order_meta <- merge(hsf1_sample_order_mat, metadata, by = "sample_base")

# generate the col annotation
col_annot_hsf1 <- hsf1_sample_order_meta %>% dplyr::select(Timepoint, Treatment)
rownames(col_annot_hsf1) <- hsf1_sample_order_meta$sample_id

# generate the row annotation
row_annot_hsf1 <- data.frame(hsf1_psi_dpsi0.3_filt$event_type)
rownames(row_annot_hsf1) <- hsf1_psi_dpsi0.3_filt$event_id
colnames(row_annot_hsf1) <- c("Event_Type")

# set the colors for the heatmap
annot_colors_hsf1 = list(
  Treatment = c(Field=TF_col, Initial=In_col, Control=Cont_col, Low = Low_col, Medium = Med_col, High = High_col),
  Timepoint = c(TF = TF_col, T0 = In_col, T1 = T1_col, T2 = T2_col, T3 = T3_col, T4 = T4_col, T5 = T5_col),
  #Genotype = c("07" = "pink", "31" = "purple", "34" = "dark green", "41" = "maroon", "48" = "black", "50" = "grey", "62" = "yellow", "CM5" = "white"),
  Event_Type  = c(alt_3prime=alt_3_col, alt_5prime=alt_5_col, exon_skip = ES_col, intron_retention = IR_col, mult_exon_skip = multi_ES_col, mutex_exons = mutex_col))

# basic heatmap
heat_hsf1<- pheatmap(mat_hsf1_psi_only, cluster_cols = F, cluster_rows = T, annotation_col = col_annot_hsf1, annotation_row = row_annot_hsf1, annotation_colors = annot_colors_hsf1, show_rownames = T, show_colnames = F, color = viridis(10))



# create a simplified upset plot (exclude Initial and Field contrasts) that includes timepoint
# filter sig events 
hist_sig_events_just_control <- host_sig_events_annot %>% dplyr::filter(!grepl("Field",vrs)) %>% dplyr::filter(!grepl("Initial",vrs))

# filter this for dPSI > 0.3
hist_sig_events_just_control_dpsi0.3 <- hist_sig_events_just_control %>% dplyr::filter(abs(dPSI) > 0.3)

# generate a binary matrix of presence absence for each contrast
sig_control_dpsi0.3_binMat<- hist_sig_events_just_control_dpsi0.3 %>% dplyr::select(event_id, vrs)
wide_sig_control_dpsi0.3_binMat <- pivot_wider(sig_control_dpsi0.3_binMat, names_from = "vrs", values_from = 'vrs', values_fill = 0, values_fn = function(x) 1)

# create a complex upset plot of overlaps of event ids between contrasts
upsetR_events_just_control_dpsi0.3 <- as.data.frame(select(wide_sig_control_dpsi0.3_binMat, -c(event_id)))
list_all_cont <- colnames(upsetR_events_just_control_dpsi0.3)
high_cont <- grep('High', list_all_cont, ignore.case = F, value = T)
med_cont <- grep('Med', list_all_cont, ignore.case = F, value = T)
events_upsetR_plot_control_dpsi0.3 <- UpSetR::upset(upsetR_events_just_control_dpsi0.3, nsets = 56, order.by = "freq", nintersects = 40)


# read in the differentially expressed (not necessarily AS) GO terms to create a you can write out of all
setwd("~/analysis/DGE_DESeq2/GO_enrichment")
GO_ldf <- list() 
GO_listcsv <- dir(pattern = ".tvs") 
for (k in 1:length(GO_listcsv)){
  GO_ldf[[k]] <- read.delim(file = GO_listcsv[k], header = T, sep = "\t")
}

names(GO_ldf) <- gsub(".tvs","",
                      dir(pattern = "*.tvs"),
                      fixed = TRUE)

GO_ldf <- lapply(names(GO_ldf), function(name) {
  GO_ldf <- GO_ldf[[name]]
  GO_ldf$contrast <- name  # Add the dataframe name as a new column
  return(GO_ldf)
})

# create one table you can write out to a file
GO_one_df <- do.call(rbind, GO_ldf)




