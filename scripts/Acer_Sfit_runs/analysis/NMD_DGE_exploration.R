# script for NMD pathway investigation
author: "Kate Stankiewicz"
date: "2024-10-31"
output: html_document

# Load libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(rstatix)
library(apeglm)
library(data.table)
library(DESeq2)
library(pheatmap)
library(viridis)


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


# read in the counts files 
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


#extract the normalized counts
counts_norm <- DESeq2::counts(dds_host, normalized = TRUE)
counts_norm <- as.data.frame(counts_norm)

# convert rownames into gene id column and split from strip transcript
counts_norm$gene_name <- row.names(counts_norm)

#convert the matrix to long format
counts_norm_long <- counts_norm %>% tidyr::pivot_longer(cols = -c("gene_name"), names_to = "sampleName", values_to = "normalized_count")


# read in the file of all DGE
DGE <- read.delim("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sumGene_all_sig_DGE_lfc_sval.csv", header = T, sep = ",")

# annot files ran emapper on annot from Nick to get GO terms
# read in the annotation files
annot_host <- read.delim("~/coral_genomes/other/Acer/Acerv.emapper.annotations", sep = "\t", header = T, skip = 4)

annot_host <- annot_host[1:20482,]

# split the transcript ID and gene ID
annot_host_clean_sep <- annot_host %>% tidyr::separate(X.query, sep = "-", into = c("gene_id", "trans_id"), remove = FALSE)

# filter to just one transcript per gene
annot_host_clean <- annot_host_clean_sep %>% dplyr::group_by(gene_id) %>% filter(trans_id == min(trans_id)) %>% distinct()

# annotate the significant events (combine annot and test file)
#AS_DGE_overlap_annot <- left_join(AS_DGE, annot_host_clean, by="gene_id")
DGE_annot <- left_join(DGE, annot_host_clean, by="gene_id")
DGE_up <- DGE_annot %>% filter(log2FoldChange > 0)
DGE_down <- DGE_annot %>% filter(log2FoldChange < 0)


# pull NMD related genes and check expression pattern

# pull out NMD genes
NMD_genes <- dplyr::filter(DGE_annot, grepl("nonsense",Description, ignore.case = T))

# filter out the field and initial contrasts
NMD_genes_controls <- NMD_genes %>% dplyr::filter(str_detect(vrs, "Field|Initial", negate = TRUE))


# select just the columns of interest
NMD_forheat <- NMD_genes_controls %>% dplyr::select(gene_id, log2FoldChange, vrs)

# convert to wide format
NMD_wide <- spread(NMD_forheat, vrs, log2FoldChange)

# set the gene ID as the rowname
NMD_wide<- NMD_wide %>% column_to_rownames(var = "gene_id") %>% dplyr::select(Low_T2_vs_Control_T2, Low_T3_vs_Control_T3, Medium_T1_vs_Control_T1, Medium_T2_vs_Control_T2, Medium_T3_vs_Control_T3, Medium_T4_vs_Control_T4, Medium_T5_vs_Control_T5, High_T1_vs_Control_T1, High_T2_vs_Control_T2, High_T3_vs_Control_T3, High_T4_vs_Control_T4, High_T5_vs_Control_T5)

# reformat as a matrix
NMD_mat <- as.matrix(NMD_wide)

#NMD_mat[is.na(NMD_mat)] <- 0

# set the breaks
bk = c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)
mycols = colorRampPalette(colors = c("blue","white","red"))(15)

# generate a heatmap
heat_NMD <- pheatmap(NMD_mat, breaks = bk, color = mycols, cluster_cols = F, cluster_rows = F, show_colnames = T)

# add metadata to normalized counts
counts_norm_meta <- merge(counts_norm_long, metadata, by = "sampleName")

# filter normalized counts by NMD gene IDs
counts_NMD <- counts_norm_meta %>% dplyr::filter(gene_name %in% NMD_genes_controls$gene_id)


# make a boxplot

# set colors for treatments
#Treatments
Cont_col="#FED976"
Low_col = "#FD8D3C"
Med_col = "#E31A1C"
High_col = "#800026"

# adjust the order
counts_NMD$Treatment <- factor(counts_NMD$Treatment, levels = c("Field", "Initial", "Control", "Low", "Medium", "High"))
counts_NMD$Timepoint <- factor(counts_NMD$Timepoint, levels = c("TF", "T0", "T1", "T2", "T3", "T4", "T5"))

# plot the boxplot
box_NMD_counts <- ggplot(counts_NMD, aes(x = Timepoint, y = normalized_count, fill = Treatment)) + geom_boxplot() + theme_bw() + scale_fill_manual(values = c(TF_col, In_col, Cont_col, Low_col, Med_col, High_col)) + facet_grid(gene_name ~ ., scales = "free")


# pull the annotations for the NMD genes
NMD_gene_IDs <- as.vector(unique(NMD_genes_controls$gene_id))

NMD_gene_annots <- annot_host_clean_sep %>% dplyr::filter(gene_id %in% NMD_gene_IDs)

NMD_gene_annots_unique <- NMD_gene_annots %>%
  distinct(gene_id, .keep_all = TRUE)  # Keeps only the first occurrence of each gene_id

# Create new column with the preferred name replacing "FUN"
NMD_gene_annots_unique <- NMD_gene_annots_unique %>%
  mutate(new_gene_id = paste(Preferred_name, sub("^FUN_", "", gene_id), sep = "_"))




