# functional enrichment DGE (and overlap with AS) Acer


# Load libraries
require(clusterProfiler)
require(AnnotationHub)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(viridis)
library(tidyverse)


# data pre-processing

# read in the file of AS and DGE overlap
AS_DGE <- read.delim("~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/sumGene_AS_DE_contrasts_overlap_lfc_sval.csv", header = T, sep = ",")

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


## Functional enrichment

# prepare the full annotation file for clusterProfiler background genes
# host
# pull relevant columns (make sure GO terms are in col 1 for clusterprofiler)
host_all_go <- annot_host_clean[, c("GOs", "gene_id")]

# put each go term in its own row
host_go_all_long <- host_all_go %>% ungroup() %>% separate_rows(GOs, sep = ",")

#remove rows with no go terms
host_go_all_noNA <- host_go_all_long[!(host_go_all_long$GOs == "" | is.na(host_go_all_long$GOs)), ]
host_go_all_noNA <- host_go_all_noNA[!(host_go_all_noNA$GOs == "-" | is.na(host_go_all_noNA$GOs)), ]
host_go_all_noNA <- na.omit(host_go_all_noNA)

#trim whitespace
host_go_all_clean <- host_go_all_noNA %>% mutate_if(is.character, str_trim)


# process the GO terms for the significantly diff spliced genes in the same way for each contrast
# host
AS_DGE_go_up <- DGE_up[, c("GOs", "gene_id", "vrs")]
AS_DGE_go_down <- DGE_down[, c("GOs", "gene_id", "vrs")]

#split the contrast id into parts
AS_DGE_split_up <- AS_DGE_go_up %>% tidyr::separate(vrs, sep = "_", into = c("var1", "timepoint1", "vs", "var2", "timepoint2"), remove = F)
AS_DGE_split_down <- AS_DGE_go_down %>% tidyr::separate(vrs, sep = "_", into = c("var1", "timepoint1", "vs", "var2", "timepoint2"), remove = F)

# concatenate the var types
AS_DGE_split_up$treatment_contrast <- paste0(AS_DGE_split_up$var1, "_vs_", AS_DGE_split_up$var2)
AS_DGE_split_down$treatment_contrast <- paste0(AS_DGE_split_down$var1, "_vs_", AS_DGE_split_down$var2)

# put each go term in its own row
host_go_vrs_long_up <- AS_DGE_split_up %>% separate_rows(GOs, sep = ",")
host_go_vrs_long_down <- AS_DGE_split_down %>% separate_rows(GOs, sep = ",")

#remove rows with no go terms
host_go_vrs_no_NA_up <- host_go_vrs_long_up[!(host_go_vrs_long_up$GOs == "" | is.na(host_go_vrs_long_up$GOs)), ]
host_go_vrs_no_NA_up <- host_go_vrs_no_NA_up[!(host_go_vrs_no_NA_up$GOs == "-" | is.na(host_go_vrs_no_NA_up$GOs)), ]
host_go_vrs_no_NA_up <- na.omit(host_go_vrs_no_NA_up)

host_go_vrs_no_NA_down <- host_go_vrs_long_down[!(host_go_vrs_long_down$GOs == "" | is.na(host_go_vrs_long_down$GOs)), ]
host_go_vrs_no_NA_down <- host_go_vrs_no_NA_down[!(host_go_vrs_no_NA_down$GOs == "-" | is.na(host_go_vrs_no_NA_down$GOs)), ]
host_go_vrs_no_NA_down <- na.omit(host_go_vrs_no_NA_down)

#trim whitespace
host_go_vrs_clean_up <- host_go_vrs_no_NA_up %>% mutate_if(is.character, str_trim)
host_go_vrs_clean_down <- host_go_vrs_no_NA_down %>% mutate_if(is.character, str_trim)

# split by contrast for testing
host_go_low_up <- host_go_vrs_clean_up %>% filter(treatment_contrast == "Low_vs_Control")
host_go_med_up <- host_go_vrs_clean_up %>% filter(treatment_contrast == "Medium_vs_Control")
host_go_high_up <- host_go_vrs_clean_up %>% filter(treatment_contrast == "High_vs_Control")

host_go_low_down <- host_go_vrs_clean_down %>% filter(treatment_contrast == "Low_vs_Control")
host_go_med_down <- host_go_vrs_clean_down %>% filter(treatment_contrast == "Medium_vs_Control")
host_go_high_down <- host_go_vrs_clean_down %>% filter(treatment_contrast == "High_vs_Control")

# get just the columns of interest
host_go_low_up <- host_go_low_up[, c("GOs", "gene_id")]
host_go_med_up <- host_go_med_up[, c("GOs", "gene_id")]
host_go_high_up <- host_go_high_up[, c("GOs", "gene_id")]

host_go_low_down <- host_go_low_down[, c("GOs", "gene_id")]
host_go_med_down <- host_go_med_down[, c("GOs", "gene_id")]
host_go_high_down <- host_go_high_down[, c("GOs", "gene_id")]



# make everything distinct
distinct_host_all <- host_go_all_clean %>% distinct()

distinct_host_go_low_up <- host_go_low_up %>% distinct()
distinct_host_go_med_up <- host_go_med_up %>% distinct()
distinct_host_go_high_up <- host_go_high_up %>% distinct()

distinct_host_go_low_down <- host_go_low_down %>% distinct()
distinct_host_go_med_down <- host_go_med_down %>% distinct()
distinct_host_go_high_down <- host_go_high_down %>% distinct()


# build go maps for host
#rename columns 
colnames(distinct_host_all) <- c("Goterms", "GeneID")

#build the go map (this will take a long time to run)
host_gomap <- clusterProfiler::buildGOmap(as.data.frame(distinct_host_all))




# now lets see what is enriched for host low

# build gene and universe arguments
genes_host_low_up <- as.vector(distinct_host_go_low_up$gene_id)
genes_host_low_down <- as.vector(distinct_host_go_low_down$gene_id)
all_genes_host <- as.vector(distinct_host_all$GeneID)

#prepare TERM2NAMES argument
names_host <- clusterProfiler::go2term(as.vector(host_gomap$GO))

# run test
host_test_low_up = clusterProfiler::enricher(gene=genes_host_low_up, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_low_down = clusterProfiler::enricher(gene=genes_host_low_down, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)



# now lets see what is enriched for host med

# build gene argument
genes_host_med_up <- as.vector(distinct_host_go_med_up$gene_id)
genes_host_med_down <- as.vector(distinct_host_go_med_down$gene_id)

# run test
host_test_med_up = clusterProfiler::enricher(gene=genes_host_med_up, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_med_down = clusterProfiler::enricher(gene=genes_host_med_down, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)



# now lets see what is enriched for host high

# build gene argument
genes_host_high_up <- as.vector(distinct_host_go_high_up$gene_id)
genes_host_high_down <- as.vector(distinct_host_go_high_down$gene_id)

# run test
host_test_high_up = clusterProfiler::enricher(gene=genes_host_high_up, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_high_down = clusterProfiler::enricher(gene=genes_host_high_down, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)



# plot the results
dotplot(host_test_low_up, showCategory = 30)
dotplot(host_test_med_up, showCategory = 30)
dotplot(host_test_high_up, showCategory = 30)

dotplot(host_test_low_down, showCategory = 30)
dotplot(host_test_med_down, showCategory = 30)
dotplot(host_test_high_down, showCategory = 30)


# run tests for each timepoint for high and medium categories


# separate out contrasts by timepoint for medium and high
host_go_med_T1 <- host_go_vrs_clean %>% filter(cont == "Medium_T1_vs_Control_T1")
host_go_med_T2 <- host_go_vrs_clean %>% filter(cont == "Medium_T2_vs_Control_T2")
host_go_med_T3 <- host_go_vrs_clean %>% filter(cont == "Medium_T3_vs_Control_T3")
host_go_med_T4 <- host_go_vrs_clean %>% filter(cont == "Medium_T4_vs_Control_T4")
host_go_med_T5 <- host_go_vrs_clean %>% filter(cont == "Medium_T5_vs_Control_T5")

host_go_high_T1 <- host_go_vrs_clean %>% filter(cont == "High_T1_vs_Control_T1")
host_go_high_T2 <- host_go_vrs_clean %>% filter(cont == "High_T2_vs_Control_T2")
host_go_high_T3 <- host_go_vrs_clean %>% filter(cont == "High_T3_vs_Control_T3")
host_go_high_T4 <- host_go_vrs_clean %>% filter(cont == "High_T4_vs_Control_T4")
host_go_high_T5 <- host_go_vrs_clean %>% filter(cont == "High_T5_vs_Control_T5")

# get just the columns of interest
host_go_med_T1 <- host_go_med_T1[, c("GOs", "gene_id")]
host_go_med_T2 <- host_go_med_T2[, c("GOs", "gene_id")]
host_go_med_T3 <- host_go_med_T3[, c("GOs", "gene_id")]
host_go_med_T4 <- host_go_med_T4[, c("GOs", "gene_id")]
host_go_med_T5 <- host_go_med_T5[, c("GOs", "gene_id")]

host_go_high_T1 <- host_go_high_T1[, c("GOs", "gene_id")]
host_go_high_T2 <- host_go_high_T2[, c("GOs", "gene_id")]
host_go_high_T3 <- host_go_high_T3[, c("GOs", "gene_id")]
host_go_high_T4 <- host_go_high_T4[, c("GOs", "gene_id")]
host_go_high_T5 <- host_go_high_T5[, c("GOs", "gene_id")]

# make everything distinct
distinct_host_go_med_T1 <- host_go_med_T1 %>% distinct()
distinct_host_go_med_T2 <- host_go_med_T2 %>% distinct()
distinct_host_go_med_T3 <- host_go_med_T3 %>% distinct()
distinct_host_go_med_T4 <- host_go_med_T4 %>% distinct()
distinct_host_go_med_T5 <- host_go_med_T5 %>% distinct()

distinct_host_go_high_T1 <- host_go_high_T1 %>% distinct()
distinct_host_go_high_T2 <- host_go_high_T2 %>% distinct()
distinct_host_go_high_T3 <- host_go_high_T3 %>% distinct()
distinct_host_go_high_T4 <- host_go_high_T4 %>% distinct()
distinct_host_go_high_T5 <- host_go_high_T5 %>% distinct()



# run the tests

# build gene arguments
genes_host_med_T1 <- as.vector(distinct_host_go_med_T1$gene_id)
genes_host_med_T2 <- as.vector(distinct_host_go_med_T2$gene_id)
genes_host_med_T3 <- as.vector(distinct_host_go_med_T3$gene_id)
genes_host_med_T4 <- as.vector(distinct_host_go_med_T4$gene_id)
genes_host_med_T5 <- as.vector(distinct_host_go_med_T5$gene_id)

genes_host_high_T1 <- as.vector(distinct_host_go_high_T1$gene_id)
genes_host_high_T2 <- as.vector(distinct_host_go_high_T2$gene_id)
genes_host_high_T3 <- as.vector(distinct_host_go_high_T3$gene_id)
genes_host_high_T4 <- as.vector(distinct_host_go_high_T4$gene_id)
genes_host_high_T5 <- as.vector(distinct_host_go_high_T5$gene_id)


# run tests
# medium
host_test_med_t1 = clusterProfiler::enricher(gene=genes_host_med_T1, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_med_t2 = clusterProfiler::enricher(gene=genes_host_med_T2, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_med_t3 = clusterProfiler::enricher(gene=genes_host_med_T3, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_med_t4 = clusterProfiler::enricher(gene=genes_host_med_T4, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_med_t5 = clusterProfiler::enricher(gene=genes_host_med_T5, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)


# high
host_test_high_t1 = clusterProfiler::enricher(gene=genes_host_high_T1, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_high_t2 = clusterProfiler::enricher(gene=genes_host_high_T2, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_high_t3 = clusterProfiler::enricher(gene=genes_host_high_T3, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_high_t4 = clusterProfiler::enricher(gene=genes_host_high_T4, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)

host_test_high_t5 = clusterProfiler::enricher(gene=genes_host_high_T5, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=all_genes_host, minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.2, TERM2GENE=host_gomap, TERM2NAME = names_host)



# make plots
# med
dotplot(host_test_med_t1, showCategory = 30)
dotplot(host_test_med_t2, showCategory = 30)
dotplot(host_test_med_t3, showCategory = 30)
dotplot(host_test_med_t4, showCategory = 30)
dotplot(host_test_med_t5, showCategory = 30)

#high
dotplot(host_test_high_t1, showCategory = 30)
dotplot(host_test_high_t2, showCategory = 30)
dotplot(host_test_high_t3, showCategory = 30)
dotplot(host_test_high_t4, showCategory = 30)
dotplot(host_test_high_t5, showCategory = 30)




# save the tests to files

# convert to data frames
sum_test_low_up <- data.frame(host_test_low_up)
sum_test_med_up <- data.frame(host_test_med_up)
sum_test_high_up <- data.frame(host_test_high_up)

sum_test_low_down <- data.frame(host_test_low_down)
sum_test_med_down <- data.frame(host_test_med_down)
sum_test_high_down <- data.frame(host_test_high_down)


# write out the tables
write.table(sum_test_low_up, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_low_up.tvs", quote = F, row.names = F, sep = "\t")

write.table(sum_test_med_up, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_med_up.tvs", quote = F, row.names = F, sep = "\t")

write.table(sum_test_high_up, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_high_up.tvs", quote = F, row.names = F, sep = "\t")

write.table(sum_test_low_down, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_low_down.tvs", quote = F, row.names = F, sep = "\t")

write.table(sum_test_med_down, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_med_down.tvs", quote = F, row.names = F, sep = "\t")

write.table(sum_test_high_down, "~/alt_splice/GS_Pilot_Acer/results/Acer_Sfit_runs/analysis/DGE_DESeq2/GO_enrichment/DGE_high_down.tvs", quote = F, row.names = F, sep = "\t")




# check which go terms are common
common_ids_up <- Reduce(intersect, list(sum_test_low_up$ID, sum_test_med_up$ID, sum_test_high_up$ID))
common_ids_down <- Reduce(intersect, list(sum_test_low_down$ID, sum_test_med_down$ID, sum_test_high_down$ID))

# now subset one of the dfs for each to look at these common term descriptions
subset_up <- subset(sum_test_low_up, ID %in% common_ids_up)
subset_down <- subset(sum_test_low_down, ID %in% common_ids_down)


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

# generate a heatmap
heat_NMD <- pheatmap(NMD_mat, cluster_cols = F, cluster_rows = F, show_colnames = T, viridis(10))



