# ......................................................
# Perform differential expression analysis
# 10/05/23
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_intron_expression/gr1234/"
mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
l2fc <- 0
alpha <- 0.05
tumor_purity <- FALSE
cohort <- "GR"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
mart <- readRDS(mart)

ambiguous_counts <- read.table("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/gr1234/ambiguous_counts/ambiguous_counts.tsv", sep = "\t")
mature_counts <- read.table("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/gr1234/mature_counts/mature_counts.tsv", sep = "\t")
nascent_counts <- read.table("/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/gr1234/nascent_counts/nascent_counts.tsv", sep = "\t")

annotation <- read_tsv(annotation, show_col_types = FALSE)
pcg <- annotation %>% filter(gene_biotype == "protein_coding")
non_coding <- annotation %>% filter(gene_biotype != "protein_coding")

pcg_ambiguous_counts <- ambiguous_counts[pcg$ensembl_gene_id, ]
pcg_mature_counts <- mature_counts[pcg$ensembl_gene_id, ]
pcg_nascent_counts <- nascent_counts[pcg$ensembl_gene_id, ]
# All but PCG
ncg_ambiguous_counts <- ambiguous_counts[setdiff(rownames(ambiguous_counts), pcg$ensembl_gene_id), ]
ncg_mature_counts <- mature_counts[setdiff(rownames(mature_counts), pcg$ensembl_gene_id), ]
ncg_nascent_counts <- nascent_counts[setdiff(rownames(nascent_counts), pcg$ensembl_gene_id), ]

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
   select(gs_name, human_ensembl_gene)

design <- read_tsv(design, show_col_types = FALSE) %>%
   mutate(Batch = factor(Batch)) %>%
   mutate(Tumor_purity = round(Tumor_purity, 2))


# ds <- c("RIGI", "MDA5", "PKR", "ZBP1", "LGP2", "NLRP1", "MAVS", "TBK1", "IRF3", "IRF7", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", 
#         "TLR6", "TLR7", "TLR8", "TLR9", "TLR10", "OAS1", "OAS2", "OAS3", "OASL", "TANK", "TRAIP", "TIFA", "TRAF1", "TRAF2", 
#         "TRAF3", "TRAF4", "TRAF5", "TRAF6", "TRAFD1", "ZFTRAF1")


# ......................................................
#
#   IDENTIFY DIFFERENTIALLY EXPRESSED INTRONS ----
#
# ......................................................
formula <- "~Response"

if ("Batch" %in% colnames(design) & length(unique(design$Batch)) >= 2) {
   design <- design %>% mutate(Batch = as.factor(Batch))
   formula <- paste0(formula, "+Batch")
} else {
   design[, "Batch"] <- rep("1", nrow(design))
   design <- design %>% mutate(Batch = as.factor(Batch))
}

if ("Biopsy_site" %in% colnames(design) & length(unique(design$Biopsy_site)) >= 2) {
   formula <- paste0(formula, "+Biopsy_site")
}

if (tumor_purity) {
   print("Add Tumor_purity as a factor in DESeq2")
   design <- design %>% mutate(Tumor_purity = round(Tumor_purity, 2))
   formula <- paste0(formula, "+Tumor_purity")
}

formula <- as.formula(formula)
contrast <- c("Response", "R", "NR")


# Get size factors from mature matrix
pcg_mature_counts <- pcg_mature_counts[, order(colnames(pcg_mature_counts))]
design <- design[order(design$Sample_ID), ]
dds <- DESeqDataSetFromMatrix(countData = pcg_mature_counts, 
                              colData = design, 
                              design = formula) 
keep <- rowSums(counts(dds)) >= 1000
dds <- dds[keep, ]
mature_size_factors <- estimateSizeFactors(dds)$sizeFactor

# DESeq2 on nascent matrix but with library size of mature matrix
pcg_nascent_counts <- pcg_nascent_counts[, order(colnames(pcg_nascent_counts))]
design <- design[order(design$Sample_ID), ]
dds <- DESeqDataSetFromMatrix(countData = pcg_nascent_counts, 
                              colData = design, 
                              design = formula)
keep <- rowSums(counts(dds)) >= 1000
dds <- dds[keep, ]
dds$sizeFactor <- mature_size_factors
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res_pcg <- results(dds, contrast = contrast, lfcThreshold = l2fc, alpha = alpha)
res_pcg_df <- filter_DEG(res_pcg, alpha, l2fc, mart) %>% mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
   inner_join(., annotation %>% select(ensembl_gene_id, gene_biotype), by = "ensembl_gene_id")
# for volcano plot and ranking (GSEA) only
shrunk <- lfcShrink(dds, coef = "Response_R_vs_NR", type = "ashr", quiet = TRUE)
shrunk_stats_coding <- as.data.frame(shrunk) %>% rownames_to_column("ensembl_gene_id") %>%
   mutate(Expression = if_else(log2FoldChange < 0, "NR", "R")) %>%
   get_gene_symbols(mart)


# ......................................................
#
#   CARACTERIZE COHORT GSEA ----
#
# ......................................................
# Coding
sorted_vector <- get_sorted_vector(shrunk_stats_coding, "log2FoldChange", "ensembl_gene_id")
hallmark_gsea <- perform_gsea(sorted_vector, hallmarks, pval_cutoff = alpha) %>%
   mutate(Condition = factor(if_else(NES > 0, "R", "NR"))) %>% 
   arrange(NES) %>%
   mutate(Cohort = cohort) %>%
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(ID = factor(ID, levels = ID))


# ......................................................
#
#   OVER-REPRESENTATION ANALYSIS OF DE INTRONS ----
#
# ......................................................
ora_nr <- enricher(res_pcg_df %>% filter(Expression == "NR") %>% pull(ensembl_gene_id) , TERM2GENE = hallmarks, pvalueCutoff = alpha) %>% as.data.frame()
ora_r <- enricher(res_pcg_df %>% filter(Expression == "R") %>% pull(ensembl_gene_id) , TERM2GENE = hallmarks, pvalueCutoff = alpha) %>% as.data.frame()
ora <- ora_nr %>% mutate(Condition = "NR") %>%
   add_row(ora_r %>% mutate(Condition = "R")) %>% 
   mutate(Cohort = "GR") %>% 
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(ID = factor(ID, levels = ID))


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(hallmark_gsea, paste0(out_dir, "/hallmark/coding_introns_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora, paste0(out_dir, "/hallmark/ORA_coding_introns_hallmarks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
