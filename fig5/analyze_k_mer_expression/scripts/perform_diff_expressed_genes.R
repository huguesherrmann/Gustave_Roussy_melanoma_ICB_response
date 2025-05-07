# ......................................................
# Perform differential expression analysis
# 07/10/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Get cutoff to classify samples as pervasives or ALR positives.")
# parser$add_argument("--design", type = "character", help = "Path to design table.")
# args <- parser$parse_args()
# 
# regulon_counts <- args$regulon_counts

counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/gene_counts_gr1234.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/pervasive_ALR_burdens/pervasive_status_design_gr1234.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/"
mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
tumor_purity <- TRUE
l2fc <- 0
alpha <- 0.05


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
   select(gs_name, gene_symbol)

if (mart == "NULL") {
   mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
   #saveRDS(mart, "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds")
} else {
   mart <- readRDS(mart)
}

annotation <- read_tsv(annotation, show_col_types = FALSE)
coding_genes <- annotation %>% filter(gene_biotype == "protein_coding")
non_coding_genes <- annotation %>% filter(ensembl_gene_id %in% setdiff(annotation$ensembl_gene_id, coding_genes$ensembl_gene_id))

design <- read_tsv(design, show_col_types = FALSE) %>%
   mutate(ALR_status = if_else(ALR_status == "1", "alr", "wt")) %>%
   mutate(Batch = as.factor(Batch))

all_counts <- read.table(counts, header = TRUE, row.names = 1) %>% select(design$Sample_ID)
if (tumor_purity) {
   all_counts <- normalize_by_tumor_purity(all_counts, design)
}

coding_counts <- all_counts[coding_genes$ensembl_gene_id, ]
non_coding_counts <- all_counts[non_coding_genes$ensembl_gene_id, ]


# ......................................................
#
#   IDENTIFY DIFFERENTIALLY EXPRESSED GENES ----
#
# ......................................................
formula <- as.formula("~ALR_status + Batch + Biopsy_site")
contrast <- c("ALR_status", "alr", "wt")

coding_deseq <- identify_diff_expr_genes(coding_counts, design, formula, contrast, l2fc, alpha)
coding_dds <- coding_deseq$dds_object
coding_deseq_contrast <- coding_deseq$contrast
coding_deg <- filter_DEG(coding_deseq_contrast, alpha, l2fc, mart) %>% mutate(Expression = if_else(log2FoldChange < 0, "wt", "alr"))

non_coding_deseq <- identify_diff_expr_genes(non_coding_counts, design, formula, contrast, l2fc, alpha)
non_coding_dds <- non_coding_deseq$dds_object
non_coding_deseq_contrast <- non_coding_deseq$contrast
non_coding_deg <- filter_DEG(non_coding_deseq_contrast, alpha, l2fc, mart) %>% mutate(Expression = if_else(log2FoldChange < 0, "wt", "alr"))

n_deg_df <- coding_deg %>% group_by(Expression) %>% 
   summarise(N = n()) %>%
   mutate(Analysis = "Coding") %>%
   rbind(non_coding_deg %>% group_by(Expression) %>% 
            summarise(N = n()) %>%
            mutate(Analysis = "Non_coding"))


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(coding_deg, paste0(out_dir, "DE/coding_deg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
coding_deg <- read_tsv(paste0(out_dir, "DE/coding_deg.tsv"), show_col_types = FALSE)


# ......................................................
#
#   IDENTIFY SIGNIFICANT HALLMARKS ----
#
# ......................................................
sorted_vector <- get_sorted_vector(coding_deg, "log2FoldChange", "hgnc_symbol")
hallmark_gsea <- perform_gsea(sorted_vector, hallmarks, pval_cutoff = alpha) %>%
   mutate(Condition = factor(if_else(NES > 0, "wt", "alr"))) %>% 
   mutate(Cohort = "GR1234") %>% 
   arrange(NES) %>%
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(ID = factor(ID, levels = ID))

xx <- head(hallmarks, 10)
hallmark_gsea <- perform_gsea(sorted_vector, xx, pval_cutoff = 0.001, min_set_size = 1)


entrez_genes <- c("RIGI", "MDA5", "PKR", "ZBP1", "LGP2", "NLRP1", "MAVS", "TBK1", "IRF3", "IRF7",
                  "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TLR10",
                  "OAS1", "OAS2", "OAS3", "OASL", "TANK", "TRAIP", "TIFA", "TRAF1", "TRAF2", "TRAF3", 
                  "TRAF4", "TRAF5", "TRAF6", "TRAFD1", "ZFTRAF1")
gomap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"), filters = "hgnc_symbol", values = entrez_genes, mart = mart)
xx_sig <- tibble(gs_name = rep("Sig", 33),
                 gene_symbol = gomap$ensembl_gene_id)
xx_sig <- tibble(gs_name = rep("Sig", 33),
                 gene_symbol = head(hallmarks$gene_symbol, 33))
sorted_vector <- get_sorted_vector(coding_deg, "log2FoldChange", "ensembl_gene_id")
hallmark_gsea2 <- perform_gsea(sorted_vector, xx_sig, pval_cutoff = alpha) #%>%
   mutate(Condition = factor(if_else(NES > 0, "wt", "alr"))) %>% 
   mutate(Cohort = "GR1234") %>% 
   arrange(NES) %>%
   mutate(ID = str_remove(ID, "HALLMARK_")) %>%
   mutate(ID = factor(ID, levels = ID))
