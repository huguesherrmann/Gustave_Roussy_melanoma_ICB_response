# ......................................................
# Repeated cross-validation differential expression analysis on non-coding genes
# 13/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Repeated cross-validation differential expression analysis on regulons.")
# parser$add_argument("--regulons", type = "character", help = "Path to regulon counts table.")
# parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
# parser$add_argument("--design", type = "character", help = "Path to design table.")
# parser$add_argument("--indice_dir", type = "character", help = "Path to directory containing all train and test partitions.")
# parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
# parser$add_argument("--l2fc", type = "double", default = 0, help = "Log2 fold change threshold for differential expression test.")
# parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for differential expression test.")
# parser$add_argument("--threads", type = "integer", default = 1, help = "Number of cores to be used.")
# args <- parser$parse_args()
# 
# regulons <- args$regulons
# correspondence <- args$correspondence
# design <- args$design
# indice_dir <- args$indice_dir
# out_dir <- args$out_dir
# l2fc <- args$l2fc
# alpha <- args$alpha
# threads <- args$threads
counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/gene_counts_gr1234.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
indice_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_cross_validation_partitioning/gr1234_20_5/"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/"
l2fc <- 0
alpha <- 0.05
threads <- 4


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE) %>%
   mutate(Batch = as.factor(Batch)) %>%
   mutate(Biopsy_site = as.factor(Biopsy_site)) %>%
   mutate(Response = as.factor(Response))

annotation <- read_tsv(annotation, show_col_types = FALSE)
non_coding_genes <- annotation %>% filter(gene_biotype != "protein_coding")
non_coding_counts <- read.table(counts, header = TRUE, row.names = 1)[non_coding_genes$ensembl_gene_id, ] %>% select(design$Sample_ID)

index_df <- load_indices(indice_dir)


# --------------------------------------------
#
#   DIFFERENTIAL ANALYSIS
#
# --------------------------------------------
formula <- as.formula("~ Response + Biopsy_site + Batch")
contrast <- c("Response", "R", "NR")

cat(paste0("Run differential analysis on ", threads, " cores...\n"))
register(MulticoreParam(threads))

n_partitions <- 50#max(index_df$Id)
results <- foreach(n = 1:n_partitions, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %do% {
   train_run <- index_df %>% filter(Id == n, Split == "train") %>%
      select(-Response) %>%
      inner_join(., design, by = "Sample_ID")
   train_feature <- non_coding_counts %>% select(train_run$Sample_ID)
   
   deseq <- identify_diff_expr_genes(train_feature, train_run, formula, contrast, l2fc, alpha, parallel = TRUE)
   dds <- deseq$dds_object
   deseq_contrast <- deseq$contrast
   deg <- deseq_contrast %>% as.data.frame() %>%
      filter(padj < alpha & abs(log2FoldChange) > l2fc) %>%
      select(log2FoldChange, padj) %>%
      rownames_to_column("Ensembl_gene_id") %>%
      mutate(Id = n) %>%
      inner_join(., annotation, by = c("Ensembl_gene_id" = "ensembl_gene_id"))
   
   
   deg
}

results <- results %>% mutate(Condition = if_else(log2FoldChange > 0, "R", "NR"))
cat("Differential analysis done.")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
write.table(results, paste0(out_dir, "stability/differential_non_coding.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
