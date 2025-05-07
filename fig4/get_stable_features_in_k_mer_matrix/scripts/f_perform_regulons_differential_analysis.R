# ......................................................
# Repeated cross-validation differential expression analysis on regulons
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
parser <- ArgumentParser(description = "Repeated cross-validation differential expression analysis on regulons.")
parser$add_argument("--regulons", type = "character", help = "Path to regulon counts table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--indice_dir", type = "character", help = "Path to directory containing all train and test partitions.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--l2fc", type = "double", default = 0, help = "Log2 fold change threshold for differential expression test.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for differential expression test.")
parser$add_argument("--threads", type = "integer", default = 1, help = "Number of cores to be used.")
parser$add_argument("--tumor_purity", type = "character", default = "FALSE", help = "Should signatures be computed with tumor purity normalization?.")
args <- parser$parse_args()

regulons <- args$regulons
correspondence <- args$correspondence
design <- args$design
indice_dir <- args$indice_dir
out_dir <- args$out_dir
l2fc <- args$l2fc
alpha <- args$alpha
threads <- args$threads
tumor_purity <- args$tumor_purity
# regulons <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/regulons/all_regulon_counts.tsv"
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/regulons/all_correspondence_contigs_regulons.tsv"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# indice_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_cross_validation_partitioning/gr1234_20_5/"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_cv_20_5/stability/"
# threads <- 1
# tumor_purity <- FALSE


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE) %>%
  mutate(Batch = as.factor(Batch)) %>%
  mutate(Biopsy_site = as.factor(Biopsy_site)) %>%
  mutate(Response = as.factor(Response))

regulons <- fread(regulons)

correspondence <- fread(correspondence)

index_df <- load_indices(indice_dir)


# --------------------------------------------
#
#   FILTER OUT INTERGENIC REGIONS WITH TOO FEW CONTIGS
#
# --------------------------------------------
cluster_regulons <- regulons[!grepl("*_0_*", Regulon), ]
cluster_correspondence <- correspondence[!grepl("*_0_*", Regulon), ]
# Keep intergenic regulons with at least 20 contigs
too_few_intergenic_regulons <- cluster_correspondence[grepl("chr", Regulon), .N, by = .(Regulon)] %>%
  .[N < 20, ]
cluster_regulons <- cluster_regulons[Regulon %ni% too_few_intergenic_regulons$Regulon, ] %>%
  column_to_rownames("Regulon")

cluster_regulons <- cluster_regulons %>% round()

  
rm(regulons)
rm(correspondence)
rm(cluster_correspondence)
gc()


# --------------------------------------------
#
#   DIFFERENTIAL ANALYSIS
#
# --------------------------------------------
formula <- "~ Response + Biopsy_site + Batch"
if (tumor_purity) {
  cat("Add Tumor_purity as a factor in DESeq2\n")
  design <- design %>% mutate(Tumor_purity = round(Tumor_purity, 2))
  formula <- paste0(formula, "+Tumor_purity")
}
formula <- as.formula(formula)
contrast <- c("Response", "R", "NR")

cat(paste0("Run differential analysis on ", threads, " cores...\n"))
register(MulticoreParam(threads))

n_partitions <- max(index_df$Id)
results <- foreach(n = 1:n_partitions, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %do% {
  train_run <- index_df %>% filter(Id == n, Split == "train") %>%
    select(-Response) %>%
    inner_join(., design, by = "Sample_ID")
  train_feature <- cluster_regulons %>% select(train_run$Sample_ID)
  
  deseq <- identify_diff_expr_genes(train_feature, train_run, formula, contrast, l2fc, alpha, parallel = TRUE)
  dds <- deseq$dds_object
  deseq_contrast <- deseq$contrast
  deg <- deseq_contrast %>% as.data.frame() %>%
    filter(padj < alpha & abs(log2FoldChange) > l2fc) %>%
    select(log2FoldChange, padj) %>%
    rownames_to_column("Regulon") %>%
    mutate(Id = n)

  
  deg
}

results <- results %>% mutate(Condition = if_else(log2FoldChange > 0, "responder", "non_responder"))
cat("Differential analysis done.")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
write.table(results, paste0(out_dir, "stability/differential_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
