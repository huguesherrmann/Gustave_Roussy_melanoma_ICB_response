# ......................................................
# Calculate correlation between pervasive / ALR burdens with all coding genes then run hallmark GSEA
# 29/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")

# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot pervasive and ALR loads with other biomarkers or Bagaev' signatures.")
parser$add_argument("--design", type = "character", help = "Path to design table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

design <- args$design
out_dir <- args$out_dir

counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/gene_counts_gr1234.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
burdens <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/pervasive_ALR_burdens/pervasive_ALR_burden.tsv"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/"
tumor_purity <- TRUE


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
if (mart == "NULL") {
   mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
   #saveRDS(mart, "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds")
} else {
   mart <- readRDS(mart)
}

design <- read_tsv(design, show_col_types = FALSE)

burdens <- read_tsv(burdens, show_col_types = FALSE)

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
   select(gs_name, gene_symbol)
c1 <- msigdbr(species = "Homo sapiens", category = "C1") %>% 
   select(gs_name, gene_symbol)

annotation <- read_tsv(annotation, show_col_types = FALSE)
coding_genes <- annotation %>% filter(gene_biotype == "protein_coding")

all_counts <- read.table(counts, header = TRUE, row.names = 1) %>% select(design$Sample_ID)
if (tumor_purity) {
   all_counts <- normalize_by_tumor_purity(all_counts, design)
}

coding_counts <- all_counts[coding_genes$ensembl_gene_id, ] #%>%
#    t_df() %>%
#    rownames_to_column("Sample_ID")
cc <- coding_counts[, ] %>%
   rownames_to_column("ensembl_gene_id")
x <- get_gene_symbols(cc, mart) %>%
   distinct(hgnc_symbol, .keep_all = TRUE) %>%
   column_to_rownames("hgnc_symbol") %>%
   select(-ensembl_gene_id) %>%
   t_df() %>%
   rownames_to_column("Sample_ID")


# --------------------------------------------
#
#   CALCULATE SPEARMAN CORRELATION WITH CODING GENES
#
# --------------------------------------------
pervasives <- x %>% inner_join(burdens %>% select(Sample_ID, Pervasive_burden), ., by = "Sample_ID") %>%
   column_to_rownames("Sample_ID")

correlation <- cor(pervasives, method = "spearman")
perv_cor <- sort(correlation[-1, 1], decreasing = TRUE) # Remove correlation with the variable itself

x2 <- perform_gsea(perv_cor, c1, pval_cutoff = 0.01)


alr <- x %>% inner_join(burdens %>% select(Sample_ID, ALR_burden), ., by = "Sample_ID") %>%
   column_to_rownames("Sample_ID")

correlation <- cor(alr, method = "spearman")
alr_cor <- sort(correlation[-1, 1], decreasing = TRUE) # Remove correlation with the variable itself

x_alr <- perform_gsea(alr_cor, hallmarks, pval_cutoff = 0.01)
