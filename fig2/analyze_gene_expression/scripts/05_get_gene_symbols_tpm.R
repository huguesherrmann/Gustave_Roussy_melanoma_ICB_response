# ......................................................
# Get gene symbols
# 05/03/22
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/functions_for_gene_expression_analysis.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Get gene symbols.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort. For output name formating.")
parser$add_argument("--tpm", type = "character", help = "Path to tpm table.")
parser$add_argument("--mart", type = "character", default = "NULL", help = "Path to ensembl mart object.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

cohort <- args$cohort
tpm <- args$tpm
design <- args$design
mart <- args$mart
out_dir <- args$out_dir
# cohort <- "gr123"
# tpm <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr123/scaled_tpm_gene_counts_gr123.tsv"
# mart <- "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
tpm <- read.table(tpm, row.names = 1, header = TRUE) %>%
  rownames_to_column("ensembl_gene_id")

if (mart == "NULL") {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  #saveRDS(mart, "/mnt/beegfs/userdata/h_herrmann/tools_software/R_mart_hsapiens_gene_ensembl/R_mart_hsapiens_gene_ensembl.rds")
} else {
  mart <- readRDS(mart)
}


# --------------------------------------------
#
#   GET GENE SYMBOLS
#
# --------------------------------------------
symbol_tpm <- get_gene_symbols(tpm, mart) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  relocate(hgnc_symbol, .before = 1) %>%
  select(-ensembl_gene_id)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
write.table(symbol_tpm, paste0(out_dir, "/Bagaev/symbols_tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
