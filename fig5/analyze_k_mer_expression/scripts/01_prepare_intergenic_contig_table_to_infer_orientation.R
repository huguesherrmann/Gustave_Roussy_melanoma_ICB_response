# ......................................................
# Prepare an intergenic contig table to infer orientation of intergenic regions
# 20/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Prepare an intergenic contig table to infer orientation of intergenic regionss.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulon table.")
parser$add_argument("--annotation", type = "character", help = "Path to regulon annotation table.")
parser$add_argument("--out_dir", type = "character", help = "Path to output results.")
args <- parser$parse_args()

differential <- args$differential
annotation <- args$annotation
out_dir <- args$out_dir
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/stability/differential_regulons.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/regulons/all_regulons_annotated.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA ----
#
# --------------------------------------------
annotation <- fread(annotation, select = c("tag", "contig", "Regulon"))

differential <- read_tsv(differential, show_col_types = FALSE)


# --------------------------------------------
#
#   PREPARE TABLE ----
#
# --------------------------------------------
intergenics <- merge(differential %>% filter(grepl("chr", Regulon)), annotation, by = "Regulon") %>%
   select(tag, contig, Regulon, Condition)


# --------------------------------------------
#
#   EXPORT ----
#
# --------------------------------------------
write.table(intergenics, paste0(out_dir, "orientation/all_intergenic_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
