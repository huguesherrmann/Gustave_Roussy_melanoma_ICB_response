# ......................................................
# Output intergenic contigs in a format chromome tag genomic_position
# 30/03/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")
options(scipen = 999)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Output intergenic contigs in a format chromome tag genomic_position.")
parser$add_argument("--annotation", type = "character", help = "Path to contig annotation table.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
args <- parser$parse_args()

annotation <- args$annotation
out_dir <- args$out_dir
#annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/mask_k_mer_matrix/gr123/annotation_merged_contig_sorted_gide_masked_kmers_6_5_kmers.tsv.gz"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
annotation <- fread(annotation)


# ......................................................
#
#   GET GENOMIC POSITION ----
#
# ......................................................
hg38_offsets <- get_chromosome_offset()

intergenics <- annotation[is_exonic == 0 & is_intronic == 0 & is.na(gene_symbol) & is.na(repeats), ] %>%
  mutate(Center = start + round((end - start) / 2)) %>%
  select(chromosome, tag, Center) %>%
  arrange(chromosome, Center)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(intergenics, paste0(out_dir, "subgraphs/intergenics_to_connect.txt"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

