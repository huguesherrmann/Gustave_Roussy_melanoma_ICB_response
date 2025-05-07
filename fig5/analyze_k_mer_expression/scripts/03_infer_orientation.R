# ......................................................
# Infer intergenic region orientation
# 20/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Infer intergenic region orientation.")
parser$add_argument("--intergenics", type = "character", help = "Path to intergenic table.")
parser$add_argument("--counts", type = "character", help = "Path to intergenic strand counts.")
parser$add_argument("--out_dir", type = "character", help = "Path to output results.")
args <- parser$parse_args()

intergenics <- args$intergenics
counts <- args$counts
out_dir <- args$out_dir
# intergenics <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/orientation/all_intergenic_table.tsv"
# counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/orientation/direct_and_indirect_kmer_counts_annotated.tsv.gz"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA ----
#
# --------------------------------------------
intergenics <- read_tsv(intergenics, show_col_types = FALSE)

counts <- read_tsv(counts, show_col_types = FALSE) 


# --------------------------------------------
#
#   INFER ORIENTATION ----
#
# --------------------------------------------
orientation <- counts %>% 
   mutate(Plus = plus_strand_direct_counts + plus_strand_indirect_counts,
          Minus = minus_strand_direct_counts + minus_strand_indirect_counts) %>%
   inner_join(intergenics, ., by = "tag") %>%
   group_by(Regulon) %>%
   summarize(Plus = sum(Plus),
             Minus = sum(Minus)) %>%
   mutate(Ratio = (Plus - Minus) / (Plus + Minus)) %>%
   mutate(Strand = if_else(Ratio > 0, "+", "-"))


# --------------------------------------------
#
#   EXPORT ----
#
# --------------------------------------------
write.table(orientation, paste0(out_dir, "orientation/predicted_intergenic_region_orientation.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
