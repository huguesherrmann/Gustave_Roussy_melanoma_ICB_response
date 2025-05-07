# ...................................................... 
# Import and summarize transcript-level exon and intron abundance estimates for transcript- and gene-level analysis
# 11/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_intron_expression/scripts/functions_for_summarizing_abundances_into_counts.R")


# ...................................................... 
#
#   PARSE ARGUMENTS ----
#
# ...................................................... 
parser <- ArgumentParser(description = "Summarize exon and intron transcript abundance at gene-level.")
parser$add_argument("--t2g", type = "character", help = "Path to tx2gene annotation file.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
args <- parser$parse_args()

t2g <- args$t2g
out_dir <- args$out_dir
correspondance <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/gr1234/correspondance.csv"
t2g <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/introns_index/t2g.txt"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/gr1234/a"


# ...................................................... 
#
#   LOAD DATA ----
#
# ......................................................
correspondance <- read_csv(correspondance, show_col_types = FALSE)

t2g <- read.table(t2g, header = FALSE, sep = "\t")


# ...................................................... 
#
#   GENERATE EXPRESSION TABLES ----
#
# ......................................................
# ---- Ambiguous
ambiguous_dir <- paste0(out_dir, "ambiguous_counts/")
ambiguous_counts <- summarize_abundances_into_counts(ambiguous_dir, correspondance)


# ---- Mature
mature_dir <- paste0(out_dir, "mature_counts/")
mature_counts <- summarize_abundances_into_counts(mature_dir, correspondance)


# ---- Nascent
nascent_dir <- paste0(out_dir, "nascent_counts/")
nascent_counts <- summarize_abundances_into_counts(nascent_dir, correspondance)


# txi_tpm <- tximport(samples, 
#                     tx2gene = tx2gene, 
#                     countsFromAbundance = "scaledTPM",
#                     type = "kallisto", 
#                     ignoreAfterBar = TRUE)


# ...................................................... 
#
#   EXPORT ----
#
# ......................................................
write.table(ambiguous_counts, paste0(ambiguous_dir, "ambiguous_counts.tsv"), sep = "\t", quote = FALSE)
write.table(mature_counts, paste0(mature_dir, "mature_counts.tsv"), sep = "\t", quote = FALSE)
write.table(nascent_counts, paste0(nascent_dir, "nascent_counts.tsv"), sep = "\t", quote = FALSE)