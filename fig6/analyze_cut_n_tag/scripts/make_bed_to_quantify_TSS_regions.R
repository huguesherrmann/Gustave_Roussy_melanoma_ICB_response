# ......................................................
# Make BED file to quantify TSS expression 
# 27/09/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_cut_n_tag/scripts/functions_for_analyzing_cut_n_tag.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Make BED file of expressed and manually-selected regions.")
parser$add_argument("--metadata", type = "character", help = "Metadata table. include chromosome, start, end, strand and sample expression.")
parser$add_argument("--offset_quantif", type = "integer", default = 1000, help = "Number of nucleotides after the TSS to quantify the region.")
parser$add_argument("--mark", type = "character", help = "Mark studied. useful for the name of the output.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

metadata <- args$metadata
offset_quantif <- args$offset_quantif
mark <- args$mark
out_dir <- args$out_dir
# metadata <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/metadata/tss_pervasive_regions.tsv"
# out_dir <-"/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"
# offset_quantif <- 1000
# mark <- "H3K4me1"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
metadata <- read_tsv(metadata, show_col_types = FALSE)


# ......................................................
#
#   MAKE BED ----
#
# ......................................................
bed_quantif <- make_bed(metadata, 
                        objective = "quantif", 
                        offset_quantif = offset_quantif)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(bed_quantif, 
            paste0(out_dir, "bed/", mark, "/bed_for_quantification_of_tss_regions.bed.saf"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
