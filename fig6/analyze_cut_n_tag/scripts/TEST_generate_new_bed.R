# ......................................................
# Generate a bed file of x bases around the center of the peaks called by macs2
# 18/12/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Generate a bed file of x bases around the center of the peaks called by macs2.")
parser$add_argument("--bed", type = "character", help = "Path to the peak bed file.")
parser$add_argument("--n_bases", type = "integer", help = "Number of bases around peaks.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

bed <- args$bed
n_bases <- args$n_bases
out_dir <- args$out_dir
# bed <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/macs2/D1332T46_H3K4me1.peaks.bed"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/test/"
# n_bases <- 500


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
bed <- read_tsv(bed, show_col_types = FALSE, col_names = FALSE)
colnames(bed) <- c("Chromosome", "Start", "End")


# ......................................................
#
#   GENERATE BED ----
#
# ......................................................
upstream_bed <- bed %>% mutate(End = Start,
                                 Start = Start - n_bases) %>%
   mutate(Region = paste0(Chromosome, ":", Start, "-", End)) %>%
   mutate(Tag = paste0(row_number(), "_1"))

downstream_bed <- bed %>% mutate(Start = End, 
                                 End = End + n_bases) %>%
   mutate(Region = paste0(Chromosome, ":", Start, "-", End)) %>%
   mutate(Tag = paste0(row_number(), "_2"))

full_bed <- upstream_bed %>% add_row(downstream_bed) %>%
   arrange(Tag) %>%
   select(Region)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(full_bed, paste0(out_dir, "bed_upstream_downstream.bed"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
