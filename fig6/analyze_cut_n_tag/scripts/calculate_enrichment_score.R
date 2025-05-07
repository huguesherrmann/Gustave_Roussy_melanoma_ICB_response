# ......................................................
# Calculate an enrichment score for each region in function of the expression of the regions and then compare the distribution of the ES
# 11/10/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Make BED file of expressed and manually-selected regions.")
parser$add_argument("--metadata", type = "character", help = "Metadata table. include chromosome, start, end, strand and sample expression.")
parser$add_argument("--mark", type = "character", help = "Mark studied. useful for the name of the output.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

metadata <- args$metadata
mark <- args$mark
out_dir <- args$out_dir
# feature_counts_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/featureCounts/"
# bam_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/"
# out_dir <-"/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"
# mark <- "H3K4me1"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
files <- list.files(paste0(bam_dir, mark), recursive = FALSE, full.names = FALSE, pattern = ".bai")[-4] # remove D437T40


# ......................................................
#
#   CALCULATE ENRICHMENT SCORES ----
#
# ......................................................
all_es <- data.frame(Sample_ID = character(),
                     Type = character(),
                     Enrichment_score = numeric())
for (file in files) {
   sample_ID <- str_split(file, "_")[[1]][1]
   
   signal_counts <- read_tsv(paste0(feature_counts_dir, mark, "/", sample_ID, "_expressed_and_TSS_regions_signal.txt"), show_col_types = FALSE, skip = 1)
   colnames(signal_counts)[7] <- "Signal_counts"
   signal_counts <- signal_counts %>% mutate(Signal_counts = if_else(Signal_counts == 0, 1, Signal_counts)) # To prevent NA in calculating ratios
   noise_counts <- read_tsv(paste0(feature_counts_dir, mark, "/", sample_ID, "_expressed_and_TSS_regions_noise.txt"), show_col_types = FALSE, skip = 1)
   colnames(noise_counts)[7] <- "Noise_counts"
   noise_counts <- noise_counts %>% mutate(Noise_counts = if_else(Noise_counts == 0, 1, Noise_counts)) # To prevent NA in calculating ratios
   
   expression_es <- cbind(signal_counts %>% select(Signal_counts), noise_counts %>% select(Noise_counts)) %>%
      mutate(Enrichment_score = Signal_counts / Noise_counts) %>%
      mutate(Sample_ID = sample_ID) %>%
      mutate(Type = "expression") %>%
      select(Sample_ID, Type, Enrichment_score)
   
   
   signal_counts <- read_tsv(paste0(feature_counts_dir, mark, "/", sample_ID, "_not_expressed_and_TSS_regions_signal.txt"), show_col_types = FALSE, skip = 1)
   colnames(signal_counts)[7] <- "Signal_counts"
   signal_counts <- signal_counts %>% mutate(Signal_counts = if_else(Signal_counts == 0, 1, Signal_counts)) # To prevent NA in calculating ratios
   noise_counts <- read_tsv(paste0(feature_counts_dir, mark, "/", sample_ID, "_not_expressed_and_TSS_regions_noise.txt"), show_col_types = FALSE, skip = 1)
   colnames(noise_counts)[7] <- "Noise_counts"
   noise_counts <- noise_counts %>% mutate(Noise_counts = if_else(Noise_counts == 0, 1, Noise_counts)) # To prevent NA in calculating ratios
   
   no_expression_es <- cbind(signal_counts %>% select(Signal_counts), noise_counts %>% select(Noise_counts)) %>%
      mutate(Enrichment_score = Signal_counts / Noise_counts) %>%
      mutate(Sample_ID = sample_ID) %>%
      mutate(Type = "no_expression") %>%
      select(Sample_ID, Type, Enrichment_score)
   
   both <- expression_es %>% add_row(no_expression_es)
   
   all_es <- all_es %>% add_row(both)
}

t.test(Enrichment_score ~ Type, data = all_es)
wilcox.test(Enrichment_score ~ Type, data = all_es)
