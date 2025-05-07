# ......................................................
# Get cutoff to classify samples as pervasives or ALR positives
# 08/07/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Get cutoff to classify samples as pervasives or ALR positives.")
# parser$add_argument("--design", type = "character", help = "Path to design table.")
# args <- parser$parse_args()
# 
# regulon_counts <- args$regulon_counts

design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
burdens <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/pervasive_and_ALR_burdens/pervasive_and_ALR_burdens.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA ----
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE)

burdens <- read_tsv(burdens, show_col_types = FALSE)
   

# --------------------------------------------
#
#   CALCULATE MAXIMALLY SELECTED CHI2 THRESHOLD ----
#
# --------------------------------------------
burdens <- burdens %>% inner_join(., design, by = "Sample_ID") %>%
   mutate(Response = if_else(Response == "NR", 0, 1))
pervasive_max_stat <- maxstat_test(Response ~ Pervasive_burden, data = burdens, teststat = "max")
pervasive_cutoff <- pervasive_max_stat@estimates$estimate$cutpoint

alr_max_stat <- maxstat_test(Response ~ ALR_burden, data = burdens, teststat = "max")
alr_cutoff <- alr_max_stat@estimates$estimate$cutpoint

design_burden <- burdens %>%
   mutate(Pervasive_status = if_else(Pervasive_burden >= pervasive_cutoff, 1, 0)) %>%
   mutate(ALR_status = if_else(ALR_burden >= alr_cutoff, 1, 0))


# --------------------------------------------
#
#   EXPORT ----
#
# --------------------------------------------
write.table(design_burden, paste0(out_dir, "pervasive_and_ALR_burdens/pervasive_and_ALR_status_design_gr1234.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
