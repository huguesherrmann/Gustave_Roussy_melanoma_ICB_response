# ......................................................
# Calculate correlation between intergenic and ALR burdens
# 24/01/25
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Calculate correlation between intergenic and ALR burdens.")
parser$add_argument("--annotation", type = "character", help = "Path to regulon annotation table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--stability", type = "character", help = "Path to stable regulons.")
parser$add_argument("--n_redundancy", type = "integer", help = "Number of top regulons to integrate for calculating burden.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

annotation <- args$annotation
correspondence <- args$correspondence
stability <- args$stability
n_redundancy <- args$n_redundancy
out_dir <- args$out_dir
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
burdens <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/pervasive_and_ALR_burdens/pervasive_and_ALR_burdens.tsv"
signatures <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_published_signatures_in_cohorts/published_signatures/published_signature_scores_all_cohorts.tsv"
trust4 <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gr1234/trust4/overall_stats.tsv"
tmb <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_wes/gr1234/all_tmb/tmb.tsv"
bagaev <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_gene_expression/gr1234/Bagaev/bagaev_subtype_proba.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE) %>% select(Sample_ID, Tumor_purity)

burdens <- read_tsv(burdens, show_col_types = FALSE)

signatures <- read_tsv(signatures, show_col_types = FALSE) %>% 
   filter(Cohort == "gr1234" & Signature %in% c("Neuronal", "Keratin", "DNA_damage_response")) %>%
   select(Sample_ID, Signature, Score) %>%
   pivot_wider(id_cols = Sample_ID, names_from = Signature, values_from = Score)
      
trust4 <- read_tsv(trust4, show_col_types = FALSE) %>%
   select(Sample_ID, BCR_clonality, BCR_diversity, TCR_clonality, TCR_diversity)

tmb <- read_tsv(tmb, show_col_types = FALSE)

bagaev <- read_tsv(bagaev, show_col_types = FALSE) %>%
   rename(Sample_ID = "...1")


all_biomarkers <- design %>%
   inner_join(., burdens, by = "Sample_ID") %>% 
   inner_join(., trust4, by = "Sample_ID") %>%
   full_join(., tmb, by = c("Sample_ID" = "Sample_ID_rna")) %>%
   inner_join(., bagaev, by = "Sample_ID") %>%
   inner_join(., signatures, by = c("Sample_ID")) %>%
   select(-P2)


# --------------------------------------------
#
#   PLOT CORRELATIONS
#
# --------------------------------------------
all_biomarkers <- all_biomarkers %>% rename(`Tumor purity` = "Tumor_purity") %>%
   rename(`Intergenic burden` = "Pervasive_burden") %>%
   rename(`ALR burden` = "ALR_burden") %>%
   rename(`BCR clonality` = "BCR_clonality") %>%
   rename(`BCR diversity` = "BCR_diversity") %>%
   rename(`TCR clonality` = "TCR_clonality") %>%
   rename(`TCR diversity` = "TCR_diversity") %>%
   rename(`Immune depleted` = "D") %>%
   rename(`Fibrotic` = "F") %>%
   rename(`Immune enriched non fibrotic` = "IE/F") %>%
   rename(`Immune enriched` = "IE") %>%
   rename(`DNA damage response` = "DNA_damage_response")
   

corr_all_biomarkers <- all_biomarkers %>% 
   select(-Sample_ID, -Neuronal, -Keratin, -`DNA damage response`) %>%
   cor(., method = "spearman", use = "complete.obs")
supp_corr_all_biomarkers <- all_biomarkers %>% 
   select(-Sample_ID) %>%
   cor(., method = "spearman", use = "complete.obs")

col <- colorRampPalette(c("#3e308e", "#fdfdfe", "#610d0e"))


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
# Main plot
png(paste0(out_dir, "correlations/correlogram.png"), width = 1000, height = 1000)
corrplot(corr_all_biomarkers,
         method = "color",
         col = col(100),
         type = "lower",
         order = "AOE",
         addgrid.col = "white",
         tl.col = "black",
         tl.srt = 30,
         tl.cex = 1.3,
         diag = FALSE)
dev.off()
pdf(paste0(out_dir, "correlations/correlogram.pdf"), width = 9.5, height = 9.5)
corrplot(corr_all_biomarkers,
         method = "color",
         col = col(100),
         type = "lower",
         order = "AOE",
         addgrid.col = "white",
         tl.col = "black",
         tl.srt = 30,
         tl.cex = 1.2,
         diag = FALSE)
dev.off()
svg(paste0(out_dir, "correlations/correlogram.svg"), width = 9.5, height = 9.5)
corrplot(corr_all_biomarkers,
         method = "color",
         col = col(100),
         type = "lower",
         order = "AOE",
         addgrid.col = "white",
         tl.col = "black",
         tl.srt = 30,
         tl.cex = 1.2,
         diag = FALSE)
dev.off()

# Supp plot
pdf(paste0(out_dir, "correlations/supp_data_correlogram.pdf"), width = 9.5, height = 9.5)
corrplot(supp_corr_all_biomarkers,
         method = "color",
         col = col(100),
         type = "lower",
         order = "AOE",
         addgrid.col = "white",
         tl.col = "black",
         tl.srt = 30,
         tl.cex = 1.2,
         diag = FALSE)
dev.off()