# ......................................................
# Plot TRUST4 results across different cohorts
# 04/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Plot stable regulons.")
# parser$add_argument("--regulons", type = "character", help = "Path to differential regulon table.")
# parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
# parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
# parser$add_argument("--n_top", type = "integer", help = "Number of regulons to show up in final plots.")
# args <- parser$parse_args()
# 
# regulons <- args$regulons
# correspondence <- args$correspondence
# out_dir <- args$out_dir
# n_top <- args$n_top
trust4 <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/"


# --------------------------------------------
#
#   PARSE DATA
#
# --------------------------------------------
trust4_variables <- c("BCR_reads", "TCR_reads", "Total_reads_mapped",
                      "BCR_ratio", "BCR_clonality", "Unique_BCR_CDR3", 
                      "TCR_ratio", "TCR_clonality", "Unique_TCR_CDR3",
                      "CPK_TCR", "CPK_BCR", "SHM_ratio")

stats_df <- data.frame(Sample_ID = character(),
                       Response = character(),
                       Cohort = character(),
                       BCR_reads = numeric(),
                       TCR_reads = numeric(),
                       Total_reads_mapped = numeric(),
                       BCR_ratio = numeric(),
                       TCR_ratio = numeric(),
                       CPK_TCR = numeric(),
                       CPK_BCR = numeric(),
                       SHM_ratio = numeric(),
                       BCR_clonality = numeric(),
                       TCR_clonality = numeric(),
                       Unique_TCR_CDR3 = numeric(),
                       Unique_BCR_CDR3 = numeric())

cohorts <- list.dirs(paste0(trust4), recursive = FALSE)
for (cohort in cohorts) {
  cohort_name <- str_split(cohort, 
                           "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire//")[[1]][2]
  
  tmp <- read_tsv(paste0(cohort, "/trust4/overall_stats.tsv"), show_col_types = FALSE) %>%
    select(c(Sample_ID, Response, all_of(trust4_variables))) %>%
    mutate(Cohort = cohort_name)
  
  stats_df <- stats_df %>% add_row(tmp)
}


# --------------------------------------------
#
#   PLOT ALL STATISTICS
#
# --------------------------------------------
longer_stats_df <- stats_df %>% pivot_longer(!c(Sample_ID, Response, Cohort), names_to = "Repertoire", values_to = "Value") %>%
  mutate(Cohort = factor(Cohort, levels = c("gr1234", "hugo", "gide", "riaz", "liu"))) %>%
  mutate(Repertoire = factor(Repertoire, levels = c("Total_reads_mapped", "SHM_ratio", 
                                                    "BCR_reads", "BCR_ratio", "BCR_clonality", "Unique_BCR_CDR3", "CPK_BCR",
                                                    "TCR_reads", "TCR_ratio", "TCR_clonality", "Unique_TCR_CDR3", "CPK_TCR")))

facet_labels <- c(Total_reads_mapped = "Total reads mapped",
                  BCR_reads = "BCR reads",
                  TCR_reads = "TCR reads",
                  BCR_ratio = "BCR ratio", 
                  BCR_clonality = "BCR clonality",
                  Unique_BCR_CDR3 = "Unique BCR CDR3",
                  TCR_ratio = "TCR ratio", 
                  TCR_clonality = "TCR clonality",
                  Unique_TCR_CDR3 = "Unique TCR CDR3",
                  CPK_TCR = "Clonotype per\nkilo reads (TCR)",
                  CPK_BCR = "Clonotype per\nkilo reads (BCR)",
                  SHM_ratio = "Somatic hypermutation\nrate")
all_stats_boxplot <- ggplot(longer_stats_df, aes(Cohort, Value, fill = Cohort)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = c("liu" = "#F0E442",  "gide" = "#009E73", "gr1234" = "#CC79A7", "riaz" = "burlywood3", "hugo" = "coral2"),
                    labels = c("FF total RNA\nGR1234", "FFPE polyA\nHugo", "FFPE exon capture\nGide", "FFPE polyA\nRiaz", "FFPE exon capture\nLiu")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(Repertoire ~ ., nrow = 3, scales = "free", labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        text = element_text(size = 14, color = "black"), 
        legend.position = "bottom",
        axis.ticks.x = element_blank()) +
  labs(y = "", x = "")


# --------------------------------------------
#
#   COMPARE R VS NR ACROSS ALL COHORT
#
# --------------------------------------------
response_stats_df <- longer_stats_df %>% filter(Repertoire != "Total_reads_mapped") %>%
  mutate(Facetting = paste0(Repertoire, " - ", Cohort)) %>%
  mutate(Response = if_else(Response == "responder", "R", "NR"))

facet_labels <- c(BCR_reads = "BCR reads",
                  TCR_reads = "TCR reads",
                  BCR_ratio = "BCR ratio", 
                  BCR_clonality = "BCR clonality",
                  Unique_BCR_CDR3 = "Unique BCR CDR3",
                  TCR_ratio = "TCR ratio", 
                  TCR_clonality = "TCR clonality",
                  Unique_TCR_CDR3 = "Unique TCR CDR3",
                  CPK_TCR = "Clonotype per\nkilo reads (TCR)",
                  CPK_BCR = "Clonotype per\nkilo reads (BCR)",
                  SHM_ratio = "Somatic hypermutation\nrate")
response_all_stats_boxplot <- ggplot(response_stats_df, aes(Response, Value)) +
  geom_boxplot(aes(fill = Cohort), position = position_dodge(0.8)) +
  scale_fill_manual(values = c("liu" = "#F0E442",  "gide" = "#009E73", "gr1234" = "#CC79A7", "riaz" = "burlywood3", "hugo" = "coral2"),
                    labels = c("FF total RNA\nGR1234", "FFPE polyA\nHugo", "FFPE exon capture\nGide", "FFPE polyA\nRiaz", "FFPE exon capture\nLiu")) +
  guides() +
  scale_y_continuous(trans = "log10") +
  facet_wrap(. ~ Facetting, nrow = 11, scales = "free_y") +
  stat_compare_means(comparisons = list(c("NR", "R")), method = "wilcox.test", na.rm = TRUE, label = "p.format", vjust = 1.35, tip.length = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14, color = "black"),
        text = element_text(size = 14, color = "black"),
        legend.position = "bottom") +
  labs(y = "", x = "")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(trust4, "comparison_of_all_cohorts.png"), all_stats_boxplot, height = 9, width = 10)
ggsave(paste0(trust4, "comparison_of_all_cohorts.pdf"), all_stats_boxplot, height = 9, width = 10)

ggsave(paste0(trust4, "comparison_of_response_across_all_cohorts.png"), response_all_stats_boxplot, height = 16, width = 12)
ggsave(paste0(trust4, "comparison_of_response_across_all_cohorts.pdf"), response_all_stats_boxplot, height = 16, width = 12)
