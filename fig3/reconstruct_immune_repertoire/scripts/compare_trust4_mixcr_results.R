# ......................................................
# Compare TRUST4 and MiXCR results
# 28/10/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/reconstruct_immune_repertoire/scripts/functions_for_reconstruct_immune_repertoire.R")


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
trust4 <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gr1234/trust4/overall_stats.tsv"
mixcr <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gr1234/mixcr/overall_stats.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gr1234/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
trust4 <- read_tsv(trust4, show_col_types = FALSE) %>%
   mutate(Algo = "TRUST4")

mixcr <- read_tsv(mixcr, show_col_types = FALSE) %>%
   mutate(Algo = "MiXCR")


# --------------------------------------------
#
#   FORMAT DATA
#
# --------------------------------------------
both <- inner_join(trust4, mixcr, by = "Sample_ID", suffix = c("_TRUST4", "_MiXCR"))
colnames(both) <- str_replace_all(colnames(both), "_", " ")


# --------------------------------------------
#
#   PLOT TRUST4 VS MIXCR
#
# --------------------------------------------
bcr_clonality_plot <- ggplot(both, aes(`BCR clonality TRUST4`, `BCR clonality MiXCR`)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, color = "turquoise4") +
   lims(x = c(0, 1), y = c(0, 1)) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   stat_cor(method = "spearman", label.x = 0.05, label.y = 0.9, size = 6)

bcr_diversity_plot <- ggplot(both, aes(`BCR diversity TRUST4`, `BCR diversity MiXCR`)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, color = "turquoise4") +
   lims(x = c(0, 850), y = c(0, 850)) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   stat_cor(method = "spearman", label.x = 0.05, label.y = 800, size = 6)

tcr_clonality_plot <- ggplot(both, aes(`TCR clonality TRUST4`, `TCR clonality MiXCR`)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, color = "turquoise4") +
   lims(x = c(0, 0.35), y = c(0, 0.35)) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   stat_cor(method = "spearman", label.x = 0.01, label.y = 0.3, size = 6)

tcr_diversity_plot <- ggplot(both, aes(`TCR diversity TRUST4`, `TCR diversity MiXCR`)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, color = "turquoise4") +
   lims(x = c(0, 1100), y = c(0, 1100)) +
   theme_classic() +
   theme(axis.text.x = element_text(size = 18, color = "black"),
         axis.text.y = element_text(size = 18, color = "black"),
         text = element_text(size = 18, color = "black"),
         strip.text = element_text(size = 18, color = "black")) +
   stat_cor(method = "spearman", label.x = 0.05, label.y = 1050, size = 6)

# multi_plot <- ggarrange(bcr_clonality_plot, 
#                         bcr_diversity_plot,
#                         tcr_clonality_plot,
#                         tcr_diversity_plot,
#                         ncol = 2, 
#                         nrow = 2,
#                         common.legend = TRUE)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_clonality.pdf"), bcr_clonality_plot, width = 4.5, height = 4.5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_clonality.svg"), bcr_clonality_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_clonality.png"), bcr_clonality_plot, height = 5, width = 5)

ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_diversity.pdf"), bcr_diversity_plot, width = 4.5, height = 4.5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_diversity.svg"), bcr_diversity_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/bcr_diversity.png"), bcr_diversity_plot, height = 5, width = 5)

ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_clonality.pdf"), tcr_clonality_plot, width = 4.5, height = 4.5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_clonality.svg"), tcr_clonality_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_clonality.png"), tcr_clonality_plot, height = 5, width = 5)

ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_diversity.pdf"), tcr_diversity_plot, width = 4.5, height = 4.5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_diversity.svg"), tcr_diversity_plot, height = 5, width = 5)
ggsave(paste0(out_dir, "/comparison_trust4_mixcr/tcr_diversity.png"), tcr_diversity_plot, height = 5, width = 5)