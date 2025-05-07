# ......................................................
# Plot TRUST4 results
# 29/05/24
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
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/reconstruct_immune_repertoire/gr1234/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
trust4 <- read_tsv(trust4, show_col_types = FALSE)

design <- read_tsv(design, show_col_types = FALSE)


# --------------------------------------------
#
#   FORMAT DATA
#
# --------------------------------------------
longer_trust4 <- trust4 %>% 
  select(Sample_ID, Response, TCR_clonality, BCR_clonality, TCR_diversity, BCR_diversity) %>%
  pivot_longer(!c(Sample_ID, Response), names_to = "Stats", values_to = "Score") %>% 
  mutate(Response = if_else(Response == "responder", "R", "NR"))


# --------------------------------------------
#
#   ALL IN ONE PLOT
#
# --------------------------------------------
# trust4_for_plot <- trust4 %>% select(Sample_ID, all_of(trust4_variables)) %>% 
#   column_to_rownames("Sample_ID") %>%
#   t() %>% as.data.frame() %>%
#   rownames_to_column("Variable")
# 
# trust4_plot <- make_violin_plot_with_pvals_and_facets(trust4_for_plot, design)


# --------------------------------------------
#
#   INDIVIDUAL PLOTS
#
# --------------------------------------------
bcr_clonality <-  ggplot(longer_trust4 %>% filter(Stats == "BCR_clonality"), aes(Response, Score, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.15, vjust = 0.5, label = "p.format", size = 6.5) +
  expand_limits(y = c(0.3, 1)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "BCR clonality", x = "") +
  guides(fill = "none")

tcr_clonality <- ggplot(longer_trust4 %>% filter(Stats == "TCR_clonality"), aes(Response, Score, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.1, vjust = 0.5, label = "p.format", size = 6.5) +
  expand_limits(y = c(0.04, 0.4)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "TCR clonality", x = "") +
  guides(fill = "none")

bcr_diversity <- ggplot(longer_trust4 %>% filter(Stats == "BCR_diversity"), aes(Response, Score, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.2, vjust = 0.5, label = "p.format", size = 6.5) +
  expand_limits(y = c(30, 1300)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "BCR diversity", x = "") +
  guides(fill = "none")

tcr_diversity <- ggplot(longer_trust4 %>% filter(Stats == "TCR_diversity"), aes(Response, Score, fill = Response)) +
  geom_boxplot() +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label.x = 1.2, vjust = 0.5, label = "p.format", size = 6.5) +
  expand_limits(y = c(30, 1300)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "bottom") +
  labs(y = "TCR diversity", x = "") +
  guides(fill = "none")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "plots/bcr_clonality.png"), bcr_clonality, width = 3, height = 4)
ggsave(paste0(out_dir, "plots/bcr_clonality.pdf"), bcr_clonality, width = 3, height = 3.5)
ggsave(paste0(out_dir, "plots/bcr_clonality.svg"), bcr_clonality, width = 3, height = 3.5)

ggsave(paste0(out_dir, "plots/tcr_clonality.png"), tcr_clonality, width = 3, height = 4)
ggsave(paste0(out_dir, "plots/tcr_clonality.pdf"), tcr_clonality, width = 3, height = 3.5)
ggsave(paste0(out_dir, "plots/tcr_clonality.svg"), tcr_clonality, width = 3, height = 3.5)

ggsave(paste0(out_dir, "plots/bcr_diversity.png"), bcr_diversity, width = 3, height = 4)
ggsave(paste0(out_dir, "plots/bcr_diversity.pdf"), bcr_diversity, width = 3, height = 3.5)
ggsave(paste0(out_dir, "plots/bcr_diversity.svg"), bcr_diversity, width = 3, height = 3.5)

ggsave(paste0(out_dir, "plots/tcr_diversity.png"), tcr_diversity, width = 3, height = 4)
ggsave(paste0(out_dir, "plots/tcr_diversity.pdf"), tcr_diversity, width = 3, height = 3.5)
ggsave(paste0(out_dir, "plots/tcr_diversity.svg"), tcr_diversity, width = 3, height = 3.5)
