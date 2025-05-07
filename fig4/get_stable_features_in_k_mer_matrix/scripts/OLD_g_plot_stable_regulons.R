# ......................................................
# Plot stable regulons
# 13/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(aplot))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_processing_regulons.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Plot stable regulons.")
parser$add_argument("--regulons", type = "character", help = "Path to differential regulon table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
parser$add_argument("--n_top", type = "integer", help = "Number of regulons to show up in final plots.")
args <- parser$parse_args()

regulons <- args$regulons
correspondence <- args$correspondence
out_dir <- args$out_dir
n_top <- args$n_top
# regulons <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_cv_20_5/stability/differential_regulons.tsv"
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_cv_20_5/regulons/all_correspondence_contigs_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/test/stability/"
# n_top <- 17


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
stable_regulons <- fread(regulons)

correspondence <- fread(correspondence)


# --------------------------------------------
#
#   GET STABLE REGULONS
#
# --------------------------------------------
arbitrary_top_stables <- round(max(stable_regulons$Id) / 2)

resp_df <- stable_regulons %>% filter(Condition == "responder") %>%
  group_by(Regulon) %>%
  summarize(N_redundancy = n(),
            Mean_l2fc = round(mean(log2FoldChange), 2),
            Low_ci_l2fc = round(quantile(log2FoldChange, 0.025), 2),
            High_ci_l2fc = round(quantile(log2FoldChange, 0.975), 2)) %>%
  arrange(desc(N_redundancy), desc(Mean_l2fc)) %>%
  setDT()
full_resp_df <- merge(resp_df, correspondence %>% select(-tag), by = "Regulon") %>%
  group_by(Regulon) %>%
  slice(1) %>%
  filter(N_redundancy >= arbitrary_top_stables) %>%
  arrange(desc(N_redundancy), desc(Mean_l2fc))
top_full_resp_df <- full_resp_df %>%
  head(n_top) %>%
  arrange(N_redundancy, Mean_l2fc) %>% # Need to re-sort df with l2fc from lower to higher to plot the in the order highest to lowest
  mutate(L2fc_stats = paste0(Mean_l2fc, " [", Low_ci_l2fc, "; ", High_ci_l2fc, "]")) %>%
  mutate(Condition = "responder") %>%
  mutate(Feature = factor(Feature, levels = .$Feature)) 
  
non_resp_df <- stable_regulons %>% filter(Condition == "non_responder") %>%
  group_by(Regulon) %>%
  summarize(N_redundancy = n(),
            Mean_l2fc = round(mean(log2FoldChange), 2),
            Low_ci_l2fc = round(quantile(log2FoldChange, 0.975), 2),
            High_ci_l2fc = round(quantile(log2FoldChange, 0.025), 2)) %>%
  filter(N_redundancy >= arbitrary_top_stables) %>%
  arrange(desc(N_redundancy), Mean_l2fc) %>%
  setDT()
full_non_resp_df <- merge(non_resp_df, correspondence %>% select(-tag), by = "Regulon") %>%
  group_by(Regulon) %>%
  slice(1) %>%
  filter(N_redundancy >= arbitrary_top_stables) %>%
  arrange(desc(N_redundancy), Mean_l2fc)
top_full_non_resp_df <- full_non_resp_df %>%
  head(n_top) %>%
  arrange(N_redundancy, desc(Mean_l2fc)) %>% # Need to re-sort df with l2fc from higher to lower to plot the in the order lowest to highest
  mutate(L2fc_stats = paste0(Mean_l2fc, " [", Low_ci_l2fc, "; ", High_ci_l2fc, "]")) %>%
  mutate(Condition = "non-responder") %>%
  mutate(Feature = factor(Feature, levels = .$Feature))


# --------------------------------------------
#
#   PLOT STABLE REGULONS
#
# --------------------------------------------
class_resp_plot <- ggplot(top_full_resp_df, aes(x = 0, y = Feature, label = Class)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(color = "black", size = 5) +
  scale_fill_identity(guide = "none") +
  theme_void() #+
  #coord_fixed(ratio = 0.2)
stable_resp_plot <- ggplot(top_full_resp_df, aes(x = N_redundancy, y = Feature)) +
  geom_segment(aes(x = 0, y = Feature, xend = N_redundancy, yend = Feature), color = "grey40") +
  geom_point(size = 2.5, color = "#00BFC4") +
  facet_wrap(~ Condition) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 14, color = "black")) + 
  labs(x = "Stable features", y = "") #+
  #coord_fixed(ratio = 5)
l2fc_resp_plot <- ggplot(top_full_resp_df, aes(x = 0, y = Feature, label = L2fc_stats)) +
  geom_tile(fill = "white", color = "white") +
  geom_text(color = "black", size = 4) +
  scale_fill_identity(guide = "none") +
  labs(x = "Log 2 fold change") +
  theme_void()

final_resp_plot <- class_resp_plot %>% insert_right(stable_resp_plot) %>% insert_right(l2fc_resp_plot)


class_non_resp_plot <- ggplot(top_full_non_resp_df, aes(x = 0, y = Feature, label = Class)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(color = "black", size = 5) +
  scale_fill_identity(guide = "none") +
  theme_void() #+
#coord_fixed(ratio = 0.2)
stable_non_resp_plot <- ggplot(top_full_non_resp_df, aes(x = N_redundancy, y = Feature)) +
  geom_segment(aes(x = 0, y = Feature, xend = N_redundancy, yend = Feature), color = "grey40") +
  geom_point(size = 2.5, color = "#F8766D") +
  facet_wrap(~ Condition) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 14, color = "black")) + 
  labs(x = "Stable features", y = "") #+
#coord_fixed(ratio = 5)
l2fc_non_resp_plot <- ggplot(top_full_non_resp_df, aes(x = 0, y = Feature, label = L2fc_stats)) +
  geom_tile(fill = "white", color = "white") +
  geom_text(color = "black", size = 4) +
  scale_fill_identity(guide = "none") +
  labs(x = "Log 2 fold change") +
  theme_void()

final_non_resp_plot <- class_non_resp_plot %>% insert_right(stable_non_resp_plot) %>% insert_right(l2fc_non_resp_plot)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "stability/responder_stable_regulons.png"), final_resp_plot, height = 6, width = 9)
ggsave(paste0(out_dir, "stability/non_responder_stable_regulons.png"), final_non_resp_plot, height = 6, width = 9)

write.table(full_resp_df, paste0(out_dir, "stability/responder_stable_regulons.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
write.table(full_non_resp_df, paste0(out_dir, "stability/non_responder_stable_regulons.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)