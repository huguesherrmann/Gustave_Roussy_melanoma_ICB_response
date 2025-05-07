# ......................................................
# Plot stable regulons
# 13/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
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
# regulons <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/stability/differential_regulons.tsv"
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/regulons/all_correspondence_contigs_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/"
# n_top <- 25


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
  group_by(Feature) %>% 
  arrange(desc(N_redundancy), desc(Mean_l2fc)) %>%
  slice_head(n = 1) %>% # remove duplicated features
  arrange(desc(N_redundancy), desc(Mean_l2fc)) %>%
  head(n_top) %>%
  arrange(Mean_l2fc) %>% # Need to re-sort df with l2fc from lower to higher to plot the in the order highest to lowest
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
  arrange(desc(N_redundancy), desc(Mean_l2fc))
top_full_non_resp_df <- full_non_resp_df %>%
  group_by(Feature) %>% 
  arrange(desc(N_redundancy), desc(Mean_l2fc)) %>%
  slice_head(n = 1) %>% # remove duplicated features
  arrange(desc(N_redundancy), desc(Mean_l2fc)) %>%
  head(n_top) %>%
  arrange(desc(Mean_l2fc)) %>% # Need to re-sort df with l2fc from higher to lower to plot the in the order lowest to highest
  mutate(L2fc_stats = paste0(Mean_l2fc, " [", Low_ci_l2fc, "; ", High_ci_l2fc, "]")) %>%
  mutate(Condition = "non-responder") %>%
  mutate(Feature = factor(Feature, levels = .$Feature))


all_stable_regulons <- full_non_resp_df %>% mutate(Condition = "non_responder") %>%
  ungroup() %>%
  add_row(full_resp_df %>% mutate(Condition = "responder") %>% ungroup())


# --------------------------------------------
#
#   PLOT STABLE REGULONS
#
# --------------------------------------------
merge <- rbind(top_full_resp_df, top_full_non_resp_df)

panel_labels <- c(`non-responder` = "NR",
                  responder = "R")

l2fc_top_stable_regulons_plot <- ggplot(merge, aes(x = Mean_l2fc, y = Feature, fill = Condition)) + 
  geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = High_ci_l2fc, xmin = Low_ci_l2fc), linewidth = 0.5, height = 0.3, color = "black") +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = c("responder" = "#56B4E9", "non-responder" = "#FF9999")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.text = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"), text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 15, color = "black")) +
  facet_wrap(~ Condition, scales = "free", labeller = as_labeller(panel_labels)) +
  labs(y = "", x = "Log2 fold change (95% CI)") +
  guides(fill = "none")


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
ggsave(paste0(out_dir, "stability/l2fc_top_stable_regulons.png"), l2fc_top_stable_regulons_plot, height = 8, width = 10)
ggsave(paste0(out_dir, "stability/l2fc_top_stable_regulons.pdf"), l2fc_top_stable_regulons_plot, height = 7, width = 9)
ggsave(paste0(out_dir, "stability/l2fc_top_stable_regulons.svg"), l2fc_top_stable_regulons_plot, height = 7, width = 9)
 
write.table(all_stable_regulons, paste0(out_dir, "stability/all_stable_regulons.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
