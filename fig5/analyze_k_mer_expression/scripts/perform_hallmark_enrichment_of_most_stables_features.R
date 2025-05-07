# ......................................................
# Perform enrichment analysis with stable features
# 13/05/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "PPerform enrichment analysis with stable features.")
parser$add_argument("--stability", type = "character", help = "Path to stable regulons.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs and stability directory.")
parser$add_argument("--alpha", type = "double", default = 0.05, help = "Type-I error threshold for enrichment test.")
args <- parser$parse_args()

stability <- args$stability
out_dir <- args$out_dir
alpha <- args$alpha
# stability <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/stability/all_stable_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234_w_purity/"
# alpha <- 0.05


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
stables <- read_tsv(stability, show_col_types = FALSE) %>%
  filter(N_redundancy >= 1)

hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  select(gs_name, gene_symbol)


# --------------------------------------------
#
#   PERFORM ENRICHMENT ANALYSIS
#
# --------------------------------------------
responder_stable_regulons <- stables %>% 
  filter(Condition == "responder") %>%
  filter(Class == "exon" | Class == "intron") %>%
  distinct(Feature) %>%
  pull(Feature)
responder_hallmarks <- enricher(responder_stable_regulons, TERM2GENE = hallmarks, pvalueCutoff = alpha) %>%
  as.data.frame() %>%
  mutate(Condition = "R")

non_responder_stable_regulons <- stables %>% 
  filter(Condition == "non_responder") %>%
  filter(Class == "exon" | Class == "intron") %>%
  distinct(Feature) %>%
  pull(Feature)
non_responder_hallmarks <- enricher(non_responder_stable_regulons, TERM2GENE = hallmarks, pvalueCutoff = alpha) %>%
  as.data.frame() %>%
  mutate(Condition = "NR")

k_mer_hallmarks <- responder_hallmarks %>% add_row(non_responder_hallmarks) %>% 
  mutate(Cohort = "GR\nk_mers") %>%
  mutate(ID = str_remove(ID, "HALLMARK_")) %>%
  mutate(NES = 3)


# --------------------------------------------
#
#   PLOT RESULTS
#
# --------------------------------------------
hallmark_dotplot <- ggplot(k_mer_hallmarks, aes(Cohort, ID, fill = Condition, size = NES)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(-1, 6), breaks = seq(1, 6, 1)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c("R" = "#56B4E9", "NR" = "#FF9999")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(size = "none") +
  labs(fill = "Association:", x = NULL, y = NULL)


# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
if (nrow(k_mer_hallmarks) >= 1) {
  ggsave(paste0(out_dir, "hallmarks_of_diff_regulons/hallmark_dotplot_of_exonic_intronic_differential_features.png"), 
         hallmark_dotplot, height = 6.5, width = 8)
  ggsave(paste0(out_dir, "hallmarks_of_diff_regulons/hallmark_dotplot_of_exonic_intronic_differential_features.pdf"), 
         hallmark_dotplot, height = 6.5, width = 8)
  ggsave(paste0(out_dir, "hallmarks_of_diff_regulons/hallmark_dotplot_of_exonic_intronic_differential_features.svg"), 
         hallmark_dotplot, height = 6.5, width = 8)
}

write.table(k_mer_hallmarks, paste0(out_dir, "hallmarks_of_diff_regulons/hallmarks_of_exonic_intronic_differential_features.tsv"), sep = "\t", row.names = FALSE)
